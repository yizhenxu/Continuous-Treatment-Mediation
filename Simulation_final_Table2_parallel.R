# Parallelized RStudio/standard-R version of Simulation_final_Table2.R
#
# How to run:
#   1. Open this file in R/RStudio.
#   2. Edit the USER SETTINGS block below.
#   3. Click Source, or run source("Simulation_final_Table2_parallel.R").

################################################################################
############################ USER SETTINGS ######################################
################################################################################
path = "~/JCI continuous treatment mediation/Continuous-Treatment-Mediation-main"
setwd(path)

n_runs <- 1000L
folds <- 3L
sample_sizes <- c(2000L, 5000L, 8000L)
bandwidth_vec <- seq(0.1, 1.0, 0.1)
out_dir <- "table2_parallel_results"

# Use one less than the number of physical cores by default.
n_cores <- 19#max(1L, parallel::detectCores(logical = FALSE) - 1L)

# Monte Carlo sample size used only to approximate the true value.
true_mc_n <- 1000000L
rdnum <- 2L
master_seed <- 0L

################################################################################
############################ PACKAGES ###########################################
################################################################################

required_packages <- c("MASS", "foreach", "doParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Please install missing packages first: install.packages(c(", paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(MASS)
  library(foreach)
  library(doParallel)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
############################ GLOBAL PARAMETERS ##################################
################################################################################

set.seed(master_seed)

a <- 4.5
a_prime <- 6
mean_X <- c(0, 0, 0, 0)
sigma_X <- diag(c(0.25, 0.1, 0.8, 0.5))
y1 <- -1
y2 <- 5

sigmoid <- function(x) {
  p <- plogis(x)
  pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
}

simulate <- function(n) {
  X <- MASS::mvrnorm(n, mean_X, sigma_X)
  A <- 5 + X[, 1] + 0.2 * X[, 1]^2 + rnorm(n)
  p_M <- sigmoid(-5 + 5 * A + 2 * X[, 2] + 10 * A * X[, 3])
  M <- rbinom(n, size = 1L, prob = p_M)
  Y <- y1 * A + 20 * M + y2 * M * X[, 1] + X[, 2] + rnorm(n)
  data.frame(Y = Y, A = A, M = M, X)
}

compute_true_value <- function(n_mc = 1000000L) {
  x <- MASS::mvrnorm(n_mc, mean_X, sigma_X)
  logistic_M <- sigmoid(-5 + 5 * a_prime + 2 * x[, 2] + 10 * a_prime * x[, 3])
  result_m1 <- (y1 * a + 20 + y2 * x[, 1] + x[, 2]) * logistic_M
  result_m0 <- (y1 * a + x[, 2]) * (1 - logistic_M)
  mean(result_m1 + result_m0)
}

true_val_mc <- compute_true_value(true_mc_n)

################################################################################
############################ INFERENCE CALCULATORS ##############################
################################################################################

outcome_regression <- function(df_inference, OR_model, eval_a) {
  df_a <- df_inference
  df_a$A <- eval_a
  predict(OR_model, newdata = df_a, type = "response")
}

propensity_score <- function(df_inference, propensity_model, eval_a) {
  mu <- predict(propensity_model, newdata = df_inference, type = "response")
  sigma <- summary(propensity_model)$sigma
  ps <- dnorm(eval_a, mean = mu, sd = sigma)
  1 / pmax(ps, .Machine$double.xmin)
}

density_ratio <- function(df_inference, mediator_model, a, a_prime) {
  df_a_prime <- df_inference; df_a_prime$A <- a_prime
  df_a <- df_inference; df_a$A <- a
  
  p_a <- predict(mediator_model, newdata = df_a, type = "response")
  p_ap <- predict(mediator_model, newdata = df_a_prime, type = "response")
  p_a <- pmin(pmax(p_a, .Machine$double.eps), 1 - .Machine$double.eps)
  p_ap <- pmin(pmax(p_ap, .Machine$double.eps), 1 - .Machine$double.eps)
  
  ifelse(df_inference$M == 1L, p_ap / p_a, (1 - p_ap) / (1 - p_a))
}

eta_hat_fn <- function(df_inference, mediator_model, OR_model, a, a_prime) {
  df_a_prime <- df_inference; df_a_prime$A <- a_prime
  p_ap <- predict(mediator_model, newdata = df_a_prime, type = "response")
  p_ap <- pmin(pmax(p_ap, .Machine$double.eps), 1 - .Machine$double.eps)
  
  df_y <- df_inference; df_y$A <- a
  df_y$M <- 1L
  y_m1 <- predict(OR_model, newdata = df_y, type = "response")
  df_y$M <- 0L
  y_m0 <- predict(OR_model, newdata = df_y, type = "response")
  
  y_m1 * p_ap + y_m0 * (1 - p_ap)
}

################################################################################
############################ ESTIMATORS #########################################
################################################################################

estimate_psi_yma <- function(df_inference, df_nuisance, bandwidth) {
  OR_model <- glm(Y ~ A + M + M:X1 + X2 - 1, family = gaussian, data = df_nuisance)
  propensity_model <- lm(A ~ X1 + I(X1^2), data = df_nuisance)
  mediator_model <- glm(M ~ A + X2 + A:X3, data = df_nuisance, family = binomial)
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  ps_a <- propensity_score(df_inference, propensity_model, eval_a = a)
  ps_ap <- propensity_score(df_inference, propensity_model, eval_a = a_prime)
  dens_ratio <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  k1 <- dnorm((df_inference$A - a) / bandwidth) / bandwidth
  k2 <- dnorm((df_inference$A - a_prime) / bandwidth) / bandwidth
  
  term1 <- k1 * ps_a * dens_ratio * (df_inference$Y - OR_prediction)
  term2 <- k2 * ps_ap * (OR_prediction - eta)
  
  influence_val <- term1 + term2 + eta
  c(psi = mean(influence_val), var = mean((influence_val - mean(influence_val))^2))
}

crossfit <- function(df, folds, bw) {
  n_rows <- nrow(df)
  row_partitions <- split(seq_len(n_rows), rep(seq_len(folds), length.out = n_rows))
  
  psi_estimate_cf <- matrix(NA_real_, nrow = folds, ncol = 2L)
  
  for (i in seq_len(folds)) {
    partition <- row_partitions[[i]]
    psi_estimate_cf[i, ] <- estimate_psi_yma(df[partition, , drop=FALSE], df[-partition, , drop=FALSE], bw)
  }
  
  psi_est <- mean(psi_estimate_cf[, 1])
  var_est <- mean(psi_estimate_cf[, 2])
  
  lb <- psi_est - qnorm(0.975) * sqrt(var_est / n_rows)
  ub <- psi_est + qnorm(0.975) * sqrt(var_est / n_rows)
  cp <- as.integer(true_val_mc > lb & true_val_mc < ub)
  
  c(psi_est, var_est, cp)
}

run_one_rep <- function(i, sample_size, folds, bandwidth_vec, rep_seeds) {
  set.seed(rep_seeds[i])
  df <- simulate(sample_size)
  
  # Return a flattened vector of length 3 * length(bandwidth_vec)
  res <- lapply(bandwidth_vec, function(bw) crossfit(df, folds, bw))
  unlist(res)
}

################################################################################
############################ PARALLEL DRIVER ####################################
################################################################################

run_sample_size_parallel <- function(sample_size, n_runs, folds, n_cores, bandwidth_vec, out_dir) {
  message(sprintf("Starting sample_size = %s, n_runs = %s, workers = %s", sample_size, n_runs, n_cores))
  
  set.seed(master_seed + sample_size)
  rep_seeds <- sample.int(.Machine$integer.max, n_runs)
  
  if (n_cores <= 1L) {
    res_list <- lapply(seq_len(n_runs), run_one_rep, sample_size = sample_size, folds = folds, bandwidth_vec = bandwidth_vec, rep_seeds = rep_seeds)
    res_mat <- do.call(rbind, res_list)
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() }, add = TRUE)
    
    res_mat <- foreach::foreach(
      i = seq_len(n_runs),
      .combine = rbind,
      .inorder = TRUE,
      .packages = "MASS",
      .export = c("a", "a_prime", "mean_X", "sigma_X", "y1", "y2", "true_val_mc",
                  "sigmoid", "simulate", "outcome_regression", "propensity_score",
                  "density_ratio", "eta_hat_fn", "estimate_psi_yma", "crossfit", 
                  "run_one_rep", "rep_seeds", "bandwidth_vec")
    ) %dopar% {
      run_one_rep(i, sample_size, folds, bandwidth_vec, rep_seeds)
    }
  }
  
  # Reshape back into 3D Array: (n_runs, 3 metrics, length(bandwidth_vec))
  psi_hat_list <- array(NA_real_, dim = c(n_runs, 3, length(bandwidth_vec)))
  for(r in seq_len(n_runs)) {
    psi_hat_list[r, , ] <- matrix(res_mat[r, ], nrow = 3, ncol = length(bandwidth_vec))
  }
  
  save_path <- file.path(out_dir, paste0("table2_", sample_size, ".RData"))
  save(psi_hat_list, file = save_path)
  message(sprintf("Saved %s", normalizePath(save_path, mustWork = FALSE)))
  
  psi_hat_list
}

make_summary_table <- function(results_by_n, true_val, bandwidth_vec, rdnum) {
  # Initialize empty lists to store rows for each metric
  bias_rows <- list()
  var_rows <- list()
  cov_rows <- list()
  
  sample_sizes <- names(results_by_n)
  
  # Calculate the metrics and assign them to the respective metric lists
  for (sample_size in sample_sizes) {
    mat <- results_by_n[[sample_size]]
    
    # 1. Absolute average bias (RMSE)
    bias_rmse <- paste0(
      round(abs(apply(mat[, 1, ] - true_val, 2, mean)), rdnum), " (",
      round(sqrt(apply(mat[, 1, ] - true_val, 2, function(x) mean(x^2))), rdnum), ")"
    )
    
    # 2. Mean (SD) of sqrt(Vhat)
    ave_var <- paste0(
      round(apply(sqrt(mat[, 2, ]), 2, mean), rdnum), " (",
      round(apply(sqrt(mat[, 2, ]), 2, sd), rdnum), ")"
    )
    
    # 3. Coverage
    coverage <- round(apply(mat[, 3, ], 2, mean), rdnum)
    
    bias_rows[[sample_size]] <- bias_rmse
    var_rows[[sample_size]] <- ave_var
    cov_rows[[sample_size]] <- coverage
  }
  
  # Combine lists into matrices
  bias_mat <- do.call(rbind, bias_rows)
  var_mat <- do.call(rbind, var_rows)
  cov_mat <- do.call(rbind, cov_rows)
  
  # Create the Metric column layout (Label on the first row, empty below)
  empty_pad <- rep("", length(sample_sizes) - 1)
  bias_metric_col <- c("Absolute Average Bias (RMSE)", empty_pad)
  var_metric_col <- c("Mean (SD) of sqrt(Vhat)", empty_pad)
  cov_metric_col <- c("Coverage", empty_pad)
  
  # Construct data blocks for each metric section
  bias_block <- data.frame(Metric = bias_metric_col, n = sample_sizes, bias_mat, stringsAsFactors = FALSE)
  var_block <- data.frame(Metric = var_metric_col, n = sample_sizes, var_mat, stringsAsFactors = FALSE)
  cov_block <- data.frame(Metric = cov_metric_col, n = sample_sizes, cov_mat, stringsAsFactors = FALSE)
  
  # Bind all blocks vertically
  final_table <- rbind(bias_block, var_block, cov_block)
  
  # Rename columns cleanly: Metric, n, and raw bandwidth values
  colnames(final_table)[3:ncol(final_table)] <- bandwidth_vec
  
  # Remove row names for cleaner print output
  row.names(final_table) <- NULL
  
  return(final_table)
}
################################################################################
############################ RUN ANALYSIS #######################################
################################################################################

results_by_n <- list()
for (sample_size in sample_sizes) {
  results_by_n[[as.character(sample_size)]] <- run_sample_size_parallel(
    sample_size = sample_size, n_runs = n_runs, folds = folds, 
    n_cores = n_cores, bandwidth_vec = bandwidth_vec, out_dir = out_dir
  )
}

summary_table <- make_summary_table(results_by_n, true_val_mc, bandwidth_vec, rdnum)
summary_path <- file.path(out_dir, "table2_summary_parallel.csv")
write.csv(summary_table, summary_path, row.names = FALSE)

print(summary_table, row.names = FALSE)
message(sprintf("True value used: %.10f", true_val_mc))
message(sprintf("Saved summary table: %s", normalizePath(summary_path, mustWork = FALSE)))