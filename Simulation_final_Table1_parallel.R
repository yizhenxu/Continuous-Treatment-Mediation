# Parallelized RStudio/standard-R version of Simulation_final_Table1.R
#
# How to run:
#   1. Open this file in R/RStudio.
#   2. Edit the USER SETTINGS block below.
#   3. Click Source, or run source("Simulation_final_Table1_parallel_RStudio.R").
#
# Required packages:
#   install.packages(c("MASS", "foreach", "doParallel"))
#
# Main changes relative to the original script:
#   - Uses foreach + doParallel instead of terminal command-line arguments.
#   - Parallelizes Monte Carlo simulation replicates within each sample size.
#   - Refactors duplicated estimator code into a single estimate_psi() function.
#   - Saves per-sample-size .RData files and a combined CSV summary table.

################################################################################
############################ USER SETTINGS ######################################
################################################################################
path = "~/JCI continuous treatment mediation/Continuous-Treatment-Mediation-main"
setwd(path)

n_runs <- 1000L
folds <- 3L
sample_sizes <- c(2000L, 5000L, 8000L)
out_dir <- "table1_parallel_results"

# Use one less than the number of physical cores by default.
# On shared servers, manually set this lower, e.g. n_cores <- 8L.
n_cores <- 19#max(1L, parallel::detectCores(logical = FALSE) - 1L)

# Monte Carlo sample size used only to approximate the true value.
true_mc_n <- 1000000L

# Rounding used in the final summary CSV.
rdnum <- 2L

# Reproducibility seed for true-value approximation and the replicate seed list.
master_seed <- 0L

################################################################################
############################ PACKAGES ###########################################
################################################################################

required_packages <- c("MASS", "foreach", "doParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Please install missing packages first: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
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

compute_true_value <- function(n_mc = 100000L) {
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
  df_a_prime <- df_inference
  df_a_prime$A <- a_prime
  df_a <- df_inference
  df_a$A <- a
  
  p_a <- predict(mediator_model, newdata = df_a, type = "response")
  p_ap <- predict(mediator_model, newdata = df_a_prime, type = "response")
  p_a <- pmin(pmax(p_a, .Machine$double.eps), 1 - .Machine$double.eps)
  p_ap <- pmin(pmax(p_ap, .Machine$double.eps), 1 - .Machine$double.eps)
  
  ifelse(df_inference$M == 1L, p_ap / p_a, (1 - p_ap) / (1 - p_a))
}

eta_hat_fn <- function(df_inference, mediator_model, OR_model, a, a_prime) {
  df_a_prime <- df_inference
  df_a_prime$A <- a_prime
  p_ap <- predict(mediator_model, newdata = df_a_prime, type = "response")
  p_ap <- pmin(pmax(p_ap, .Machine$double.eps), 1 - .Machine$double.eps)
  
  df_y <- df_inference
  df_y$A <- a
  df_y$M <- 1L
  y_m1 <- predict(OR_model, newdata = df_y, type = "response")
  df_y$M <- 0L
  y_m0 <- predict(OR_model, newdata = df_y, type = "response")
  
  y_m1 * p_ap + y_m0 * (1 - p_ap)
}

################################################################################
############################ MODEL SPECIFICATIONS ###############################
################################################################################

model_specs <- list(
  YMA = list(
    y_formula = Y ~ A + M + M:X1 + X2 - 1,
    a_formula = A ~ X1 + I(X1^2),
    m_formula = M ~ A + X2 + A:X3
  ),
  YM = list(
    y_formula = Y ~ A + M + M:X1 + X2 - 1,
    a_formula = A ~ 1,
    m_formula = M ~ A + X2 + A:X3
  ),
  MA = list(
    y_formula = Y ~ M + X1 + X2,
    a_formula = A ~ X1 + I(X1^2),
    m_formula = M ~ A + X2 + A:X3
  ),
  YA = list(
    y_formula = Y ~ A + M + M:X1 + X2 - 1,
    a_formula = A ~ X1 + I(X1^2),
    m_formula = M ~ 1
  ),
  Incorr = list(
    y_formula = Y ~ A + X2,
    a_formula = A ~ 1,
    m_formula = M ~ A + A:X3
  )
)

################################################################################
############################ ESTIMATORS #########################################
################################################################################

estimate_psi <- function(df_inference, df_nuisance, bandwidth_silverman, spec) {
  OR_model <- glm(spec$y_formula, family = gaussian, data = df_nuisance)
  propensity_model <- lm(spec$a_formula, data = df_nuisance)
  mediator_model <- glm(spec$m_formula, data = df_nuisance, family = binomial)
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  ps_a <- propensity_score(df_inference, propensity_model, eval_a = a)
  ps_ap <- propensity_score(df_inference, propensity_model, eval_a = a_prime)
  dens_ratio <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  k1 <- dnorm((df_inference$A - a) / bandwidth_silverman) / bandwidth_silverman
  k2 <- dnorm((df_inference$A - a_prime) / bandwidth_silverman) / bandwidth_silverman
  
  term1 <- k1 * ps_a * dens_ratio * (df_inference$Y - OR_prediction)
  term2 <- k2 * ps_ap * (OR_prediction - eta)
  influence_plus_eta <- term1 + term2 + eta
  
  psi_hat <- mean(influence_plus_eta)
  v_hat <- mean((influence_plus_eta - psi_hat)^2)
  
  ww <- k1 * ps_a * dens_ratio
  denom <- sum(ww)
  hest <- if (is.finite(denom) && abs(denom) > .Machine$double.eps) {
    sum(ww * df_inference$Y) / denom
  } else {
    NA_real_
  }
  
  c(psi = psi_hat, var = v_hat, eta = mean(eta), huber = hest)
}

crossfit <- function(df, folds) {
  n_rows <- nrow(df)
  row_partitions <- split(seq_len(n_rows), rep(seq_len(folds), length.out = n_rows))
  spec_names <- names(model_specs)
  
  by_spec <- matrix(NA_real_, nrow = folds, ncol = length(spec_names) * 4L)
  colnames(by_spec) <- as.vector(outer(c("psi", "var", "eta", "huber"), spec_names, paste, sep = "_"))
  
  for (fold_id in seq_len(folds)) {
    partition <- row_partitions[[fold_id]]
    df_inf <- df[partition, , drop = FALSE]
    df_nui <- df[-partition, , drop = FALSE]
    bw <- density(df_inf$A, bw = "nrd0")$bw
    
    fold_out <- lapply(model_specs, function(spec) estimate_psi(df_inf, df_nui, bw, spec))
    for (s in spec_names) {
      by_spec[fold_id, paste0("psi_", s)] <- fold_out[[s]]["psi"]
      by_spec[fold_id, paste0("var_", s)] <- fold_out[[s]]["var"]
      by_spec[fold_id, paste0("eta_", s)] <- fold_out[[s]]["eta"]
      by_spec[fold_id, paste0("huber_", s)] <- fold_out[[s]]["huber"]
    }
  }
  
  psi_estimates <- colMeans(by_spec[, paste0("psi_", spec_names), drop = FALSE], na.rm = TRUE)
  var_estimates <- colMeans(by_spec[, paste0("var_", spec_names), drop = FALSE], na.rm = TRUE)
  eta_estimates <- colMeans(by_spec[, paste0("eta_", spec_names), drop = FALSE], na.rm = TRUE)
  huber_estimates <- colMeans(by_spec[, paste0("huber_", spec_names), drop = FALSE], na.rm = TRUE)
  
  lb <- psi_estimates - qnorm(0.975) * sqrt(var_estimates / n_rows)
  ub <- psi_estimates + qnorm(0.975) * sqrt(var_estimates / n_rows)
  cp <- as.integer(true_val_mc > lb & true_val_mc < ub)
  
  names(psi_estimates) <- paste0("psi_", spec_names)
  names(var_estimates) <- paste0("var_", spec_names)
  names(cp) <- paste0("CP_", spec_names)
  names(eta_estimates) <- paste0("eta_", spec_names)
  names(huber_estimates) <- paste0("huber_", spec_names)
  
  c(psi_estimates, var_estimates, cp, eta_estimates, huber_estimates)
}

run_one_rep <- function(i, sample_size, folds, rep_seeds) {
  set.seed(rep_seeds[i])
  df <- simulate(sample_size)
  crossfit(df, folds)
}

################################################################################
############################ PARALLEL DRIVER ####################################
################################################################################

run_sample_size_parallel <- function(sample_size, n_runs, folds, n_cores, out_dir) {
  message(sprintf("Starting sample_size = %s, n_runs = %s, workers = %s", sample_size, n_runs, n_cores))
  
  set.seed(master_seed + sample_size)
  rep_seeds <- sample.int(.Machine$integer.max, n_runs)
  
  if (n_cores <= 1L) {
    res_list <- lapply(seq_len(n_runs), run_one_rep, sample_size = sample_size, folds = folds, rep_seeds = rep_seeds)
    psi_hat_list <- do.call(rbind, res_list)
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    psi_hat_list <- foreach::foreach(
      i = seq_len(n_runs),
      .combine = rbind,
      .inorder = TRUE,
      .packages = "MASS",
      .export = c(
        "a", "a_prime", "mean_X", "sigma_X", "y1", "y2", "true_val_mc",
        "sigmoid", "simulate", "outcome_regression", "propensity_score",
        "density_ratio", "eta_hat_fn", "model_specs", "estimate_psi",
        "crossfit", "run_one_rep", "rep_seeds"
      )
    ) %dopar% {
      run_one_rep(i, sample_size = sample_size, folds = folds, rep_seeds = rep_seeds)
    }
  }
  
  psi_hat_list <- as.matrix(psi_hat_list)
  storage.mode(psi_hat_list) <- "double"
  
  save_path <- file.path(out_dir, paste0("table1_M1_", sample_size, ".RData"))
  save(psi_hat_list, file = save_path)
  message(sprintf("Saved %s", normalizePath(save_path, mustWork = FALSE)))
  
  psi_hat_list
}

make_summary_table <- function(results_by_n, true_val, rdnum) {
  spec_names <- c("YMA", "YM", "MA", "YA", "Incorr")
  display_spec_names <- c("YMA", "YM", "MA", "YA", "None")
  estimator_blocks <- list(
    Proposal = paste0("psi_", spec_names),
    `CV Eta` = paste0("eta_", spec_names),
    `CV Huber` = paste0("huber_", spec_names)
  )
  
  rows <- list()
  for (sample_size in names(results_by_n)) {
    mat <- results_by_n[[sample_size]]
    for (estimator_name in names(estimator_blocks)) {
      cols <- estimator_blocks[[estimator_name]]
      vals <- vapply(cols, function(col) {
        err <- mat[, col] - true_val
        paste0(
          round(abs(mean(err, na.rm = TRUE)), rdnum),
          " (",
          round(sqrt(mean(err^2, na.rm = TRUE)), rdnum),
          ")"
        )
      }, character(1L))
      rows[[length(rows) + 1L]] <- data.frame(
        n = as.integer(sample_size),
        estimator = estimator_name,
        as.list(setNames(vals, display_spec_names)),
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

################################################################################
############################ RUN ANALYSIS #######################################
################################################################################

results_by_n <- list()
for (sample_size in sample_sizes) {
  results_by_n[[as.character(sample_size)]] <- run_sample_size_parallel(
    sample_size = sample_size,
    n_runs = n_runs,
    folds = folds,
    n_cores = n_cores,
    out_dir = out_dir
  )
}

summary_table <- make_summary_table(results_by_n, true_val_mc, rdnum)
summary_path <- file.path(out_dir, "table1_summary_parallel.csv")
write.csv(summary_table, summary_path, row.names = FALSE)

print(summary_table, row.names = FALSE)
message(sprintf("True value used for bias/RMSE: %.10f", true_val_mc))
message(sprintf("Saved summary table: %s", normalizePath(summary_path, mustWork = FALSE)))
