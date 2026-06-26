################################################################################
##################### Data Application Figure 2: 2 Model Options ################
##################### Parallelized over d_list within each mod_option ############
################################################################################
# Parallelization follows the same RStudio/standard-R strategy used in
# Simulation_final_Table1_parallel.R: foreach + doParallel with n_cores = 19.
# This script assumes that JCdata has already been loaded and cleaned as in
# data application main.R, including the variables m.shifted, d.log, and y.binary.
# It also assumes that path is defined. If path is not defined, the working
# directory is used for output.
#
# Current model options:
#   mod_option = 1: parametric nuisance models
#   mod_option = 2: RKHS/CME nuisance-density models
#
# The middle tree-based branch has been removed. The RKHS/CME branch is now
# relabeled as mod_option = 2.

set.seed(2025)

if (!exists("JCdata")) {
  stop("JCdata is not available. Run the data-cleaning part of data application main.R first.")
}
if (!exists("crossfit") || !is.function(crossfit)) {
  stop("crossfit() is not available. Source the revised functions script before running this script.")
}
if (!all(c("m.shifted", "d.log", "y.binary") %in% colnames(JCdata))) {
  stop("JCdata must contain m.shifted, d.log, and y.binary before running this script.")
}
if (!exists("path")) {
  path <- getwd()
}

# Required packages by option.
if (!requireNamespace("betareg", quietly = TRUE)) {
  stop("Package 'betareg' is required for mod_option = 1.")
}
if (!requireNamespace("kernlab", quietly = TRUE)) {
  stop("Package 'kernlab' is required for mod_option = 2.")
}

# Required for the parallel driver below, following Simulation_final_Table1_parallel.R.
required_parallel_packages <- c("foreach", "doParallel")
missing_parallel_packages <- required_parallel_packages[
  !vapply(required_parallel_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_parallel_packages) > 0L) {
  stop(
    "Please install missing parallel packages first: install.packages(c(",
    paste(sprintf('"%s"', missing_parallel_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

################################################################################
##################### Analysis parameters ######################################
################################################################################

folds <- 3
# d_list <- c(100, 1000, 2000)
d_list <- seq(100, 2000, 100)
d_prime <- 60

# Increase nsim_m for final runs if Monte Carlo noise in eta_hat_fn() is visible.
nsim_m <- 1000
cc <- 0.01

# Use the same core setting strategy as Simulation_final_Table1_parallel.R.
# On shared servers, manually set this lower, e.g. n_cores <- 8L.
n_cores <- 19L # max(1L, parallel::detectCores(logical = FALSE) - 1L)

# Reproducibility seed for parallel treatment-level jobs.
master_seed <- 2025L

mod_options <- 1:2
n_mod <- length(mod_options)
n_d <- length(d_list)

# Arrays: rows = treatment levels, columns = CI lower, CI upper, variance, estimate;
# third dimension = model option.
resNDE <- array(NA_real_, dim = c(n_d, 4, n_mod),
                dimnames = list(paste0("d_", d_list),
                                c("lb", "ub", "vhat", "psi_hat"),
                                paste0("M", mod_options)))
resNIE <- array(NA_real_, dim = c(n_d, 4, n_mod),
                dimnames = list(paste0("d_", d_list),
                                c("lb", "ub", "vhat", "psi_hat"),
                                paste0("M", mod_options)))

################################################################################
##################### Run cross-fitting in parallel over d_list #################
################################################################################

run_one_d_application <- function(i, d_list, d_prime, mod_option_here,
                                  JCdata, folds, nsim_m, cc, d_seeds) {
  set.seed(d_seeds[i])
  
  d_here <- d_list[i]
  a <- log(d_here)
  a_prime <- log(d_prime)
  
  tmp <- crossfit(df = JCdata,
                  folds = folds,
                  a = a,
                  a_prime = a_prime,
                  mod_option = mod_option_here,
                  nsim_m = nsim_m,
                  cc = cc)
  
  list(i = i, d = d_here, NDE = tmp$NDE, NIE = tmp$NIE)
}

parallel_packages <- unique(c("betareg", "kernlab"))

run_d_list_for_one_model <- function(m_ind, mod_option_here) {
  message("Running mod_option = ", mod_option_here)
  for (i in seq_along(d_list)) {
    message("  Queued treatment a = ", d_list[i], ", reference a' = ", d_prime)
  }
  
  set.seed(master_seed + as.integer(mod_option_here))
  d_seeds <- sample.int(.Machine$integer.max, n_d)
  
  if (n_cores <= 1L) {
    out_list <- lapply(
      seq_along(d_list),
      run_one_d_application,
      d_list = d_list,
      d_prime = d_prime,
      mod_option_here = mod_option_here,
      JCdata = JCdata,
      folds = folds,
      nsim_m = nsim_m,
      cc = cc,
      d_seeds = d_seeds
    )
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    # Export all currently defined global objects, including hidden helper
    # functions whose names start with a dot. Plain ls() omits these names,
    # so all.names = TRUE is essential for PSOCK workers.
    parallel_export_names <- ls(envir = .GlobalEnv, all.names = TRUE)
    parallel_export_names <- setdiff(parallel_export_names, c(".Random.seed", ".Last.value"))
    parallel_export_names <- parallel_export_names[
      vapply(parallel_export_names, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
    ]
    parallel::clusterExport(cl, varlist = parallel_export_names, envir = .GlobalEnv)
    
    out_list <- foreach::foreach(
      i = seq_along(d_list),
      .inorder = TRUE,
      .packages = parallel_packages,
      .export = parallel_export_names,
      .errorhandling = "stop"
    ) %dopar% {
      run_one_d_application(i = i,
                            d_list = d_list,
                            d_prime = d_prime,
                            mod_option_here = mod_option_here,
                            JCdata = JCdata,
                            folds = folds,
                            nsim_m = nsim_m,
                            cc = cc,
                            d_seeds = d_seeds)
    }
  }
  
  out_list
}

for (m_ind in seq_along(mod_options)) {
  mod_option_here <- mod_options[m_ind]
  out_list <- run_d_list_for_one_model(m_ind = m_ind, mod_option_here = mod_option_here)
  
  for (out in out_list) {
    resNDE[out$i, , m_ind] <- out$NDE
    resNIE[out$i, , m_ind] <- out$NIE
  }
}

save(resNDE, resNIE, d_list, d_prime, folds, mod_options, n_cores,
     file = file.path(path, "CM_application_Figure2_mod_options_M1_M2_results.RData"))

################################################################################
##################### Plot helper ##############################################
################################################################################

plot_effect <- function(d_list, res_mat, effect_name, model_index,
                        ylim = NULL, xlab = "Treatment a") {
  est <- res_mat[, "psi_hat"]
  lb <- res_mat[, "lb"]
  ub <- res_mat[, "ub"]
  
  if (is.null(ylim)) {
    ylim <- range(c(lb, ub, 0), na.rm = TRUE)
  }
  
  plot(d_list, est, type = "n",
       ylab = effect_name,
       xlab = xlab,
       ylim = ylim,
       main = bquote(M[.(model_index)] ~ .(effect_name)))
  polygon(c(d_list, rev(d_list)),
          c(lb, rev(ub)),
          col = "grey75", border = FALSE)
  lines(d_list, est, lwd = 2)
  points(d_list, est, pch = 16)
  abline(h = 0, lty = 2)
}

# Use shared y-limits across both model options within each estimand so the panels
# are directly comparable.
ylim_nde <- range(c(resNDE[, "lb", ], resNDE[, "ub", ], 0), na.rm = TRUE)
ylim_nie <- range(c(resNIE[, "lb", ], resNIE[, "ub", ], 0), na.rm = TRUE)

################################################################################
##################### 2 x 2 Figure #############################################
################################################################################

pdf(file.path(path, "CM_application_Hajek_NDE_NIE_mod_options_M1_M2_2x2.pdf"),
    width = 10, height = 8)

op <- par(
  mfrow = c(2, 2),
  mar = c(4.2, 4.2, 3.2, 1.2),
  oma = c(0, 0, 1, 0)
)

# Row 1: M1 parametric specification.
plot_effect(d_list, resNDE[, , 1], "Natural Direct Effect", 1, ylim = ylim_nde)
plot_effect(d_list, resNIE[, , 1], "Natural Indirect Effect", 1, ylim = ylim_nie)

# Row 2: M2 RKHS/CME specification.
plot_effect(d_list, resNDE[, , 2], "Natural Direct Effect", 2, ylim = ylim_nde)
plot_effect(d_list, resNIE[, , 2], "Natural Indirect Effect", 2, ylim = ylim_nie)

mtext(
  "Natural direct and indirect effect estimates under different nuisance-model specifications",
  outer = TRUE,
  cex = 1.1,
  font = 2
)

par(op)
dev.off()

################################################################################
##################### Optional PNG output ######################################
################################################################################

# png(file.path(path, "CM_application_Hajek_NDE_NIE_mod_options_M1_M2_2x2.png"),
#     width = 2000, height = 1600, res = 200)
# 
# op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.2), oma = c(0, 0, 1, 0))
# 
# plot_effect(d_list, resNDE[, , 1], "Natural Direct Effect", 1, ylim = ylim_nde)
# plot_effect(d_list, resNIE[, , 1], "Natural Indirect Effect", 1, ylim = ylim_nie)
# plot_effect(d_list, resNDE[, , 2], "Natural Direct Effect", 2, ylim = ylim_nde)
# plot_effect(d_list, resNIE[, , 2], "Natural Indirect Effect", 2, ylim = ylim_nie)
# 
# mtext("Natural direct and indirect effect estimates under different nuisance-model specifications",
#       outer = TRUE, cex = 1.1, font = 2)
# 
# par(op)
# dev.off()
