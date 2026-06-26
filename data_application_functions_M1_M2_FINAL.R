source("fit_cme_conditional_density_cv.R")

################################################################################
##################### Inference Calculators ####################################
################################################################################

epanechnikovKernel <- function(u) {
  ifelse(abs(u) <= 1, (3/4) * (1 - u^2), 0)
}

################################################################################
##################### Model wrappers ###########################################
################################################################################
# mod_option:
#   1 = current parametric models: logistic GLM + Gaussian GLM + beta regression
#   2 = RKHS/Gaussian-kernel working models with CME conditional densities
#
# The middle tree-based branch has been removed. The RKHS/CME branch is now
# relabeled as mod_option = 2.
#
# For mod_option = 2:
#   * kernlab::ksvm is used for the outcome regression and working mean models.
#   * fit_cme_conditional_density() is used for f(d.log | X).
#   * fit_cme_conditional_density() is used for f(logit(m.shifted) | d.log, X).
#   * Mediator density ratios are evaluated on the logit scale. The logit
#     Jacobian cancels in the numerator/denominator ratio.

safe_logit <- function(x, eps = 1e-6) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log(x / (1 - x))
}

safe_expit <- function(z) {
  1 / (1 + exp(-z))
}

################################################################################
##################### Exact RKHS/CME conditional density tools ##################
################################################################################
# The exact RKHS/CME conditional-density engine is sourced from:
#
#   fit_cme_conditional_density_cv.R
#
# That file defines fit_cme_conditional_density(), with internal 3-fold CV for
# lambda by default, and returns predict_density_at() for pointwise conditional
# density evaluation.

fit_nuisance_models <- function(df_nuisance, mod_option = 1,
                                cme_lambda = 1e-3,
                                cme_lambda_grid = c(1e-1, 3e-2, 1e-2,
                                                    3e-3, 1e-3, 3e-4,
                                                    1e-4),
                                cme_cv_lambda = TRUE,
                                cme_cv_folds = 3,
                                cme_cv_seed = 1,
                                cme_pred_chunk_size = 1000) {
  OR_formula <- y.binary ~ . - y - m - d
  propensity_formula <- d.log ~ . - y - m - m.shifted - d - y.binary
  mediator_formula <- m.shifted ~ . - y - y.binary - d - m
  mediator_logit_formula <- m.logit ~ . - y - y.binary - d - m - m.shifted
  
  if (mod_option == 1) {
    OR_model <- glm(OR_formula, data = df_nuisance, family = "binomial")
    
    propensity_model <- glm(propensity_formula, family = gaussian,
                            data = df_nuisance)
    
    mediator_model <- betareg::betareg(mediator_formula, data = df_nuisance)
    mediator_phi <- summary(mediator_model)$coefficients$precision[1]
    
    return(list(mod_option = mod_option,
                OR_model = OR_model,
                propensity_model = propensity_model,
                mediator_model = mediator_model,
                mediator_phi = mediator_phi))
  }
  
  if (mod_option == 2) {
    if (!requireNamespace("kernlab", quietly = TRUE)) {
      stop("Package 'kernlab' is required for mod_option = 2.")
    }
    
    df_or <- df_nuisance
    df_or$y.binary <- factor(df_or$y.binary, levels = c(0, 1))
    OR_model <- kernlab::ksvm(OR_formula, data = df_or,
                              type = "C-svc", kernel = "rbfdot",
                              prob.model = TRUE)
    
    # Keep the SVR working model for compatibility with predict_propensity_mu()
    # and diagnostics. It is not used for f(d.log | X).
    propensity_model <- kernlab::ksvm(propensity_formula, data = df_nuisance,
                                      type = "eps-svr", kernel = "rbfdot")
    prop_mu_train <- as.numeric(predict(propensity_model, newdata = df_nuisance))
    propensity_sigma <- sd(df_nuisance$d.log - prop_mu_train)
    
    df_med <- df_nuisance
    df_med$m.logit <- safe_logit(df_med$m.shifted)
    
    # Keep the SVR working model for compatibility with predict_mediator_mean_or_logit()
    # and diagnostics. The eta Monte Carlo draw and mediator density ratio both use CME.
    mediator_model <- kernlab::ksvm(mediator_logit_formula, data = df_med,
                                    type = "eps-svr", kernel = "rbfdot")
    med_mu_train <- as.numeric(predict(mediator_model, newdata = df_med))
    mediator_sigma <- sd(df_med$m.logit - med_mu_train)
    
    # Strict RKHS/CME estimator for f(d.log | X).
    prop_x_cols <- setdiff(names(df_nuisance),
                           c("d.log", "y", "m", "m.shifted",
                             "d", "y.binary", "m.logit"))
    
    propensity_cme_model <- fit_cme_conditional_density(
      df = df_nuisance,
      y_col = "d.log",
      x_cols = prop_x_cols,
      lambda = cme_lambda,
      lambda_grid = cme_lambda_grid,
      cv_lambda = cme_cv_lambda,
      cv_folds = cme_cv_folds,
      cv_seed = cme_cv_seed,
      verbose = FALSE
    )
    
    # Strict RKHS/CME estimator for f_Z(z | d.log, X), where Z = logit(m.shifted).
    # Fitting on the logit scale avoids boundary leakage from Gaussian density
    # kernels on (0, 1). The density ratio on the original M-scale is identical
    # because the logit Jacobian cancels.
    med_x_cols <- setdiff(names(df_med),
                          c("m.logit", "m.shifted", "y",
                            "y.binary", "d", "m"))
    
    mediator_cme_model <- fit_cme_conditional_density(
      df = df_med,
      y_col = "m.logit",
      x_cols = med_x_cols,
      lambda = cme_lambda,
      lambda_grid = cme_lambda_grid,
      cv_lambda = cme_cv_lambda,
      cv_folds = cme_cv_folds,
      cv_seed = cme_cv_seed,
      verbose = FALSE
    )
    
    return(list(mod_option = mod_option,
                OR_model = OR_model,
                propensity_model = propensity_model,
                propensity_sigma = max(propensity_sigma, 1e-6),
                mediator_model = mediator_model,
                mediator_sigma = max(mediator_sigma, 1e-6),
                mediator_logit = TRUE,
                propensity_cme_model = propensity_cme_model,
                mediator_cme_model = mediator_cme_model,
                propensity_cme_lambda_selected = propensity_cme_model$lambda_selected,
                mediator_cme_lambda_selected = mediator_cme_model$lambda_selected,
                propensity_cme_cv_table = propensity_cme_model$cv_table,
                mediator_cme_cv_table = mediator_cme_model$cv_table,
                cme_pred_chunk_size = cme_pred_chunk_size))
  }
  
  stop("mod_option must be either 1 or 2.")
}

predict_outcome_prob <- function(fit, newdata) {
  if (fit$mod_option == 1) {
    p <- predict(fit$OR_model, newdata = newdata, type = "response")
  } else if (fit$mod_option == 2) {
    pp <- predict(fit$OR_model, newdata = newdata, type = "probabilities")
    p <- pp[, "1"]
  } else {
    stop("Unsupported mod_option.")
  }
  as.numeric(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

predict_propensity_mu <- function(fit, newdata) {
  if (fit$mod_option == 1) {
    as.numeric(predict(fit$propensity_model, newdata = newdata, type = "response"))
  } else if (fit$mod_option == 2) {
    as.numeric(predict(fit$propensity_model, newdata = newdata))
  } else {
    stop("Unsupported mod_option.")
  }
}

predict_mediator_mean_or_logit <- function(fit, newdata) {
  if (fit$mod_option == 1) {
    as.numeric(predict(fit$mediator_model, newdata = newdata, type = "response"))
  } else if (fit$mod_option == 2) {
    as.numeric(predict(fit$mediator_model, newdata = newdata))
  } else {
    stop("Unsupported mod_option.")
  }
}

outcome_regression <- function(df, fit, a) {
  df_a <- df
  df_a$d.log <- a
  predict_outcome_prob(fit, df_a)
}

hajek_propensity_score_corrected <- function(df_inference, df_nuisance, fit, a, cc = 0.01) {
  if (fit$mod_option == 2) {
    # RKHS/CME estimate of f(d.log = a | X).
    ps <- fit$propensity_cme_model$predict_density_at(
      new_df = df_inference,
      y_eval = rep(a, nrow(df_inference)),
      pred_chunk_size = fit$cme_pred_chunk_size
    )
    
    ps_nc <- fit$propensity_cme_model$predict_density_at(
      new_df = df_nuisance,
      y_eval = rep(a, nrow(df_nuisance)),
      pred_chunk_size = fit$cme_pred_chunk_size
    )
  } else if (fit$mod_option == 1) {
    mu <- predict_propensity_mu(fit, df_inference)
    sigma <- sqrt(var(residuals(fit$propensity_model)))
    ps <- dnorm(a, mu, sigma)
    
    mu_nc <- predict_propensity_mu(fit, df_nuisance)
    ps_nc <- dnorm(a, mu_nc, sigma)
  } else {
    stop("Unsupported mod_option.")
  }
  
  bw <- density(df_nuisance$d.log, bw = "nrd0")$bw
  hajek_kernel <- (1 / bw) * dnorm((df_nuisance$d.log - a) / bw, mean = 0, sd = 1)
  normalizing_constant <- mean(hajek_kernel / pmax(ps_nc, 1e-12))
  
  normalized_ps <- ps * normalizing_constant
  print(mean(normalized_ps < cc))
  normalized_ps[normalized_ps < cc] <- cc
  return(1 / normalized_ps)
}

density_ratio <- function(df, fit, a, a_prime) {
  df_a_prime <- df
  df_a_prime$d.log <- a_prime
  df_a <- df
  df_a$d.log <- a
  
  if (fit$mod_option == 1) {
    mu_a <- predict_mediator_mean_or_logit(fit, df_a)
    mu_a_prime <- predict_mediator_mean_or_logit(fit, df_a_prime)
    phi <- fit$mediator_phi
    
    # Extra safety for numerical stability in dbeta.
    mu_a <- pmin(pmax(mu_a, 1e-6), 1 - 1e-6)
    mu_a_prime <- pmin(pmax(mu_a_prime, 1e-6), 1 - 1e-6)
    
    return(dbeta(df$m.shifted, mu_a_prime * phi, (1 - mu_a_prime) * phi) /
             pmax(dbeta(df$m.shifted, mu_a * phi, (1 - mu_a) * phi), 1e-12))
  } else if (fit$mod_option == 2) {
    # RKHS/CME estimate of
    #   f_Z(logit(m.shifted) | d.log = a_prime, X) /
    #   f_Z(logit(m.shifted) | d.log = a,       X),
    # where Z = logit(M). The logit Jacobian cancels in the M-density ratio.
    z <- safe_logit(df$m.shifted)
    
    f_z_aprime <- fit$mediator_cme_model$predict_density_at(
      new_df = df_a_prime,
      y_eval = z,
      pred_chunk_size = fit$cme_pred_chunk_size
    )
    
    f_z_a <- fit$mediator_cme_model$predict_density_at(
      new_df = df_a,
      y_eval = z,
      pred_chunk_size = fit$cme_pred_chunk_size
    )
    
    return(f_z_aprime / pmax(f_z_a, 1e-12))
  } else {
    stop("Unsupported mod_option.")
  }
}

.sample_from_density_grid <- function(grid, density, nsim, eps = 1e-12) {
  grid <- as.numeric(grid)
  density <- pmax(as.numeric(density), eps)
  
  if (length(grid) != length(density)) {
    stop("grid and density must have the same length.")
  }
  if (length(grid) < 2) {
    stop("grid must have length at least 2.")
  }
  
  # Approximate probability mass on the grid. The common grid-spacing factor
  # cancels after normalization for an evenly spaced grid.
  prob <- density / sum(density)
  if (!all(is.finite(prob)) || sum(prob) <= 0) {
    prob <- rep(1 / length(grid), length(grid))
  }
  
  sample(grid, size = nsim, replace = TRUE, prob = prob)
}

monte_carlo_integration <- function(row, fit, a, nsim_m = 1000,
                                    eta_density_n_grid = 300,
                                    eta_density_grid_probs = c(0.001, 0.999)) {
  # row should already have d.log set to the mediator-distribution treatment
  # level a_prime. This function then simulates M from that conditional mediator
  # distribution and evaluates the outcome model at d.log = a.
  row_df <- as.data.frame(row, stringsAsFactors = FALSE)
  row_df[] <- lapply(row_df, function(x) type.convert(as.character(x), as.is = TRUE))
  
  if (fit$mod_option == 1) {
    m_samples <- rbeta(
      nsim_m,
      row_df$mu * fit$mediator_phi,
      (1 - row_df$mu) * fit$mediator_phi
    )
  } else if (fit$mod_option == 2) {
    # CME draw from f_Z(z | d.log = a_prime, X), where Z = logit(M),
    # then transform back to the mediator scale. This avoids boundary leakage
    # from fitting Gaussian density kernels directly on m.shifted in (0, 1).
    dens_fit <- fit$mediator_cme_model$predict_density_grid(
      new_df = row_df,
      n_grid = eta_density_n_grid,
      grid_probs = eta_density_grid_probs,
      pred_chunk_size = fit$cme_pred_chunk_size,
      normalize = TRUE,
      nonnegative = TRUE
    )
    
    z_grid <- if (!is.null(dens_fit$y_grid)) dens_fit$y_grid else dens_fit$A_grid
    
    z_samples <- .sample_from_density_grid(
      grid = z_grid,
      density = dens_fit$density[1, ],
      nsim = nsim_m
    )
    
    m_samples <- safe_expit(z_samples)
  } else {
    stop("Unsupported mod_option.")
  }
  
  rep.row_df <- row_df[rep(1, each = nsim_m), ]
  rep.row_df$m.shifted <- m_samples
  rep.row_df$d.log <- a
  
  OR_predictions <- predict_outcome_prob(fit, rep.row_df)
  mean(OR_predictions)
}

eta_hat_fn <- function(df, fit, a, a_prime, nsim_m = 1000,
                       eta_density_n_grid = 300,
                       eta_density_grid_probs = c(0.001, 0.999)) {
  # Estimate eta(a, a', X) = E[ h(a, M, X) | A = a', X ].
  # For mod_option = 2, mediator draws come from the fitted CME conditional
  # density estimator, not from a Gaussian residual approximation.
  df_a_prime <- df
  df_a_prime$d.log <- a_prime
  
  if (fit$mod_option == 1) {
    mu <- predict_mediator_mean_or_logit(fit, df_a_prime)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    df_a_prime$mu <- mu
  }
  
  eta_hat <- vapply(
    seq_len(nrow(df_a_prime)),
    function(i) {
      monte_carlo_integration(
        row = df_a_prime[i, , drop = FALSE],
        fit = fit,
        a = a,
        nsim_m = nsim_m,
        eta_density_n_grid = eta_density_n_grid,
        eta_density_grid_probs = eta_density_grid_probs
      )
    },
    numeric(1)
  )
  
  eta_hat
}

################################################################################
##################### Estimators ###############################################
################################################################################

estimate_nde_nie <- function(df_inference, df_nuisance, a, a_prime,
                             bdw, cc,
                             mod_option = 1,
                             nsim_m = 1000,
                             eta_density_n_grid = 300,
                             eta_density_grid_probs = c(0.001, 0.999),
                             cme_lambda = 1e-3,
                             cme_lambda_grid = c(1e-1, 3e-2, 1e-2,
                                                 3e-3, 1e-3, 3e-4,
                                                 1e-4),
                             cme_cv_lambda = TRUE,
                             cme_cv_folds = 3,
                             cme_cv_seed = 1,
                             cme_pred_chunk_size = 1000) {
  
  bw_silverman <- bdw
  
  fit <- fit_nuisance_models(df_nuisance,
                             mod_option = mod_option,
                             cme_lambda = cme_lambda,
                             cme_lambda_grid = cme_lambda_grid,
                             cme_cv_lambda = cme_cv_lambda,
                             cme_cv_folds = cme_cv_folds,
                             cme_cv_seed = cme_cv_seed,
                             cme_pred_chunk_size = cme_pred_chunk_size)
  
  OR_prediction <- outcome_regression(df = df_inference, fit = fit, a = a)
  
  ps_a <- hajek_propensity_score_corrected(df_inference, df_nuisance,
                                           fit = fit, a = a, cc = cc)
  ps_a_prime <- hajek_propensity_score_corrected(df_inference, df_nuisance,
                                                 fit = fit, a = a_prime, cc = cc)
  
  densratio_prediction <- density_ratio(df = df_inference, fit = fit,
                                        a = a, a_prime = a_prime)
  
  eta <- eta_hat_fn(df = df_inference, fit = fit,
                    a = a, a_prime = a_prime, nsim_m = nsim_m,
                    eta_density_n_grid = eta_density_n_grid,
                    eta_density_grid_probs = eta_density_grid_probs)
  
  if_kernel_1 <- (1 / bw_silverman) *
    dnorm((df_inference$d.log - a) / bw_silverman, mean = 0, sd = 1)
  
  if_term_1 <- if_kernel_1 * ps_a *
    densratio_prediction * (df_inference$y.binary - OR_prediction)
  
  if_kernel_2 <- (1 / bw_silverman) *
    dnorm((df_inference$d.log - a_prime) / bw_silverman, mean = 0, sd = 1)
  
  if_term_2 <- if_kernel_2 * ps_a_prime * (OR_prediction - eta)
  
  psi_a_aprime <- mean(if_term_1 + if_term_2 + eta)
  
  eta_a_prime <- eta_hat_fn(df = df_inference, fit = fit,
                            a = a_prime, a_prime = a_prime, nsim_m = nsim_m,
                            eta_density_n_grid = eta_density_n_grid,
                            eta_density_grid_probs = eta_density_grid_probs)
  
  if_kernel_3 <- if_kernel_2
  if_term_3 <- if_kernel_3 * ps_a_prime *
    (df_inference$y.binary - eta_a_prime)
  
  psi_aprime <- mean(if_term_3 + eta_a_prime)
  
  eta_a <- eta_hat_fn(df = df_inference, fit = fit,
                      a = a, a_prime = a, nsim_m = nsim_m,
                      eta_density_n_grid = eta_density_n_grid,
                      eta_density_grid_probs = eta_density_grid_probs)
  
  if_kernel_4 <- if_kernel_1
  if_term_4 <- if_kernel_4 * ps_a *
    (df_inference$y.binary - eta_a)
  
  psi_a <- mean(if_term_4 + eta_a)
  
  nde <- psi_a_aprime - psi_aprime
  vhat_nde <- mean(((if_term_1 + if_term_2 + eta) -
                      (if_term_3 + eta_a_prime) - nde)^2)
  
  nie <- psi_a - psi_a_aprime
  vhat_nie <- mean(((if_term_4 + eta_a) -
                      (if_term_1 + if_term_2 + eta) - nie)^2)
  
  return(c(nde, vhat_nde, nie, vhat_nie))
}

crossfit <- function(df, folds, a, a_prime,
                     mod_option = 1,
                     nsim_m = 1000,
                     eta_density_n_grid = 300,
                     eta_density_grid_probs = c(0.001, 0.999),
                     cc = 0.01,
                     cme_lambda = 1e-3,
                     cme_lambda_grid = c(1e-1, 3e-2, 1e-2,
                                         3e-3, 1e-3, 3e-4,
                                         1e-4),
                     cme_cv_lambda = TRUE,
                     cme_cv_folds = 3,
                     cme_cv_seed = 1,
                     cme_pred_chunk_size = 1000) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds,
                                        labels = FALSE))
  ndenie_hat_folds <- matrix(NA, nrow = folds, ncol = 4)
  
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    dinf <- df[partition, ]
    bw_silverman <- density(dinf$d.log, bw = "nrd0")$bw
    
    ndenie_hat_folds[i, ] <- estimate_nde_nie(
      df_inference = dinf,
      df_nuisance = df[-partition, ],
      a = a,
      a_prime = a_prime,
      bdw = bw_silverman,
      cc = cc,
      mod_option = mod_option,
      nsim_m = nsim_m,
      eta_density_n_grid = eta_density_n_grid,
      eta_density_grid_probs = eta_density_grid_probs,
      cme_lambda = cme_lambda,
      cme_lambda_grid = cme_lambda_grid,
      cme_cv_lambda = cme_cv_lambda,
      cme_cv_folds = cme_cv_folds,
      cme_cv_seed = cme_cv_seed,
      cme_pred_chunk_size = cme_pred_chunk_size
    )
  }
  
  NDE_psi_hat <- mean(ndenie_hat_folds[, 1])
  NDE_v_hat <- mean(ndenie_hat_folds[, 2])
  
  NDE_confidence_interval_lb <- NDE_psi_hat - qnorm(0.975) * sqrt(NDE_v_hat / n_rows)
  NDE_confidence_interval_ub <- NDE_psi_hat + qnorm(0.975) * sqrt(NDE_v_hat / n_rows)
  
  NIE_psi_hat <- mean(ndenie_hat_folds[, 3])
  NIE_v_hat <- mean(ndenie_hat_folds[, 4])
  
  NIE_confidence_interval_lb <- NIE_psi_hat - qnorm(0.975) * sqrt(NIE_v_hat / n_rows)
  NIE_confidence_interval_ub <- NIE_psi_hat + qnorm(0.975) * sqrt(NIE_v_hat / n_rows)
  
  return(list(NDE = c(NDE_confidence_interval_lb, NDE_confidence_interval_ub,
                      NDE_v_hat, NDE_psi_hat),
              NIE = c(NIE_confidence_interval_lb, NIE_confidence_interval_ub,
                      NIE_v_hat, NIE_psi_hat)))
}

crossfit_sensitivity <- function(df, folds, a, a_prime, bdw, cc,
                                 mod_option = 1,
                                 nsim_m = 1000,
                                 eta_density_n_grid = 300,
                                 eta_density_grid_probs = c(0.001, 0.999),
                                 cme_lambda = 1e-3,
                                 cme_lambda_grid = c(1e-1, 3e-2, 1e-2,
                                                     3e-3, 1e-3, 3e-4,
                                                     1e-4),
                                 cme_cv_lambda = TRUE,
                                 cme_cv_folds = 3,
                                 cme_cv_seed = 1,
                                 cme_pred_chunk_size = 1000) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds,
                                        labels = FALSE))
  ndenie_hat_folds <- matrix(NA, nrow = folds, ncol = 4)
  
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    ndenie_hat_folds[i, ] <- estimate_nde_nie(
      df_inference = df[partition, ],
      df_nuisance = df[-partition, ],
      a = a,
      a_prime = a_prime,
      bdw = bdw,
      cc = cc,
      mod_option = mod_option,
      nsim_m = nsim_m,
      eta_density_n_grid = eta_density_n_grid,
      eta_density_grid_probs = eta_density_grid_probs,
      cme_lambda = cme_lambda,
      cme_lambda_grid = cme_lambda_grid,
      cme_cv_lambda = cme_cv_lambda,
      cme_cv_folds = cme_cv_folds,
      cme_cv_seed = cme_cv_seed,
      cme_pred_chunk_size = cme_pred_chunk_size
    )
  }
  
  NDE_psi_hat <- mean(ndenie_hat_folds[, 1])
  NDE_v_hat <- mean(ndenie_hat_folds[, 2])
  
  NDE_confidence_interval_lb <- NDE_psi_hat - qnorm(0.975) * sqrt(NDE_v_hat / n_rows)
  NDE_confidence_interval_ub <- NDE_psi_hat + qnorm(0.975) * sqrt(NDE_v_hat / n_rows)
  
  NIE_psi_hat <- mean(ndenie_hat_folds[, 3])
  NIE_v_hat <- mean(ndenie_hat_folds[, 4])
  
  NIE_confidence_interval_lb <- NIE_psi_hat - qnorm(0.975) * sqrt(NIE_v_hat / n_rows)
  NIE_confidence_interval_ub <- NIE_psi_hat + qnorm(0.975) * sqrt(NIE_v_hat / n_rows)
  
  return(list(NDE = c(NDE_confidence_interval_lb, NDE_confidence_interval_ub,
                      NDE_v_hat, NDE_psi_hat),
              NIE = c(NIE_confidence_interval_lb, NIE_confidence_interval_ub,
                      NIE_v_hat, NIE_psi_hat)))
}
