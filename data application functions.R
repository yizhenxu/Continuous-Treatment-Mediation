
################################################################################
##################### Inference Calculators ####################################
################################################################################
epanechnikovKernel <- function(u) {
  # Check if |u| <= 1
  ifelse(abs(u) <= 1, (3/4) * (1 - u^2), 0)
}

outcome_regression <- function(df, model, a){
  df_a <- df
  df_a$d.log <- a
  y.pred <- predict(model, newdata = df_a, type = "response")
  return(y.pred)
}


# Propensity score with parametric density estimation
hajek_propensity_score <- function(df, model, a){
  mu <- predict(model, newdata = df, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  bw <- density(df$d.log, bw = "nrd0")$bw
  hajek_kernel <- (1/bw) * dnorm((df$d.log - a)/bw, mean = 0, sd = 1)
  normalizing_constant <- mean(hajek_kernel/ps)
  normalized_ps <- ps * normalizing_constant
  
  if(1){ # clipping
    cc = 0.01
    print(mean(normalized_ps < cc))
    
    normalized_ps[normalized_ps < cc] <- cc
    return(1/normalized_ps)
  }
  
  if(0){ # shifting
    print(summary(normalized_ps))
    normalized_ps = normalized_ps + 0.02
    return(1/normalized_ps)
  }
  
}

hajek_propensity_score_corrected <- function(df_inference, df_nuisance, model, a, cc = 0.01){
  # propensity based on df_inference I(l)
  mu <- predict(model, newdata = df_inference, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  
  # calculate normalizing constant based on df_nuisance I(-l)
  mu_nc <- predict(model, newdata = df_nuisance, type = "response")
  ps_nc <- dnorm(a, mu_nc, sigma)
  bw <- density(df_nuisance$d.log, bw = "nrd0")$bw
  hajek_kernel <- (1/bw) * dnorm((df_nuisance$d.log - a)/bw, mean = 0, sd = 1)
  normalizing_constant <- mean(hajek_kernel/ps_nc)
  
  normalized_ps <- ps * normalizing_constant
  
  #cc = 0.01 # clipping threshold
  print(mean(normalized_ps < cc))
  
  normalized_ps[normalized_ps < cc] <- cc
  return(1/normalized_ps)
  
}

# Propensity score with parametric density estimation
propensity_score <- function(df, model, a, cc = 0.01){
  mu <- predict(model, newdata = df, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  
  ##cc = 0.01
  print(mean(ps < cc))
  
  ps[ps < cc] <- cc # clipping
  
  # ps <- ps + 2/sqrt(nrow(df)) # shifting
  # ps[ps < 0.05] <- 0.05 # clipping
  return(1/ps)
}

# Density ratio
density_ratio <- function(df, model, phi, a, a_prime){
  df_a_prime <- df
  df_a_prime$d.log <- a_prime
  df_a <- df
  df_a$d.log <- a
  
  mu_a <- predict(model, newdata = df_a, type = "response")
  mu_a_prime <- predict(model, newdata = df_a_prime, type = "response")
  
  return(dbeta(df$m.shifted, mu_a_prime*phi, (1 - mu_a_prime)*phi)/
           dbeta(df$m.shifted, mu_a*phi, (1 - mu_a)*phi))
}

monte_carlo_integration <- function(row, OR_model, mediator_model,  phi, a){
  # Convert current row into df
  row_df <- as.data.frame(t(row), stringsAsFactors = FALSE)
  # Get predicted mu_i from beta regression
  #mu_i <- predict(mediator_model, new_data = row_df, type = "response")
  # Sample mediator values from beta distribution for corresponding row
  m_samples <- rbeta(1000, row_df$mu*phi, (1 - row_df$mu) * phi)
  
  #rep.row_df <- do.call(rbind, replicate(1000, row_df, simplify = FALSE))
  # Update replicated data frame
  rep.row_df <- row_df[rep(1, each = 1000), ]
  rep.row_df$m.shifted <- m_samples
  rep.row_df$d.log <- a
  
  OR_predictions <- predict(OR_model, newdata = rep.row_df, type = "response")
  # Take average of OR to obtain eta estimate
  mean(OR_predictions)
}

eta_hat_fn <- function(df, mediator_model, mediator_phi, OR_model, a, a_prime){
  df_a_prime <- df
  df_a_prime$d.log <- a_prime
  mu <- predict(mediator_model, newdata = df_a_prime, type = "response")
  df_a_prime$mu <- mu
  eta_hat <- apply(df_a_prime, 1, monte_carlo_integration, 
                   mediator_model = mediator_model, 
                   OR_model = OR_model,
                   a = a, 
                   phi = mediator_phi)
  eta_hat
}

################################################################################
##################### Estimators ###############################################
################################################################################


estimate_nde_nie <- function(df_inference, df_nuisance, a, a_prime, bdw, cc) {

  bw_silverman = bdw
  #bw_silverman <- density(df_inference$d.log, bw = "nrd0")$bw * nrow(df_inference)^(-1 * delta/5)
  
  OR_model <-  glm(y.binary ~ . - y - m - d, data = df_nuisance, 
                   family = "binomial")
  
  propensity_model <- glm(d.log ~ . - y - m - m.shifted - d - y.binary,
                          family = gaussian, data = df_nuisance)
  
  #propensity_model <- glm(d.log ~ .^2, family = gaussian, data = df_nuisance[,c(4:68, 70)])
  
  mediator_model <- betareg(m.shifted ~ .- y - y.binary - d - m,
                            data = df_nuisance)
  mediator_phi <- summary(mediator_model)$coefficients$precision[1]
  
  # Mediation predictions
  OR_prediction <- outcome_regression(df = df_inference,
                                      model = OR_model, a = a)
  
  # Hajek
  ps_a <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a, cc)
  ps_a_prime <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a_prime, cc)
  # Non-Hajek
  #ps_a <- propensity_score(df = df_inference, model = propensity_model, a = a, cc)
  #ps_a_prime <- propensity_score(df = df_inference, model = propensity_model, a = a_prime, cc)
  
  densratio_prediction <- density_ratio(df = df_inference,
                                        model = mediator_model,
                                        phi = mediator_phi,
                                        a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df = df_inference, mediator_model = mediator_model,
                    mediator_phi = mediator_phi, OR_model = OR_model,
                    a = a, a_prime = a_prime)
  
  # Mediator counterfactual psi(a,a')
  if_kernel_1 = (1/bw_silverman) *
    dnorm((df_inference$d.log - a)/bw_silverman, mean = 0, sd = 1)
  
  if_term_1 = if_kernel_1 * ps_a *
    densratio_prediction * (df_inference$y.binary - OR_prediction)
  
  if_kernel_2 = (1/bw_silverman) *
    dnorm((df_inference$d.log - a_prime)/bw_silverman, mean = 0, sd = 1)
  
  if_term_2 = if_kernel_2 * ps_a_prime * (OR_prediction - eta)
  
  psi_a_aprime = mean(if_term_1 + if_term_2 + eta)
  
  # ATE counterfactual psi(a', a')
  eta_a_prime <- eta_hat_fn(df = df_inference, mediator_model = mediator_model,
                            mediator_phi = mediator_phi, OR_model = OR_model,
                            a = a_prime, a_prime = a_prime)
  if_kernel_3 = if_kernel_2
  #(1/bw_silverman) * 
  #dnorm((df_inference$d.log - a_prime)/bw_silverman, mean = 0, sd = 1)
  
  if_term_3 = if_kernel_3 * ps_a_prime *
    (df_inference$y.binary - eta_a_prime)
  
  psi_aprime = mean(if_term_3 + eta_a_prime)
  
  # ATE counterfactual psi(a, a)
  eta_a <- eta_hat_fn(df = df_inference, mediator_model = mediator_model,
                      mediator_phi = mediator_phi, OR_model = OR_model,
                      a = a, a_prime = a)
  
  if_kernel_4 = if_kernel_1
  
  if_term_4 = if_kernel_4 * ps_a *
    (df_inference$y.binary - eta_a)
  
  psi_a = mean(if_term_4 + eta_a)
  
  # NDE, NIE
  nde = psi_a_aprime - psi_aprime
  vhat_nde =  mean(( (if_term_1 + if_term_2 + eta) - (if_term_3 + eta_a_prime) - nde)^2) 
  
  nie = psi_a - psi_a_aprime
  vhat_nie =  mean(( (if_term_4 + eta_a) - (if_term_1 + if_term_2 + eta)  - nie)^2) 
  #cancel *bw_silverman with denominator later
  
  return(c(nde, vhat_nde, nie, vhat_nie))
}

# Cross fitting scheme -- optimal bandwidth (Silverman's) and 0.01 clipping
crossfit <- function(df, folds, a, a_prime) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds, 
                                        labels = FALSE))
  ndenie_hat_folds <- matrix(NA, nrow = folds, ncol = 4)
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    
    # Calculate bandwidth using Silverman's rule of thumb 
    dinf = df[partition, ]
    bw_silverman <- density(dinf$d.log, bw = "nrd0")$bw
    
    ndenie_hat_folds[i,] <- estimate_nde_nie(df_inference = dinf, 
                                             df_nuisance = df[-partition, ],
                                             a = a, 
                                             a_prime = a_prime, 
                                             bdw = bw_silverman, 
                                             cc = 0.01)
  }
  
  NDE_psi_hat <- mean(ndenie_hat_folds[,1])
  NDE_v_hat <- mean(ndenie_hat_folds[,2])
  
  NDE_confidence_interval_lb <- NDE_psi_hat - qnorm(0.975) * sqrt(NDE_v_hat/n_rows)
  NDE_confidence_interval_ub <- NDE_psi_hat + qnorm(0.975) * sqrt(NDE_v_hat/n_rows)
  
  NIE_psi_hat <- mean(ndenie_hat_folds[,3])
  NIE_v_hat <- mean(ndenie_hat_folds[,4])
  
  NIE_confidence_interval_lb <- NIE_psi_hat - qnorm(0.975) * sqrt(NIE_v_hat/n_rows)
  NIE_confidence_interval_ub <- NIE_psi_hat + qnorm(0.975) * sqrt(NIE_v_hat/n_rows)
  
  return(list(NDE = c(NDE_confidence_interval_lb, NDE_confidence_interval_ub, NDE_v_hat, NDE_psi_hat),
              NIE = c(NIE_confidence_interval_lb, NIE_confidence_interval_ub, NIE_v_hat, NIE_psi_hat)) )
}

# Cross fitting scheme -- bandwidth bdw and clipping cc
crossfit_sensitivity <- function(df, folds, a, a_prime, bdw, cc) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds, 
                                        labels = FALSE))
  ndenie_hat_folds <- matrix(NA, nrow = folds, ncol = 4)
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    ndenie_hat_folds[i,] <- estimate_nde_nie(df_inference = df[partition, ], 
                                             df_nuisance = df[-partition, ],
                                             a = a, 
                                             a_prime = a_prime, bdw, cc)
  }
  
  NDE_psi_hat <- mean(ndenie_hat_folds[,1])
  NDE_v_hat <- mean(ndenie_hat_folds[,2])
  
  NDE_confidence_interval_lb <- NDE_psi_hat - qnorm(0.975) * sqrt(NDE_v_hat/n_rows)
  NDE_confidence_interval_ub <- NDE_psi_hat + qnorm(0.975) * sqrt(NDE_v_hat/n_rows)
  
  NIE_psi_hat <- mean(ndenie_hat_folds[,3])
  NIE_v_hat <- mean(ndenie_hat_folds[,4])
  
  NIE_confidence_interval_lb <- NIE_psi_hat - qnorm(0.975) * sqrt(NIE_v_hat/n_rows)
  NIE_confidence_interval_ub <- NIE_psi_hat + qnorm(0.975) * sqrt(NIE_v_hat/n_rows)
  
  return(list(NDE = c(NDE_confidence_interval_lb, NDE_confidence_interval_ub, NDE_v_hat, NDE_psi_hat),
              NIE = c(NIE_confidence_interval_lb, NIE_confidence_interval_ub, NIE_v_hat, NIE_psi_hat)) )
}

