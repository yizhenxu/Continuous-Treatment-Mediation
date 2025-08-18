
# Table 1: Compare our estimator to mean eta and Huber's, under different sample sizes n = 500, 1000, 2000.
#          Use Silverman bandwidth, under different mis-specification settings, report average mean bias and RMSE

# Table 2: sensitivity analysis over sample sizes and bandwidths, 
#          report absolute average bias, mean of $\sqrt{\hat{V}}$, and coverage


library(MASS)
library(cubature)
library(mvtnorm)
library(betareg)
set.seed(0)

sigmoid <- function(x){
  return(exp(x)/(1 + exp(x)))
}

logit <- function(x){
  return(log(x/(1 - x)))
}



# Set simulation parameters
a = 4.5
a_prime = 6
mean_X <- c(0, 0, 0, 0)
sigma_X <- diag(c(0.25, 0.1,0.8, 0.5)) # diag(c(1, 2, 1.5, 0.5))

y1 = -1 
y2 = 5

# Monte Carlo integration to calculate true value of parameter
x <- mvrnorm(100000, mean_X, sigma_X)

logistic_M <- sigmoid(-5 + 5*a_prime + 2*x[,2] + 10*a_prime*x[,3])
result_m1 <- (y1 * a +20+ y2 * 1 * x[,1] + x[,2] ) * logistic_M
result_m0 <- (y1 * a +0+ y2 * 0 * x[,1] + x[,2] ) * (1 - logistic_M)
true_val_mc <- mean(result_m1 + result_m0)

simulate <- function(n) {
  # Simulate  covariates
  X <- mvrnorm(n, mean_X, sigma_X)
  
  # Simulate treatment variable
  A =  5 + X[, 1] + 0.2*X[,1]^2 + rnorm(n)
  
  # Simulate mediator
  p_M <- sigmoid(-5+ 5*A + 2*X[,2] + 10*A*X[,3])
  M <- rbinom(n, size = 1, p = p_M)
  
  # Simulate outcome variable
  Y <- y1 * A +20*M+ y2 * M * X[, 1] + X[, 2]  + rnorm(n)
  df <- data.frame(Y, A, M, X)
  return(df)
}


################################################################################
##################### Inference Calculators ####################################
################################################################################
outcome_regression <- function(df_inference, OR_model, eval_a){
  df_inference_a <- df_inference
  df_inference_a$A <- eval_a
  y_predictions <- predict(OR_model, newdata = df_inference_a, type = "response")
  return(y_predictions)
}

# Propensity score with parametric density estimation
propensity_score <- function(df_inference, df_nuisance, propensity_model, eval_a){
  df_inference_a <- df_inference
  df_inference_a$A <- eval_a
  ps <- dnorm(a, predict(propensity_model, newdata = df_inference_a, type = "response"), summary(propensity_model)$sigma)
  
  return(1/ps)
}

# Density ratio
density_ratio <- function(df_inference, mediator_model, a, a_prime){
  df_a_prime <- df_inference
  df_a_prime["A"] <- a_prime
  df_a <- df_inference
  df_a["A"] <- a
  density_a <- predict(mediator_model, newdata = df_a, type = "response")
  density_a_prime <- predict(mediator_model, newdata = df_a_prime, type = "response")
  density_ratio <- (df_a$M == 1) * density_a_prime/density_a + 
    (df_a$M == 0) * (1 - density_a_prime)/(1 - density_a)
  return(density_ratio)
}

eta_hat_fn <- function(df_inference, mediator_model, OR_model, a, a_prime){
  df_a_prime <- df_inference
  df_a_prime["A"] <- a_prime
  density_a_prime = predict(mediator_model, newdata = df_a_prime, type = "response")
  
  df_inference_a <- df_inference
  df_inference_a$A <- a
  df_inference_a$M <- 1
  y_predictions <- predict(OR_model, newdata = df_inference_a, type = "response")  
  eta_hat <- y_predictions * density_a_prime 
  
  df_inference_a$M <- 0
  y_predictions <- predict(OR_model, newdata = df_inference_a, type = "response")
  eta_hat <- eta_hat + y_predictions * (1 - density_a_prime) 
  
  return(eta_hat)
}

################################################################################
##################### Estimators ###############################################
################################################################################

estimate_psi_if_yma <- function(df_inference, df_nuisance, bandwidth_silverman) {
  # Calculate bandwidth using Silverman's rule of thumb 
  # bandwidth_silverman <- density(df_inference$A, bw = "nrd0")$bw
  
  OR_model <- glm(Y ~ A + M + M:X1 + X2 - 1, family = gaussian, data = df_nuisance)
  
  propensity_model <- lm(A ~ X1 + I(X1^2) , data = df_nuisance)

  mediator_model <- glm(M ~ A + X2 + A:X3 , 
                        data = df_nuisance, family = binomial)

  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  propensity_prediction_a <- propensity_score(df_inference, df_nuisance, propensity_model,  a)
  propensity_prediction_a_prime <- propensity_score(df_inference, df_nuisance, propensity_model,  a_prime)
  densratio_prediction <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  if_kernel_1 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_1 = if_kernel_1 * propensity_prediction_a * 
    densratio_prediction * (df_inference$Y - OR_prediction)
  
  if_kernel_2 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a_prime)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_2 = if_kernel_2 * propensity_prediction_a_prime *
    (OR_prediction - eta)
  
  psi_hat <- mean(if_term_1 + if_term_2 + eta)
  v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  
  ww = if_kernel_1 * propensity_prediction_a * densratio_prediction
  hest = sum(ww * df_inference$Y) / sum(ww)
  
  return(c(psi_hat, v_hat, mean(eta), hest))
}

estimate_psi_if_ym <- function(df_inference, df_nuisance, bandwidth_silverman) {
  # Calculate bandwidth using Silverman's rule of thumb
  # bandwidth_silverman <- density(df_inference$A, bw = "nrd0")$bw
  
  OR_model <- glm(Y ~ A + M + M:X1 + X2 - 1, family = gaussian, data = df_nuisance)
  
  propensity_model <- lm(A ~ 1, data = df_nuisance)#~1
  
  mediator_model <- glm(M ~ A + X2 + A:X3 , 
                        data = df_nuisance, family = binomial)
  
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  propensity_prediction_a <- propensity_score(df_inference, df_nuisance, propensity_model,  a)
  propensity_prediction_a_prime <- propensity_score(df_inference, df_nuisance, propensity_model,  a_prime)
  densratio_prediction <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  if_kernel_1 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_1 = if_kernel_1 * propensity_prediction_a * densratio_prediction *
    (df_inference$Y - OR_prediction)
  
  if_kernel_2 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a_prime)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_2 = if_kernel_2 * propensity_prediction_a_prime * 
    (OR_prediction - eta)
  
  psi_hat <- mean(if_term_1 + if_term_2 + eta)
  v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  
  ww = if_kernel_1 * propensity_prediction_a * densratio_prediction
  hest = sum(ww * df_inference$Y) / sum(ww)
  
  return(c(psi_hat, v_hat, mean(eta), hest))

}

estimate_psi_if_ma <- function(df_inference, df_nuisance, bandwidth_silverman) {
  # Calculate bandwidth using Silverman's rule of thumb
  # bandwidth_silverman <- density(df_inference$A, bw = "nrd0")$bw
  
  OR_model <- glm(Y ~ M + X1 +X2 , family = gaussian, data = df_nuisance) # M + X1 -1
  
  propensity_model <- lm(A ~ X1 + I(X1^2) , data = df_nuisance)
  
  mediator_model <- glm(M ~ A + X2 + A:X3 , 
                        data = df_nuisance, family = binomial)
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  propensity_prediction_a <- propensity_score(df_inference, df_nuisance, propensity_model,  a)
  propensity_prediction_a_prime <- propensity_score(df_inference, df_nuisance, propensity_model, a_prime)
  densratio_prediction <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  if_kernel_1 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_1 = if_kernel_1 * propensity_prediction_a * densratio_prediction * 
    (df_inference$Y - OR_prediction)
  
  if_kernel_2 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a_prime)/bandwidth_silverman, 
          mean = 0, sd = 1)
  if_term_2 = if_kernel_2 * propensity_prediction_a_prime * (OR_prediction - eta)
  
  
  psi_hat <- mean(if_term_1 + if_term_2 + eta)
  v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  
  ww = if_kernel_1 * propensity_prediction_a * densratio_prediction
  hest = sum(ww * df_inference$Y) / sum(ww)
  
  return(c(psi_hat, v_hat, mean(eta), hest))
}

estimate_psi_if_ya <- function(df_inference, df_nuisance, bandwidth_silverman) {
  # Calculate bandwidth using Silverman's rule of thumb
  # bandwidth_silverman <- density(df_inference$A, bw = "nrd0")$bw
  
  OR_model <- glm(Y ~ A + M + M:X1 + X2 - 1, family = gaussian, data = df_nuisance)
  
  propensity_model <- lm(A ~ X1 + I(X1^2) , data = df_nuisance)
  
  mediator_model <- glm(M ~  1, 
                        data = df_nuisance, family = binomial) #A + X1 + X3
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  propensity_prediction_a <- propensity_score(df_inference, df_nuisance, propensity_model,  a)
  propensity_prediction_a_prime <- propensity_score(df_inference, df_nuisance, propensity_model, a_prime)
  densratio_prediction <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  #OR_prediction = eta = 0
  if_kernel_1 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_1 = if_kernel_1 * propensity_prediction_a * 
    densratio_prediction * (df_inference$Y - OR_prediction)
  
  if_kernel_2 = (1/bandwidth_silverman) * dnorm((df_inference$A - a_prime)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_2 = if_kernel_2 * propensity_prediction_a_prime * (OR_prediction - eta)
  
  psi_hat <- mean(if_term_1 + if_term_2 + eta)
  v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  
  ww = if_kernel_1 * propensity_prediction_a * densratio_prediction
  hest = sum(ww * df_inference$Y) / sum(ww)
  
  return(c(psi_hat, v_hat, mean(eta), hest))
}

estimate_psi_if_incorrect <- function(df_inference, df_nuisance, bandwidth_silverman) {
  # Calculate bandwidth using Silverman's rule of thumb
  # bandwidth_silverman <- density(df_inference$A, bw = "nrd0")$bw
  
  OR_model <- glm(Y ~ A  +X2, family = gaussian, data = df_nuisance)
  
  propensity_model <- lm(A ~ 1, data = df_nuisance)

  mediator_model <- glm(M ~ A +A:X3 , 
                        data = df_nuisance, family = binomial)
  
  OR_prediction <- outcome_regression(df_inference, OR_model, eval_a = a)
  propensity_prediction_a <- propensity_score(df_inference, df_nuisance, propensity_model,  a)
  propensity_prediction_a_prime <- propensity_score(df_inference, df_nuisance, propensity_model, a_prime)
  densratio_prediction <- density_ratio(df_inference, mediator_model, a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df_inference, mediator_model, OR_model, a = a, a_prime = a_prime)
  
  if_kernel_1 = (1/bandwidth_silverman) * dnorm((df_inference$A - a)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_1 = if_kernel_1 * propensity_prediction_a * densratio_prediction * 
    (df_inference$Y - OR_prediction)
  
  if_kernel_2 = (1/bandwidth_silverman) * 
    dnorm((df_inference$A - a_prime)/bandwidth_silverman, mean = 0, sd = 1)
  if_term_2 = if_kernel_2 * propensity_prediction_a_prime * (OR_prediction - eta)
  
  psi_hat <- mean(if_term_1 + if_term_2 + eta)
  v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  
  ww = if_kernel_1 * propensity_prediction_a * densratio_prediction
  hest = sum(ww * df_inference$Y) / sum(ww)
  
  return(c(psi_hat, v_hat, mean(eta), hest))
}


# Crossfitting sceheme
crossfit <- function(df, folds) {
  n_rows <- nrow(df)
  
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds, labels = FALSE))
  psi_estimate_cf_yma <- matrix(NA, nrow = folds, ncol = 4)
  psi_estimate_cf_ym <- matrix(NA, nrow = folds, ncol = 4)
  psi_estimate_cf_ma <- matrix(NA, nrow = folds, ncol = 4)
  psi_estimate_cf_ya <- matrix(NA, nrow = folds, ncol = 4)
  psi_estimate_cf_incorrect <- matrix(NA, nrow = folds, ncol = 4)
  
  for (i in 1:folds) {
    parition <- row_partitions[[i]]
    bw <- density(df[parition, ]$A, bw = "nrd0")$bw # bandwidth from df_inference (I_\ell)
    
    psi_hat_yma <- estimate_psi_if_yma(df[parition, ], df[-parition, ], bw)
    psi_estimate_cf_yma[i,] <- psi_hat_yma
    
    psi_hat_ym <- estimate_psi_if_ym(df[parition, ], df[-parition, ], bw)
    psi_estimate_cf_ym[i,] <- psi_hat_ym
    
    psi_hat_ma <- estimate_psi_if_ma(df[parition, ], df[-parition, ], bw)
    psi_estimate_cf_ma[i,] <- psi_hat_ma
    
    psi_hat_ya <- estimate_psi_if_ya(df[parition, ], df[-parition, ], bw)
    psi_estimate_cf_ya[i,] <- psi_hat_ya
    
    psi_hat_incorrect <- estimate_psi_if_incorrect(df[parition, ], df[-parition, ], bw)
    psi_estimate_cf_incorrect[i,] <- psi_hat_incorrect
  }
  psi_estimates <- c(mean(psi_estimate_cf_yma[,1]), mean(psi_estimate_cf_ym[,1]),
                     mean(psi_estimate_cf_ma[,1]), mean(psi_estimate_cf_ya[,1]),
                     mean(psi_estimate_cf_incorrect[,1]))
  var_estimates <- c(mean(psi_estimate_cf_yma[,2]), mean(psi_estimate_cf_ym[,2]),
                     mean(psi_estimate_cf_ma[,2]), mean(psi_estimate_cf_ya[,2]),
                     mean(psi_estimate_cf_incorrect[,2]))
  eta_estimates <- c(mean(psi_estimate_cf_yma[,3]), mean(psi_estimate_cf_ym[,3]),
                     mean(psi_estimate_cf_ma[,3]), mean(psi_estimate_cf_ya[,3]),
                     mean(psi_estimate_cf_incorrect[,3]))
  huber_estimates <- c(mean(psi_estimate_cf_yma[,4]), mean(psi_estimate_cf_ym[,4]),
                     mean(psi_estimate_cf_ma[,4]), mean(psi_estimate_cf_ya[,4]),
                     mean(psi_estimate_cf_incorrect[,4]))
  
  confidence_interval_lb <- psi_estimates - qnorm(0.975) * sqrt(var_estimates/n_rows)
  confidence_interval_ub <- psi_estimates + qnorm(0.975) * sqrt(var_estimates/n_rows)
  
  cp <- 1 * ( (true_val_mc > confidence_interval_lb) & (true_val_mc < confidence_interval_ub) )
  
  return(c(psi_estimates, var_estimates, cp, eta_estimates, huber_estimates ))
}

### Perform 1000 Monte Carlo simulations
n_runs <- 1000
folds <- 3

ssvec <- c(2000, 5000, 8000)

for(sample_size in ssvec){
  print(sample_size)
  
  psi_hat_list <- matrix(NA, nrow = n_runs, ncol = length(c("psi", "var", "CP", "eta","huber"))*5)
  colnames(psi_hat_list) = c(paste0("psi_", c("YMA", "YM", "MA", "YA", "Incorr")),
                             paste0("var_", c("YMA", "YM", "MA", "YA", "Incorr")),
                             paste0("CP_", c("YMA", "YM", "MA", "YA", "Incorr")),
                             paste0("eta_", c("YMA", "YM", "MA", "YA", "Incorr")),
                             paste0("huber_", c("YMA", "YM", "MA", "YA", "Incorr")))

  for(i in 1:n_runs){
    print(i)
    set.seed(i)
    df <- simulate(sample_size)
    psi_hat_list[i,] <- crossfit(df, folds)
    
  }

  save(psi_hat_list, file = paste0("table1_M1_",sample_size,".RData"))
  
  
}


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

rdnum = 2 # round results to 2 decimal 
tabshow = c()
for(sample_size in ssvec){
  load(paste0("/Users/yizhenxu/table1_M1_",sample_size,".RData"))
  tab = matrix(NA, ncol = 5, nrow = 3)
  rownames(tab) = c("Proposal", "CV Eta", "CV Huber")
  colnames(tab) = c("YMA", "YM", "MA", "YA", "None") # correct model is -
  psi_est = psi_hat_list[,1:5]
  tab[1,] = paste0(round(abs(apply(psi_est - true_val_mc, 2, mean)),rdnum), " (",
                   round(sqrt(apply(psi_est - true_val_mc, 2, function(x) mean(x^2))),rdnum), ")")
  psi_est = psi_hat_list[,16:20]
  tab[2,] = paste0(round(abs(apply(psi_est - true_val_mc, 2, mean)),rdnum), " (",
                   round(sqrt(apply(psi_est - true_val_mc, 2, function(x) mean(x^2))),rdnum), ")")
  psi_est = psi_hat_list[,21:25]
  tab[3,] = paste0(round(abs(apply(psi_est - true_val_mc, 2, mean)),rdnum), " (",
                   round(sqrt(apply(psi_est - true_val_mc, 2, function(x) mean(x^2))),rdnum), ")")
  tabshow = rbind(tabshow, tab)
}
tabshow = cbind(c(ssvec[1], rep("",2),ssvec[2], rep("",2),ssvec[3], rep("",2)), tabshow)
noquote(tabshow)
