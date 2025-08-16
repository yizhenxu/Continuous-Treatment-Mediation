# Hajek or no?

# install.packages("reshape2")
# install.packages("Hmisc")
rm(list=ls())
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(betareg)
library(car)
library(sn)
library(randomForest)
################################################################################
##################### Data Cleaning ############################################
################################################################################
# Load and subset data with positive d
#load("hhll/JCdata.RData")
load("JCdata.RData")
#load("/Users/yizhenxu/Downloads/hhll/JCdata.RData")
JCdata = JCdata[JCdata$d > 0,]

# Transform m for analysis
JCdata$m.shifted <- JCdata$m/100
JCdata$m.shifted <- (JCdata$m.shifted * (nrow(JCdata) - 1) + 0.5)/nrow(JCdata)

# Transform d for analysis 
JCdata$d.log <- log(JCdata$d)

# Binarize y - AUC with logistic regression model is 0.7663
JCdata$y.binary <- 1*(JCdata$y > 0)

# Apply box-cox trasnformation
# optimal_lambda <- boxcox(d ~ 1, data = JCdata)$x[which.max(
#                                                 boxcox(d ~ 1, data = JCdata)$y)]
# JCdata$d.log <- bcPower(JCdata$d, optimal_lambda)

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

hajek_propensity_score_corrected <- function(df_inference, df_nuisance, model, a){
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
  
  cc = 0.01 # clipping threshold
  print(mean(normalized_ps < cc))
  
  normalized_ps[normalized_ps < cc] <- cc
  return(1/normalized_ps)
  
}

# Propensity score with parametric density estimation
propensity_score <- function(df, model, a){
  mu <- predict(model, newdata = df, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  
  cc = 0.01
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


estimate_nde <- function(df_inference, df_nuisance, a, a_prime) {
  # Calculate bandwidth using Silverman's rule of thumb 
  #bw_silverman <- 0.5
  bw_silverman <- density(df_inference$d.log, bw = "nrd0")$bw
  #delta = 0.25 
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
  
  #ps_a <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a)
  #ps_a_prime <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a_prime)
  ps_a <- propensity_score(df = df_inference, model = propensity_model, a = a)
  ps_a_prime <- propensity_score(df = df_inference, model = propensity_model,
                                  a = a_prime)
  
  
  
  if(0){
    ps_a <- hajek_propensity_score(df = df_inference, model = propensity_model, a = a)
    ps_a_prime <- hajek_propensity_score(df = df_inference, model = propensity_model,
                                         a = a_prime)
    
    ps_a1 <- propensity_score(df = df_inference, model = propensity_model, a = a)
    ps_a_prime1 <- propensity_score(df = df_inference, model = propensity_model,
                                    a = a_prime)
    
    ps_a2 <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a)
    ps_a_prime2 <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a_prime)
    
    par(mfrow = c(1,3)); hist(ps_a, 30, main='Hajek');hist(ps_a2, 30, main='Hajek Corrected'); hist(ps_a1, 30, main='No Hajek')
    par(mfrow = c(1,3)); hist(ps_a_prime, 30, main='Hajek');hist(ps_a_prime2, 30, main='Hajek Corrected'); hist(ps_a_prime1, 30, main='No Hajek')
    
  }

  densratio_prediction <- density_ratio(df = df_inference,
                                        model = mediator_model,
                                        phi = mediator_phi,
                                        a = a, a_prime = a_prime)
  eta <- eta_hat_fn(df = df_inference, mediator_model = mediator_model,
                    mediator_phi = mediator_phi, OR_model = OR_model,
                    a = a, a_prime = a_prime)
  
  # ATE predictions
  OR_prediction_a_prime <- outcome_regression(df = df_inference, 
                                              model = OR_model, a = a_prime)
  eta_a_prime <- eta_hat_fn(df = df_inference, mediator_model = mediator_model,
                            mediator_phi = mediator_phi, OR_model = OR_model,
                            a = a_prime, a_prime = a_prime)
  
  # Throw away data with unstable propensity weights
  # indices <- which(ps_a_prime < 1000000)
  # OR_prediction <- OR_prediction[indices]
  # ps_a <- ps_a[indices]
  # ps_a_prime <- ps_a_prime[indices]
  # densratio_prediction <- densratio_prediction[indices]
  # eta <- eta[indices]
  # OR_prediction_a_prime <- OR_prediction_a_prime[indices]
  # eta_a_prime <- eta_a_prime[indices]
  # df_inference <- df_inference[indices,]
  
  
  # Mediator counterfactual
  if_kernel_1 = (1/bw_silverman) *
    dnorm((df_inference$d.log - a)/bw_silverman, mean = 0, sd = 1)
  
  if_term_1 = if_kernel_1 * ps_a *
    densratio_prediction * (df_inference$y.binary - OR_prediction)
  
  if_kernel_2 = (1/bw_silverman) *
    dnorm((df_inference$d.log - a_prime)/bw_silverman, mean = 0, sd = 1)
  
  if_term_2 = if_kernel_2 * ps_a_prime * (OR_prediction - eta)
  #mean(if_term_2)
  
  # ATE counterfactual
  if_kernel_3 = (1/bw_silverman) * 
    dnorm((df_inference$d.log - a_prime)/bw_silverman, mean = 0, sd = 1)
  
  if_term_3 = if_kernel_3 * ps_a_prime *
    (df_inference$y.binary - eta_a_prime)
  
  # if_kernel_4 = (1/bw_silverman) *
  #   dnorm((df_inference$d.log - a_prime)/bw_silverman, mean = 0, sd = 1)
  
  # if_term_4 = if_kernel_4 * ps_a_prime * (OR_prediction_a_prime - eta_a_prime)
  # mean(if_term_3 + if_term_4)
  # 
  psi_a_aprime = mean(if_term_1 + if_term_2 + eta)
  psi_aprime = mean(if_term_3 + eta_a_prime)
  nde = psi_a_aprime - psi_aprime
  
  vhat =  mean((if_term_1 + if_term_2 + eta -if_term_3 - eta_a_prime - nde)^2) #cancel *bw_silverman with denominator later
  #psi_hat <- mean(if_term_1 + if_term_2 + eta)
  #v_hat <- mean((if_term_1 + if_term_2 + eta - psi_hat)^2)
  return(c(nde, vhat))
}

# Cross fitting scheme
crossfit <- function(df, folds, a, a_prime) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds, 
                                        labels = FALSE))
  nde_hat_folds <- matrix(NA, nrow = folds, ncol = 2)
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    nde_hat <- estimate_nde(df_inference = df[partition, ], 
                            df_nuisance = df[-partition, ],
                            a = a, 
                            a_prime = a_prime)
    nde_hat_folds[i,] <- nde_hat
  }
  #mean(nde_hat_folds)
  psi_hat <- mean(nde_hat_folds[,1])
  v_hat <- mean(nde_hat_folds[,2])
  
  confidence_interval_lb <- psi_hat - qnorm(0.975) * sqrt(v_hat/n_rows)
  confidence_interval_ub <- psi_hat + qnorm(0.975) * sqrt(v_hat/n_rows)
  return(c(confidence_interval_lb, confidence_interval_ub, v_hat, psi_hat))
}


################################################################################
####################### Main Inference Code ####################################
################################################################################

folds <- 3
d_list <- seq(100, 2000, 100)
d_prime <- 40

res <- matrix(0, nrow = length(d_list), ncol = 4)
for (i in 1:length(d_list)) {
  a <- log(d_list[i])#bcPower(d_list[i], optimal_lambda)
  a_prime <- log(d_prime) #bcPower(d_prime, optimal_lambda)
  res[i,] <- crossfit(JCdata, folds, a, a_prime)
}

#save(res, file = "CM_application_p1_Hajek.RData") # screen 81 3 folds
# clipped away 0.1% and 0% for a, 4.7%, 5%, 9.7% for a_prime

save(res, file = "CM_application_p1_NotHajek.RData") # screen 12 3 folds
# clipped away 0.07% and 0% for a,  39%, 30%, 38% for a_prime



load("~/Downloads/CM_application_p1_Hajek.RData")
load("~/Downloads/CM_application_p1_NotHajek.RData")

pdf("/Users/yizhenxu/Downloads/CM_application_Hajek.pdf", width = 5.5, height = 4)
pdf("/Users/yizhenxu/Downloads/CM_application_NotHajek.pdf", width = 5.5, height = 4)

#plot(d_list, res[,4], ylim = range(res[,c(1,2,4)]), type = "l", ylab = "NDE", main = "clipping at 0.02, 3 folds, undersmooth")
plot(d_list, res[,4], ylim = range(res[,c(1,2,4)]), type = "l", ylab = "Natural Direct Effect", xlab = "Treatment a")
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(d_list,rev(d_list)),c(res[,1],rev(res[,2])),col = "grey75", border = FALSE)
lines(d_list, res[,4], lwd = 2)
#add red lines on borders of polygon
lines(d_list, res[,2], col="red",lty=2)
lines(d_list, res[,1], col="red",lty=2)
abline(h=0)
dev.off()