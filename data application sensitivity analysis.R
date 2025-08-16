# sensitivity to bandwidth and clipping (cc)
# compare a = 40 vs a' = 1500

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

if(0){# delete
  tmp = matrix(c(apply(JCdata, 2, function(x) mean(x, na.rm=T)), apply(JCdata, 2, function(x) max(x, na.rm=T)) ), ncol=2)
  rownames(tmp) = colnames(JCdata)
  View(tmp)
  
  xdat = JCdata[, c(4:68, 2,3,71)]
  misind = c(rep(0,5), 1,2, 
             rep(0,11),1,2, 
             0, 1,2,1,2,1,2,1,2,
             rep(0,4),
             1,2 ,0,0,
             1,2,1,2,0,1,2,1,2,
             0,0,1,2,1,2,
             0,1,2,1,2,
             0,0,1,2,0,
             0,0,0, 
             0,0,0)
  xname = c("female", "age", "white", "black", "Hispanic", "years of education", "GED diploma", "high school diploma",
    "native English", "divorced", "separated", "cohabiting", "married", "has children", "ever worked",
    "average weekly earnings in USD", "is household head", "household size", 
    "designated for nonresidential slot", "total household gross income", "total personal gross income", "mum's years of education", "dad's years of education",
    "dad did not work at 14", "received AFDC per month", "received public assistance per month", "received food stamps",
    "welfare receipt during childhood", "poor/fair health", "physical/emotional problems",
    "extent of marijuana use", "extent of hallucinogen use", "ever used other illegal drugs", "extent of smoking", "extent of alcohol consumption",
    "ever arrested", "times in prison", "time spent by Job Corps recruiter", "extent of recruiter support",
    "idea about wished training", "expected hourly wage after training", "expected improvement in maths",
    "expected improvement in reading skills","expected improvement in reading skills", "expected to be training for a job", "worried about training",
    "1st contact with recruiter by phone","1st contact with recruiter in office", "expected stay in training",
    "total training hours in year 1", "proportion of weeks employed in year 2", "any arrests in year 4")
  xvar = c("female", "age", "race_white", "race_black", "race_hispanic", "hgrd", "educ_geddiploma", "educ_hsdiploma",
    "ntv_engl", "marstat_divorced", "marstat_separated", "marstat_livetogunm",  "marstat_married",  "haschldY0",  "everwkd",
    "mwearn", "hohhd0", "peopleathome", 
    "nonres", "g10", "g12","hgrd_mum", "hgrd_dad",
    "work_dad_didnotwork", "g2", "g5", "g7",
    "welfare_child","h1_fair_poor", "h2",
    "h10", "h25", "h29", "h5", "h7",
    "i1", "i10", "e12", "e16",
    "e21", "e24usd", "e30",
    "e31", "e32", "e35", "e37",
    "e6_byphone", "e8_recruitersoffice", "e9ef", 
    "d" , "m", "y.binary")
   
   tab = matrix("", nrow = length(xname), ncol = 9)
   colnames(tab) = c("Missing", "Median (IQR)","0","1", "2", "3", "4", "5", "6")
   rownames(tab) = xname
   cc = 1
   for(j in 1:length(xname)){
     xn = xname[j]
     xind = which(colnames(xdat) == xvar[j])
     x = xdat[ , xind]
     
     if(misind[xind] == 1){ # exists missing
       xmisdum = xdat[, xind+1]
       obx = x[xmisdum == 0] # observed x
       ux = unique(obx) #unique x
       tab[cc, 1] = round(mean(xmisdum) * 100, 2) # missing proportion
     } 
     if(misind[xind] == 0){ # no missing
       obx = x
       ux = unique(x)
     }
     if(length(ux) == 2){ # binary
       tab[cc, 4] = round(mean(obx)*100, 2)
     } else if(length(ux) > 7){ # numeric
       qx = round(quantile(obx, c(0.5,0.25, 0.75)),2)
       tab[cc, 2] = paste0(qx[1], " (", qx[2], "-", qx[3], ")" )
     } else { # categorical
       tt = round(table(obx)*100/length(obx) ,2 )
       ttn = dimnames(tt)$obx
       for(k in 1:length(tt)){
         tab[cc, which(colnames(tab) == ttn[k])] = tt[k]
       }
     }
     cc = cc + 1
     
   }
   noquote(tab)
}

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
hajek_propensity_score <- function(df, model, a, cc){
  mu <- predict(model, newdata = df, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  bw <- density(df$d.log, bw = "nrd0")$bw
  hajek_kernel <- (1/bw) * dnorm((df$d.log - a)/bw, mean = 0, sd = 1)
  normalizing_constant <- mean(hajek_kernel/ps)
  normalized_ps <- ps * normalizing_constant
  
  if(1){ # clipping
    #cc = 0.01
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

hajek_propensity_score_corrected <- function(df_inference, df_nuisance, model, a, cc){
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
propensity_score <- function(df, model, a, cc){
  mu <- predict(model, newdata = df, type = "response")
  sigma <- sqrt(var(residuals(model)))
  ps <- dnorm(a, mu, sigma)
  
  #cc = 0.01
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


estimate_nde <- function(df_inference, df_nuisance, a, a_prime,  bdw, cc) {
  # Calculate bandwidth using Silverman's rule of thumb 
  bw_silverman <- bdw#0.5
  #bw_silverman <- density(df_inference$d.log, bw = "nrd0")$bw
  #delta = 0.25 
  #bw_undersmooth <- density(df_inference$d.log, bw = "nrd0")$bw * nrow(df_inference)^(-1 * delta/5)
  
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
  ps_a <- hajek_propensity_score(df = df_inference, model = propensity_model, a = a, cc)
  ps_a_prime <- hajek_propensity_score(df = df_inference, model = propensity_model,
                                       a = a_prime, cc)
  ps_a <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a, cc)
  ps_a_prime <- hajek_propensity_score_corrected(df_inference, df_nuisance, model = propensity_model, a = a_prime, cc)
  #ps_a <- propensity_score(df = df_inference, model = propensity_model, a = a, cc)
  #ps_a_prime <- propensity_score(df = df_inference, model = propensity_model, a = a_prime, cc)
  
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
crossfit <- function(df, folds, a, a_prime, bdw, cc) {
  n_rows <- nrow(df)
  row_partitions <- split(1:n_rows, cut(seq_along(1:n_rows), folds, 
                                        labels = FALSE))
  nde_hat_folds <- matrix(NA, nrow = folds, ncol = 2)
  for (i in 1:folds) {
    partition <- row_partitions[[i]]
    nde_hat <- estimate_nde(df_inference = df[partition, ], 
                            df_nuisance = df[-partition, ],
                            a = a, 
                            a_prime = a_prime, bdw, cc)
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
#d_list <- seq(100, 2000, 100)
d_list <- 1500
d_prime <- 40

bandwidth_vec <- seq(0.1, 0.8, 0.05)
clipping_vec <- seq(0.005, 0.1, 0.005)
book <- cbind(rep(bandwidth_vec, each = length(clipping_vec)), rep(clipping_vec, length(bandwidth_vec)))
colnames(book) = c("bandwidth", "clipping_threshold")

res <- matrix(0, nrow = nrow(book), ncol = 4)
for (i in 1:nrow(book)) {
  a <- log(d_list)#bcPower(d_list[i], optimal_lambda)
  a_prime <- log(d_prime) #bcPower(d_prime, optimal_lambda)
  res[i,] <- crossfit(JCdata, folds, a, a_prime, book[i,1], book[i,2])
}

#save(book, res, file = "CM_application_sensitivity.RData") # screen 12
#save(book, res, file = "CM_application_sensitivity_NotHajek.RData") # screen 81
save(book, res, file = "CM_application_sensitivity_Hajek.RData") # screen 11

################################################################################
####################### Visualize           ####################################
################################################################################
type = "_NotHajek"
load(paste0("/Users/yizhenxu/Downloads/CM_application_sensitivity",type,".RData"))
test <- as.data.frame(book)
test$Mean <- round(res[,4], 2)
test$SD = round(sqrt(res[,3]/nrow(JCdata)), 2)  

test = test[test$bandwidth <= 0.6 & test$clipping_threshold <= 0.05,]

### black white contour plots for paper
#par(mfrow = c(2,1), mai = c(1, 1, 0.1, 0.1))
data.loess <- loess(Mean ~ bandwidth*clipping_threshold, data = test)
pmean <-  predict(data.loess, newdata = test)
#MeanMat = matrix(test$Mean, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
MeanMat = matrix(pmean, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
pdf(paste0("/Users/yizhenxu/Downloads/CM_application_Mean",type,".pdf"), width = 5.5, height = 4)
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(MeanMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold",
        labcex = 1, lwd = 2, lty = 1) 
dev.off()

data.loess <- loess(SD ~ bandwidth*clipping_threshold, data = test)
psd <-  predict(data.loess, newdata = test)
SDMat = matrix(psd, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
pdf(paste0("/Users/yizhenxu/Downloads/CM_application_SD",type,".pdf"), width = 5.5, height = 4)
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(SDMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold",
        labcex = 1, lwd = 2, lty = 1) 
dev.off()



### color heatplot
library(ggplot2)
library(gridExtra)

p1 = ggplot(test, aes(x = bandwidth, y = clipping_threshold, fill = Mean)) +
  geom_tile(color = "black") +coord_equal(ratio = 5)+
  #geom_text(aes(label = Mean), color = "black", size = 3) +
  scale_fill_gradient2(limits = range(test$Mean),
                       low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000")    +
  theme(text = element_text(size = 10)) +ylab("Clipping Threshold")+xlab("Bandwidth")

p2 = ggplot(test, aes(x = bandwidth, y = clipping_threshold, fill = SD)) +
  geom_tile(color = "black") +coord_equal(ratio = 5)+
  #geom_text(aes(label = SD), color = "black", size = 3) +
  scale_fill_gradient2(limits = range(test$SD),
                       low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000")    +
  theme(text = element_text(size = 10))+ylab("Clipping Threshold")+xlab("Bandwidth")

ggsave(paste0("/Users/yizhenxu/Downloads/full", type, ".pdf"), arrangeGrob(p1, p2))

#png("/Users/yizhenxu/Downloads/foo.png",  width = 800, height=800)
#grid.arrange(p1, p2, nrow=2)
#dev.off()
