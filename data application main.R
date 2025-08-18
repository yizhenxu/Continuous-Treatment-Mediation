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
setwd("/users/yxu2/delete")
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

source("data applciation functions.R")

################################################################################
##################### Data Summary #############################################
################################################################################
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
################################################################################
############################ Data Analysis #####################################
################################################################################

source("data application Figure 2.R")

################################################################################
##################### Sensitivity Analysis #####################################
################################################################################
# sensitivity to bandwidth and clipping (cc)
# compare a' = 40 vs a = 1500

source("data application Figure 3.R")
