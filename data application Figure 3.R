

################################################################################
####################### Main Inference Code ####################################
################################################################################

folds <- 3
#d_list <- seq(100, 2000, 100)
d_list <- 1500
d_prime <- 60

bandwidth_vec <- seq(0.1, 0.8, 0.05)
clipping_vec <- seq(0.005, 0.1, 0.005)
book <- cbind(rep(bandwidth_vec, each = length(clipping_vec)), rep(clipping_vec, length(bandwidth_vec)))
colnames(book) = c("bandwidth", "clipping_threshold")

sresNDE = sresNIE = matrix(0, nrow = nrow(book), ncol = 4)
for (i in 1:nrow(book)) {
  a <- log(d_list)
  a_prime <- log(d_prime) #bcPower(d_prime, optimal_lambda)
  tmp = crossfit_sensitivity(JCdata, folds, a, a_prime, book[i,1], book[i,2])
  sresNDE[i,] <- tmp[[1]]
  sresNIE[i,] <- tmp[[2]]
}

save(book, sresNDE, sresNIE, file = "CM_application_sensitivity_Hajek_60.RData") # screen 2169

################################################################################
####################### Visualize           ####################################
################################################################################

load("~/Downloads/CM_application_sensitivity_Hajek_60.RData")

pdf(paste0("/Users/yizhenxu/Downloads/CM_application_Contours.pdf"), width = 11, height = 7)
par(mfrow = c(2,2))

type = "_NDE"
test <- as.data.frame(book)
test$Mean <- round(sresNDE[,4], 2)
test$SD = round(sqrt(sresNDE[,3]/nrow(JCdata)), 2)  

test = test[test$bandwidth <= 0.6 & test$clipping_threshold <= 0.05,]

### black white contour plots for paper
data.loess <- loess(Mean ~ bandwidth*clipping_threshold, data = test)
pmean <-  predict(data.loess, newdata = test)
MeanMat = matrix(pmean, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(MeanMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold", main = "(a) NDE Mean",
        labcex = 1, lwd = 2, lty = 1) 

data.loess <- loess(SD ~ bandwidth*clipping_threshold, data = test)
psd <-  predict(data.loess, newdata = test)
SDMat = matrix(psd, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(SDMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold", main = "(b) NDE SD",
        labcex = 1, lwd = 2, lty = 1) 

########

type = "_NIE"
test <- as.data.frame(book)
test$Mean <- round(sresNIE[,4], 4)
test$SD = round(sqrt(sresNIE[,3]/nrow(JCdata)), 4)  

test = test[test$bandwidth <= 0.6 & test$clipping_threshold <= 0.05,]

data.loess <- loess(Mean ~ bandwidth*clipping_threshold, data = test)
pmean <-  predict(data.loess, newdata = test)
MeanMat = matrix(pmean, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(MeanMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold", main = "(c) NIE Mean",
        labcex = 1, lwd = 2, lty = 1) 

data.loess <- loess(SD ~ bandwidth*clipping_threshold, data = test)
psd <-  predict(data.loess, newdata = test)
SDMat = matrix(psd, nrow = length(unique(test$clipping_threshold)), ncol = length(unique(test$bandwidth)))
contour(x = unique(test$bandwidth), y = unique(test$clipping_threshold), z = t(SDMat), 
        xlab = "Bandwidth", ylab = "Clipping Threshold", main = "(d) NIE SD",
        labcex = 1, lwd = 2, lty = 1) 
dev.off()

