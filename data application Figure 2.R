# Hajek

################################################################################
####################### Main Inference Code ####################################
################################################################################

set.seed(2025)
folds <- 3
d_list <- seq(100, 2000, 100)
d_prime <- 60

resNDE = resNIE = matrix(0, nrow = length(d_list), ncol = 4)
nn = nrow(JCdata)
for (i in 1:length(d_list)) {
  a <- log(d_list[i])
  a_prime <- log(d_prime)
  tmp = crossfit(JCdata, folds, a, a_prime)
  resNDE[i,] <- tmp[[1]]
  resNIE[i,] <- tmp[[2]]
}
save(resNDE, resNIE, file = "CM_application_p1_Hajek_60.RData") # screen 2161 

# clipped away < 5%

load("~/Downloads/CM_application_p1_Hajek_60.RData")

pdf("/Users/yizhenxu/Downloads/CM_application_Hajek_NDE_NIE.pdf", width = 10, height = 4)

par(mfrow = c(1,2))
plot(d_list, resNDE[,4], type = "l", ylab = "Natural Direct Effect", xlab = "Treatment a", ylim =c(-0.11,0.11))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(d_list,rev(d_list)),c(resNDE[,1],rev(resNDE[,2])),col = "grey75", border = FALSE)
lines(d_list, resNDE[,4], lwd = 2)
abline(h=0)

plot(d_list, resNIE[,4], type = "l", ylab = "Natural Indirect Effect", xlab = "Treatment a", ylim =c(-0.011,0.011))
polygon(c(d_list,rev(d_list)),c(resNIE[,1],rev(resNIE[,2])),col = "grey75", border = FALSE)
lines(d_list, resNIE[,4], lwd = 2)
abline(h=0)
dev.off()