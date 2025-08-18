# Hajek

################################################################################
####################### Main Inference Code ####################################
################################################################################

folds <- 3
d_list <- seq(100, 2000, 100)
d_prime <- 40

resNDE = resNIE = matrix(0, nrow = length(d_list), ncol = 4)
for (i in 1:length(d_list)) {
  a <- log(d_list[i])
  a_prime <- log(d_prime)
  tmp = crossfit(JCdata, folds, a, a_prime)
  resNDE[i,] <- tmp[[1]]
  resNIE[i,] <- tmp[[2]]
}

save(resNDE, resNIE, file = "CM_application_p1_Hajek.RData") # screen 30 3 folds
# clipped away 0.1% and 0% for a, 4.7%, 5%, 9.7% for a_prime


load("~/Downloads/CM_application_p1_Hajek.RData")

pdf("/Users/yizhenxu/Downloads/CM_application_Hajek_NDE_NIE.pdf", width = 10, height = 4)

par(mfrow = c(1,2))
#plot(d_list, resNDE[,4], ylim = range(resNDE[,c(1,2,4)]), type = "l", ylab = "NDE", main = "clipping at 0.02, 3 folds, undersmooth")
plot(d_list, resNDE[,4], ylim = range(resNDE[,c(1,2,4)]), type = "l", ylab = "Natural Direct Effect", xlab = "Treatment a")
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(d_list,rev(d_list)),c(resNDE[,1],rev(resNDE[,2])),col = "grey75", border = FALSE)
lines(d_list, resNDE[,4], lwd = 2)
#add red lines on borders of polygon
lines(d_list, resNDE[,2], col="red",lty=2)
lines(d_list, resNDE[,1], col="red",lty=2)
abline(h=0)

#plot(d_list, resNIE[,4], ylim = range(resNIE[,c(1,2,4)]), type = "l", ylab = "NDE", main = "clipping at 0.02, 3 folds, undersmooth")
plot(d_list, resNIE[,4], ylim = range(resNIE[,c(1,2,4)]), type = "l", ylab = "Natural Indirect Effect", xlab = "Treatment a")
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(d_list,rev(d_list)),c(resNIE[,1],rev(resNIE[,2])),col = "grey75", border = FALSE)
lines(d_list, resNIE[,4], lwd = 2)
#add red lines on borders of polygon
lines(d_list, resNIE[,2], col="red",lty=2)
lines(d_list, resNIE[,1], col="red",lty=2)
abline(h=0)
dev.off()