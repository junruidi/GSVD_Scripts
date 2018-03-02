################################################################
##  Univariate Properties for standardized data               ##
##                    J Di 02/12/2018                         ##
################################################################

#1. univariate properties
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = scale(hr50[,-1],center = F, scale = F)
moment_p = function(x,p){
  return(mean(x^p))
}
Y2nd = apply(Y,2,moment_p,p = 2)
Y3rd = apply(Y,2,moment_p,p = 3)
Y4th = apply(Y,2,moment_p,p = 4)

pdf("Write Up/plots/univariate_moments.pdf",width = 10,height = 10)
par(mfrow = c(3,1))

plot(Y2nd,main = "Variance",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y3rd,main = "Skewness",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y4th,main = "Kurtosis",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()

