################################################################
##  Univariate Properties for standardized data               ##
##                    J Di 02/12/2018                         ##
################################################################

#1. moments
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

plot(Y2nd,main = "Ex^2",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y3rd,main = "Ex^3",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y4th,main = "Ex^4",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()


#2. center moments
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = scale(hr50[,-1],center = T, scale = F)
moment_p = function(x,p){
  return(mean(x^p))
}
Y2nd = apply(Y,2,moment_p,p = 2)
Y3rd = apply(Y,2,moment_p,p = 3)
Y4th = apply(Y,2,moment_p,p = 4)

pdf("Write Up/plots/univariate_centermoments.pdf",width = 10,height = 10)
par(mfrow = c(3,1))

plot(Y2nd,main = "E(x-mu)^2",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y3rd,main = "E(x-mu)^3",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(Y4th,main = "E(x-mu)^4",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()

#3. cumulant
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = scale(hr50[,-1],center = T, scale = F)
moment_p = function(x,p){
  return(mean(x^p))
}
cumulant_4 = function(x){
  a = moment_p(x,4) - 3*moment_p(x,2)^2
  return(a)
}

k2 = apply(Y,2,moment_p,p = 2)
k3 = apply(Y,2,moment_p,p = 3)
k4 = apply(Y,2,cumulant_4)

pdf("Write Up/plots/univariate_cumulant.pdf",width = 10,height = 10)
par(mfrow = c(3,1))

plot(k2,main = "K2",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(k3,main = "K3",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(k4,main = "K4",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()


#4. sd moment (no center)
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = scale(hr50[,-1],center = F, scale = T)
sd_moment_p = function(x,p){
  return(mean(x^p))
}


k2 = apply(Y,2,sd_moment_p,p = 2)
k3 = apply(Y,2,sd_moment_p,p = 3)
k4 = apply(Y,2,sd_moment_p, p = 4)

pdf("Write Up/plots/univariate_sdmoment.pdf",width = 10,height = 10)
par(mfrow = c(3,1))

plot(k2,main = "sdM2",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(k3,main = "sdM3",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(k4,main = "sdM4",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()
