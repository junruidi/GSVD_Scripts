rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/RobPTI.R")
library(MASS)
Y = scale(hr50[,-1],center = T, scale = F)
u = svd(Y)$u
v = svd(Y)$v
s = diag(svd(Y)$d)
w = v %*% solve(s)

Yp = Y %*% w



order3_y = Rob_TSM(Y, L = 30, N = 30, order = 3)
phi3_y = order3_y$eigenv
ev3_y = order3_y$eigenl

order3_yp = Rob_TSM(Yp, L = 30, N = 30, order = 3)
phi3_yp = order3_yp$eigenv
ev3_yp = order3_yp$eigenl

transfer = ginv(t(w)) %*% phi3_yp
normalization = function(x){return(x/sqrt(sum(x^2)))}
norm_v = function(x){return(sqrt(sum(x^2)))}
transfer2 = apply(transfer, 2, normalization)
norm_values = apply(transfer, 2, norm_v)

ev3_yp * norm_values^3


TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)


pdf("phiy.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi3_y[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  
}
dev.off()

pdf("phiy_transfer.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(transfer2[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  
}
dev.off()

#robust power tensor iteration
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/RobPTI.R")
library(MASS)
Y = scale(hr50[,-1],center = T, scale = F)
u = svd(Y)$u
v = svd(Y)$v
s = diag(svd(Y)$d)
w = v %*% solve(s)

Yp = Y %*% w


#1. 2nd order PC
pca = prcomp(Y, center = F, scale. = F)
phi2 = pca$rotation
sc2 = pca$x

summary(pca) # 2 50%+, 3 60%+, 5 70%+, 9 80%+, 16 90%+


#2. 3rd order PC and U scores

order3 = Rob_TSM(Yp, L = 30, N = 30, order = 3)
phi3 = order3$eigenv
ev3 = order3$eigenl
t_phi3 = ginv(t(w)) %*% phi3
pdf("transforme.pdf",width = 10,height = 10)
par(mfrow = c(4,1))
for(i in 1:32){
  plot(t_phi3_1[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)

}
dev.off()
