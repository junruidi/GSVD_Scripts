#Generalized SVD
rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/GSVD/")
source("GSVDScripts/applyall.R")
load("Data/hr50.rda")
Y = scale(hr50[,-1], center = T, scale = F)

#1. 2nd order
SVD2 = svd(Y)
U2 = SVD2$u
V2 = SVD2$v
S2 = diag(SVD2$d)

#2. 3rd 
moment3 = MGT3(Y)
V3 = hoevd(moment3,rank = 32)$u
core_v3 = hoevd(moment3,rank = 32)$z
unfold_core_v3_1 = k_unfold(as.tensor(core_v3),m=1)@data
unfold_core_v3_2 = k_unfold(as.tensor(core_v3),m=2)@data
unfold_core_v3_3 = k_unfold(as.tensor(core_v3),m=3)@data

core_v3_1 = core_v3[1,,]
core_v3_2 = core_v3[2,,]
core_v3_3 = core_v3[3,,]


U3 = Gram3_hosvd(Y)



S3 = t(U3) %*% Y %*% V3

util3 = svd(S3)$u
vtil3 = svd(S3)$v

u3u3t = U3 %*% util3
v3v3t = V3 %*% vtil3

#3. 4th
moment4 = MGT4(Y)
V4= hoevd(moment4,rank = 32)
U4 = Gram4_hosvd(Y)
S4 = t(U4) %*% Y %*% V4

util4 = svd(S4)$u
vtil4 = svd(S4)$v


u4u4t = U4 %*% util4
v4v4t = V4 %*% vtil4







pdf(file = "result/0117/corplot_U.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(U),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()


pdf(file = "result/0117/center moment/corplot_US.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(US),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()



load("data/cov50.rda")
cov = subset(cov50, select = c(ID,Male,MobilityProblem,mortstat,cancer, diabetes))

col11 = rgb(0, 0, 238, alpha=100, maxColorValue=255)
col22 = rgb(205,15,20,maxColorValue=255)

u2 = U[,c(1:10)]
u3 = U[,c(33:42)]
u4 = U[,c(65:74)]

u = cbind(u2,u3,u4)

pdf(file = "result/0117/heamaps_S.pdf", width = 15, height = 15)
par(mfrow = c(3,3))
heatmap(diag(SVD2$d),Rowv = NA, Colv = NA)
heatmap(S3,Rowv = NA, Colv = NA)
heatmap(S4,Rowv = NA, Colv = NA)
dev.off()


pdf(file = "result/0117/Upairs_mortality.pdf", width = 15, height = 15)
par(mfrow = c(3,3))
for(i in 1:29){
  for(j in (i+1):30 ){
    plot(u[,i],u[,j],xlab = names(u)[i],ylab = names(u)[j],col=c(col11,col22)[as.factor(cov$mortstat)],pch = 19)
  }
}
dev.off()

pdf(file = "result/0117/Upairs_mobility.pdf", width = 15, height = 15)
par(mfrow = c(3,3))
for(i in 1:29){
  for(j in (i+1):30 ){
    plot(u[,i],u[,j],xlab = names(u)[i],ylab = names(u)[j],col=c(col11,col22)[as.factor(cov$MobilityProblem)],pch = 19)
  }
}
dev.off()


library(qdap)
library(timeDate)
library(lubridate)
TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)


for(i in 1:ncol(V2)){
  
  sign2 = sign(V2[1,i])
  sign23 = sign(V3[1,i])
  sign234 = sign(V4[1,i])
  if(sign2 != sign23){
    V3[,i] = -V3[,i]
  }
  if(sign2 != sign234){
    V4[,i] = -V4[,i]
  }
}



pdf("result/0117/V_comparison.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(V2[,i],main = paste0("V2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(V3[,i],main = paste0("V3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(V4[,i],main = paste0("V4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

pdf(file = "result/0117/corplot_V.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(V),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

