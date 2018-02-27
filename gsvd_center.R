rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/applyall.R")

Y = scale(hr50[,-1], center = T, scale = F)

#1. 2nd order
SVD2 = svd(Y)
U2 = SVD2$u
V2 = SVD2$v

d2 = SVD2$d
pct_ev2 = cumsum(d2)/sum(d2) # 25


#2. 3rd 
moment3 = MGT3(Y)
hosvd3_v = hoevd(moment3,rank = 32)
hosvd3_u = Gram3_hosvd(Y)

V3 = hosvd3_v$u
U3 = hosvd3_u$u
S3 = t(U3) %*% Y %*% V3

core3_v = hosvd3_v$z
unfold_3_1_v = k_unfold(as.tensor(core3_v),1)@data
svd_3_v = svd(unfold_3_1_v)$d
pct_ev3_v = cumsum(svd_3_v)/sum(svd_3_v) # 15

core3_u = hosvd3_u$middle
svd_3_u = svd(core3_u)$d
pct_ev3_u = cumsum(svd_3_u)/sum(svd_3_u) # 17


# 3. 4th
moment4 = MGT4(Y)
hosvd4_v = hoevd(moment4,rank = 32)
hosvd4_u = Gram4_hosvd(Y)

V4 = hosvd4_v$u
U4 = hosvd4_u$u
S4 = t(U4) %*% Y %*% V4

core4_v = hosvd4_v$z
unfold_4_1_v = k_unfold(as.tensor(core4_v),1)@data
svd_4_v = svd(unfold_4_1_v)$d
pct_ev4_v = cumsum(svd_4_v)/sum(svd_4_v) # 12

core4_u = hosvd4_u$middle
svd_4_u = svd(core4_u)$d
pct_ev4_u = cumsum(svd_4_u)/sum(svd_4_u) # 17


V = as.data.frame(cbind(V2,V3,V4))
names(V) = paste0("V",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

U = as.data.frame(cbind(U2,U3,U4))
names(U) = paste0("U",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# plots 
pdf(file = "Write Up/cor_gsvd_center.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(U),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/S_center.pdf", width = 15, height = 15)
par(mfrow = c(3,3))
heatmap(diag(SVD2$d),Rowv = NA, Colv = NA)
heatmap(S3,Rowv = NA, Colv = NA)
heatmap(S4,Rowv = NA, Colv = NA)
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



pdf("Write Up/gsvd_V_center.pdf",width = 10,height = 10)
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

load("Data/cov50.rda")

survival_gsvd_center = cbind(cov50,U)

library(survey)
svydata_gsvdc = svydesign(id=~SDMVPSU,
                        strat=~SDMVSTRA,
                        weight=~wt4yr_norm,
                        nest=TRUE,
                        data=survival_gsvd_center)

save(svydata_gsvdc, file = "Data/survival_gsvdc.dta")
