rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/applyall.R")

Y = scale(hr50[,-1],center = T, scale = F)
u = svd(Y)$u
v = svd(Y)$v
s = diag(svd(Y)$d)

Yp = Y %*% v %*% solve(s)

#1. 2nd order
pca = prcomp(Y, center = F, scale. = F)
phi2 = pca$rotation
sc2 = pca$x


#2. 3rd order
moment3 = MGT3(Yp)
hosvd3_v = hoevd(moment3,rank = 32)
phi3 = hosvd3_v$u

core3_v = hosvd3_v$z
unfold_3_1_v = k_unfold(as.tensor(core3_v),2)@data
svd_3_v = svd(unfold_3_1_v)$d
pct_ev3_v = cumsum(svd_3_v)/sum(svd_3_v) 

sc3 = Yp %*% phi3

# 3. 4th
moment4 = MGT4(Yp)
hosvd4_v = hoevd(moment4,rank = 32)
phi4 = hosvd4_v$u


core4_v = hosvd4_v$z
unfold_4_1_v = k_unfold(as.tensor(core4_v),1)@data
svd_4_v = svd(unfold_4_1_v)$d
pct_ev4_v = cumsum(svd_4_v)/sum(svd_4_v) # 12


sc4 = Yp %*% phi4

phi = as.data.frame(cbind(phi2, phi3, phi4))
names(phi) = paste0("phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

sc = as.data.frame(cbind(sc2,sc3,sc4))
names(sc) = paste0("SC",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


# plots 
pdf(file = "Write Up/cor_gsvd_scaled_center.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(sc),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()


library(qdap)
library(timeDate)
library(lubridate)
TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)
for(i in 1:ncol(phi2)){
  
  sign2 = sign(phi2[1,i])
  sign23 = sign(phi3[1,i])
  sign234 = sign(phi4[1,i])
  if(sign2 != sign23){
    phi3[,i] = -phi3[,i]
  }
  if(sign2 != sign234){
    phi4[,i] = -phi4[,i]
  }
}



pdf("Write Up/gsvd_phi_scaled_center.pdf",width = 10,height = 10)
par(mfrow = c(3,1))
for(i in 1:32){
  plot(phi2[,i],main = paste0("phi2 - ",i),type = "l",xaxt = "n",ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi3[,i],main = paste0("phi3 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
  plot(phi4[,i],main = paste0("phi4 - ",i),type = "l",xaxt = "n", ylab = "RSV")
  axis(1, at = c(seq(1,32,2)),labels = TIME)
  abline(h = 0,lty = 3)
}
dev.off()

load("Data/cov50.rda")

survival_gsvd_center = cbind(cov50,SC)


save(svydata_gsvdc, file = "Data/survival_gsvdc_centerscale.rda")


