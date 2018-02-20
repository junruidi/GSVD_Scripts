rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/applyall.R")

#1. 2nd order
pc.hr = prcomp(x = hr50[,-1], scale. = F, center = T)
PC2 = pc.hr$rotation[,1:16] # up to 90%
SC2 = pc.hr$x[,1:16]

#2. 3rd order
Y = scale(hr50[,-1],center = T,scale = F)
CUM3 = third_cumulant(Y)
hosvd3 = hoevd(CUM3, 32)

PC3 = hosvd3$u
core3 = hosvd3$z
unfold_3_1 = k_unfold(as.tensor(core3),1)@data
svd_3 = svd(unfold_3_1)$d
pct_ev3 = cumsum(svd_3)/sum(svd_3) #15

PC3 = PC3[,1:15]
SC3 = Y %*% PC3

#or looking at TCA
CUM3_2sc = third_cumulant(SC2)
hosvd3_2sc = hoevd(CUM3_2sc,16)

PC3_2sc = hosvd3_2sc$u
core3_2sc = hosvd3_2sc$z
unfold_3_1_2sc = k_unfold(as.tensor(core3_2sc),1)@data
svd_3_2sc = svd(unfold_3_1_2sc)$d
pct_ev3_2sc = cumsum(svd_3_2sc)/sum(svd_3_2sc) #10
PC3_2sc = PC3_2sc[,1:10]

PC3_p = PC2 %*% PC3_2sc
SC3_p = Y %*% PC3_p


#3. 4th order
CUM4 = four_cumulants(Y)
hosvd4 = hoevd(CUM4, 32)

PC4 = hosvd4$u
core4 = hosvd4$z
unfold_4_1 = k_unfold(as.tensor(core4),1)@data
svd_4 = svd(unfold_4_1)$d
pct_ev4 = cumsum(svd_4)/sum(svd_4) #12

PC4 = PC4[,1:12]
SC4 = Y %*% PC4

#or looking at TCA
CUM4_3sc = four_cumulants(SC3_p)
hosvd4_3sc = hoevd(CUM4_3sc,10)

PC4_3sc = hosvd4_3sc$u
core4_3sc = hosvd4_3sc$z
unfold_4_1_3sc = k_unfold(as.tensor(core4_3sc),1)@data
svd_4_3sc = svd(unfold_4_1_3sc)$d
pct_ev4_3sc = cumsum(svd_4_3sc)/sum(svd_4_3sc) #6
PC4_3sc = PC4_3sc[,1:6]

PC4_p = PC3_p %*% PC4_3sc
SC4_p = Y %*% PC4_p


