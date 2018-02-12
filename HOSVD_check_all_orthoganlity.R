rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
source("GSVDScripts/applyall.R")
load("Data/hr50.rda")
Y = scale(hr50[,-1], center = T, scale = F)
Y[, 32] = Y[,31]*5
moment3 = MGT3(Y)
rm(list = setdiff(ls(),c("moment3","MGT3")))

hosvd3 = hosvd(as.tensor(moment3),c(32,32,32))
S = hosvd3$Z@data
for(i in 2:32){
  print(sum(S[1,,]*S[i,,]))
}

S1 = k_unfold(as.tensor(S),m = 1)@data


#####
x = matrix(rnorm(32000,mean = 1000, sd = 1),1000,32)
x = scale(x, center = T, scale = F)
moment3 = MGT3(x)
hosvd3 = hosvd(as.tensor(moment3),c(32,32,32))
S = hosvd3$Z@data
for(i in 2:32){
  print(sum(S[1,,]*S[i,,]))
}

x = matrix(rnorm(32000,mean = 1000, sd = 100),1000,32)
x = scale(x, center = T, scale = F)
moment3 = MGT3(x)
hosvd3 = hosvd(as.tensor(moment3),c(32,32,32))
S = hosvd3$Z@data
for(i in 2:32){
  print(sum(S[1,,]*S[i,,]))
}

####
tnsr <- rand_tensor(c(6,7,8))
hosvdD <-hosvd(tnsr)
z = hosvdD$Z@data
sum(z[1,,] * z[2,,])
sum(z[1,,] * z[3,,])
sum(z[1,,] * z[4,,])
