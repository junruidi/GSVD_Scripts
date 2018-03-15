rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/cov50.rda")

library(nhanesaccel)
cov50$yr5_mort = as.integer(ifelse(cov50$mort_yr <= 5 & cov50$mortstat == 1, 1,0))
cov50$permth_exm = NULL

include.column = which(names(cov50)=="include")
c1 = subset(cov50,yr==34)
d1 = subset(cov50,yr==56)

c2 = nhanes.accel.reweight(acceldata = c1, wave = 1, seqn.column = 1, 
                           include.column = include.column)
d2 = nhanes.accel.reweight(acceldata = d1, wave = 2, seqn.column = 1, 
                           include.column = include.column)
cov50$wt4yr = c(c2$wtmec2yr_adj/2,d2$wtmec2yr_adj/2)
normalize=function(v){return(v/sum(v))}
cov50$wt4yr_norm= normalize(cov50$wt4yr)*nrow(cov50)

save(cov50, file = "Data/cov50.rda")
