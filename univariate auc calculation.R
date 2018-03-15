rm(list = ls())
# setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
# load("Data/survival_gsvdc_centerscale.rda")
# source("GSVDScripts/bootauc_mort.R")
load("survival_gsvdc_centerscale.rda")
source("bootauc_mort.R")
library(survey)
library(pROC)

base_vars = c("age")

aucs.dat = data.frame()
for(i in 33:224){
  new_vars = names(survival_gsvd_center_scale)[i]
  auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 100 ,seed.start = 1243,
                         base_vars = base_vars, new_vars = new_vars)
  aucs.dat = rbind(aucs.dat,auc.all)
}

aucs.dat = as.data.frame(t(aucs.dat))
names(aucs.dat) = names(survival_gsvd_center_scale)[33:224]
row.names(aucs.dat) = NULL

save(aucs.dat, file = "aucs_age_univariate.rda")


