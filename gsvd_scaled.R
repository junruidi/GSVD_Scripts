# rm(list = ls())
# setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
# load("Data/hr50.rda")
# source("GSVDScripts/applyall.R")
# library(rTensor)
# Y = scale(hr50[,-1],center = T, scale = F)
# u = svd(Y)$u
# v = svd(Y)$v
# s = diag(svd(Y)$d)
# w = v %*% solve(s)
# 
# Yp = Y %*% w
# # a = NULL
# # for(i in 1:3390){
# # yp1 = Yp[i,]
# # check1 = t(w) %*% Y[i,]
# # a = c(a,cor(yp1, check1))
# # }
# 
# 
# mmy = MGT3(Y)
# tnsr_y = as.tensor(mmy)
# unfold_y = k_unfold(tnsr_y,m = 1)@data
# 
# uy = svd(unfold_y)$u
# sy = svd(unfold_y)$d
# vy = svd(unfold_y)$v
# 
# mmyp = MGT3(Yp)
# tnsr_p = as.tensor(mmyp)
# unfold_yp = k_unfold(tnsr_p,m = 1)@data
# 
# uyp = svd(unfold_yp)$u
# syp = svd(unfold_yp)$d
# vyp = svd(unfold_yp)$v
# 
# 
# check = v %*% uyp
# 
# 
# more = t(w) %*% uy %*% diag(sy) %*% t(vy) %*% (w %x% w)

rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/applyall.R")
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
moment3 = MGT3(Yp)
hosvd3_v = hoevd(moment3,rank = 32)
phi3 = hosvd3_v$u

core3_v = hosvd3_v$z
unfold_3_1_v = k_unfold(as.tensor(core3_v),2)@data
svd_3_v = svd(unfold_3_1_v)$d
pct_ev3_v = cumsum(svd_3_v)/sum(svd_3_v)  #11 50%+, 14 60%+, 17 70%+, 22 80%+, 26 90%+

sc3 = Yp %*% phi3 

hosvd3_u = Gram3_hosvd(Yp)

#3. 4th order PC and U scores
moment4 = MGT4(Yp)
hosvd4_v = hoevd(moment4,rank = 32)
phi4 = hosvd4_v$u

core4_v = hosvd4_v$z
unfold_4_1_v = k_unfold(as.tensor(core4_v),1)@data
svd_4_v = svd(unfold_4_1_v)$d
pct_ev4_v = cumsum(svd_4_v)/sum(svd_4_v) # 8 50%+, 11 60%+, 15 70%+, 19 80%+, 25 90%+


sc4 = Yp %*% phi4 
hosvd4_u = Gram4_hosvd(Yp)


phi3 = ginv(t(w)) %*% phi3
phi4 = ginv(t(w)) %*% phi4


phi = as.data.frame(cbind(phi2, phi3, phi4))
names(phi) = paste0("phi",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

sc = as.data.frame(cbind(sc2,sc3,sc4))
names(sc) = paste0("SC",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))

u_score = as.data.frame(cbind(u,hosvd3_u$u,hosvd4_u$u))
names(u_score) = paste0("US",rep(c(2,3,4),each = 32),"_",rep(c(1:32),3))


score_all = cbind(sc, u_score)

pdf(file = "Write Up/cor_gsvd_scaled_center.pdf", width = 56, height = 56)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(score_all),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/cor_gsvd_scaled_center_pcscore.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(sc),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()

pdf(file = "Write Up/cor_gsvd_scaled_center_mdsscore.pdf", width = 28, height = 28)
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot::corrplot(cor(u_score),cl.pos = "b",cl.cex = 2,tl.cex = 1.6,cl.align.text = "r",na.label = "-")
dev.off()



library(qdap)
library(timeDate)
library(lubridate)
TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)
# for(i in 1:ncol(phi2)){
#   
#   sign2 = sign(phi2[1,i])
#   sign23 = sign(phi3[1,i])
#   sign234 = sign(phi4[1,i])
#   if(sign2 != sign23){
#     phi3[,i] = -phi3[,i]
#   }
#   if(sign2 != sign234){
#     phi4[,i] = -phi4[,i]
#   }
# }



pdf("Write Up/gsvd_phi_scaled_center_transfer.pdf",width = 10,height = 10)
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

survival_gsvd_center_scale = cbind(cov50,score_all)


save(survival_gsvd_center_scale, file = "Data/survival_gsvdc_centerscale.rda")



#########################################
#prediction model
rm(list = ls())
setwd("/Users/junruidi/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/survival_gsvdc_centerscale.rda")
library(survey)
library(pROC)
source("GSVDScripts/bootauc_mort.R")
base_vars = c("age", "Male", "BMXBMI", 
              "diabetes", "CHF", "stroke", "cancer", "CHD", 
              "MobilityProblem", "Mexican", "HispanicOther", "Black", "OtherRace", 
              "LessThanHS", "HighSchool", "MissingEduc","CurrentSmoker", "FormerSmoker", 
              "FormerDrinker", "ModerateDrinker", "HeavyDrinker", "MissingAlcohol")
base_vars = c("age")
svydata = svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt4yr_norm, data = survival_gsvd_center_scale, nest = TRUE)

orders = c(2, 3, 4)
p50 = c(2, 11, 8)
p60 = c(3, 14, 11)
p70 = c(5, 17, 15)
p80 = c(9, 22, 19)
p90 = c(16, 26, 25)
npc = cbind.data.frame(orders, p50, p60, p70, p80, p90)

#1. model 1: all pcs

new_vars = paste0("SC2_",c(1:npc$p60[1]))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_1_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())
  
summary(model_1_allpcs, df.resid=Inf) #pick up 1

new_vars = paste0("SC2_",1)
form = paste(c(base_vars, new_vars), collapse=" + ")
model1 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())
summary(model1, df.resid=Inf)

auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all) #0.8041292

#2. model 2 add in 3rd oder
new_vars = c(paste0("SC2_",c(1:npc$p60[1])), paste0("SC3_",c(1:npc$p60[2])))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_2_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())

summary(model_2_allpcs, df.resid=Inf) #pick up 1
new_vars = c(paste0("SC2_",1), paste0("SC3_",c(8,13,14)))
form = paste(c(base_vars, new_vars), collapse=" + ")
model2 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                design=svydata, family=quasibinomial())
summary(model2, df.resid=Inf)
auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all)  #0.8113471


#3. model 3 add in 4th oder
new_vars = c(paste0("SC2_",c(1:npc$p60[1])), paste0("SC4_",c(1:npc$p60[3])))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_3_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())

summary(model_3_allpcs, df.resid=Inf) #pick up 1
new_vars = c(paste0("SC2_",1), paste0("SC3_",c(4,9)))
form = paste(c(base_vars, new_vars), collapse=" + ")
model3 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                design=svydata, family=quasibinomial())
summary(model3, df.resid=Inf)
auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all) #0.8099764


##########
#prediction model
rm(list = ls())
setwd("/Users/junruidi/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/survival_gsvdc_centerscale.rda")
library(survey)
library(pROC)
source("GSVDScripts/bootauc_mort.R")
base_vars = c("age")
svydata = svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt4yr_norm, data = survival_gsvd_center_scale, nest = TRUE)

#1. model 1: all pcs

new_vars = paste0("SC2_",c(1:32))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_1_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())

summary(model_1_allpcs, df.resid=Inf) #pick up 1

new_vars = paste0("SC2_",1)
form = paste(c(base_vars, new_vars), collapse=" + ")
model1 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                design=svydata, family=quasibinomial())
summary(model1, df.resid=Inf)

auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all) #0.8041292
#0.7868283


#2. model 2 add in 3rd oder
new_vars= paste0("SC3_",c(1:32))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_2_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())

summary(model_2_allpcs, df.resid=Inf) #pick up 1
new_vars = c(paste0("SC2_",1), paste0("SC3_",c(8,13,14)))
form = paste(c(base_vars, new_vars), collapse=" + ")
model2 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                design=svydata, family=quasibinomial())
summary(model2, df.resid=Inf)
auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all)  #0.8113471
#0.7085179

#3. model 3 add in 4th oder
new_vars = paste0("SC4_",c(1:32))
form = paste(c(base_vars, new_vars), collapse=" + ")
model_3_allpcs = svyglm(as.formula(paste("yr5_mort ~", form)), 
                        design=svydata, family=quasibinomial())

summary(model_3_allpcs, df.resid=Inf) #pick up 1
new_vars = c(paste0("SC2_",1), paste0("SC3_",c(4,9)))
form = paste(c(base_vars, new_vars), collapse=" + ")
model3 = svyglm(as.formula(paste("yr5_mort ~", form)), 
                design=svydata, family=quasibinomial())
summary(model3, df.resid=Inf)
auc.all = bootauc_mort(dat = survival_gsvd_center_scale,nboot = 1000,seed.start = 1243,
                       base_vars = base_vars, new_vars = new_vars)
mean(auc.all) #0.596319


