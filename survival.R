# 1. Survival TCA ---------------------------------------------------------
rm(list = ls())
library(survey)
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/survival_TCA.dta")
dat = svydata_TCA$variables
varname = NULL
result_tca = data.frame()
for(i in 32:74){
  dat.i = dat[,c(1:31,i)]
  var.i = names(dat.i)[32]
  names(dat.i)[32] = "x"
  svy.i = svydesign(id=~SDMVPSU,
                              strat=~SDMVSTRA,
                              weight=~wt4yr_norm,
                              nest=TRUE,
                              data=dat.i)
  
  cox= svycoxph(Surv(survyr,mortstat) ~ age + Male  + BMXBMI +
                  diabetes + CHF + stroke + cancer + CHD+
                  MobilityProblem + Mexican + HispanicOther + Black + OtherRace + 
                  LessThanHS + HighSchool + MissingEduc +CurrentSmoker + FormerSmoker + 
                  FormerDrinker + ModerateDrinker+ HeavyDrinker+ MissingAlcohol+
                  x,design=svy.i)
  tbl = summary(cox)$coefficients[23,]
  
  varname = c(varname, var.i)
  result_tca = rbind(result_tca, tbl)
}
result_tca = cbind(varname, result_tca)
names(result_tca) = c("var",names(tbl)) 

# 2. Survival gsvd center ---------------------------------------------------------
rm(list = ls())
library(survey)
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/survival_gsvdc.dta")
dat = svydata_gsvdc$variables
varname = NULL
result_gsvdc = data.frame()
for(i in 32:127){
  dat.i = dat[,c(1:31,i)]
  var.i = names(dat.i)[32]
  names(dat.i)[32] = "x"
  svy.i = svydesign(id=~SDMVPSU,
                    strat=~SDMVSTRA,
                    weight=~wt4yr_norm,
                    nest=TRUE,
                    data=dat.i)
  
  cox= svycoxph(Surv(survyr,mortstat) ~ age + Male  + BMXBMI +
                  diabetes + CHF + stroke + cancer + CHD+
                  MobilityProblem + Mexican + HispanicOther + Black + OtherRace + 
                  LessThanHS + HighSchool + MissingEduc +CurrentSmoker + FormerSmoker + 
                  FormerDrinker + ModerateDrinker+ HeavyDrinker+ MissingAlcohol+
                  x,design=svy.i)
  tbl = summary(cox)$coefficients[23,]
  
  varname = c(varname, var.i)
  result_gsvdc = rbind(result_gsvdc, tbl)
}
result_gsvdc = cbind(varname, result_gsvdc)
names(result_gsvdc) = c("var",names(tbl)) 