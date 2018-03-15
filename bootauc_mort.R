#boot auc for mortality models
bootauc_mort = function(dat, nboot, seed.start,base_vars, new_vars){
 
  N = nrow(dat)
  
  auc.all = NULL
  for(i in 1:nboot){
    set.seed(i + seed.start)
    
    ## use same bootstrapped data set for each iteration of the forward selection procedure
    inx      = sample(1:N, replace=TRUE, size=N, prob=dat$wt4yr_norm)
    inx_pred = sample(1:N, replace=TRUE, size=N, prob=dat$wt4yr_norm)
    
    dat_tmp = dat[inx,]
    dat_tmp_pred = dat[inx_pred,]
    
    form = paste(c(base_vars, new_vars), collapse=" + ")
    fit_tmp = glm(as.formula(paste("yr5_mort ~", form)), data=dat_tmp,family=quasibinomial())
    
    auc.i = roc(dat_tmp_pred$yr5_mort,predict(fit_tmp, newdata=dat_tmp_pred,type='response'))$auc[1]
    
    rm(list=c("fit_tmp","form"))
    auc.all = c(auc.all, auc.i)
  }
  return(auc.all)
}