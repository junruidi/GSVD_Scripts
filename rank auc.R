rm(list = ls())
setwd("D:/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/aucs_all.rda")
library(qdap)
library(qdapRegex)

age = data.frame(var = names(auc.age), auc = colMeans(auc.age))
age = age[order(age$auc,decreasing = T),]
noage = data.frame(var = names(auc.noage), auc = colMeans(auc.noage))
noage = noage[order(noage$auc,decreasing = T),] 

age2 = age[-which(grepl("u[0-9]",age$var)),]
age2 = age2[-which(grepl("convert",age2$var)),]


age2$method = beg2char(age2$var,"_")
age2$order = rm_between(age2$var,"_","_",extract = T)
age2$pc = as.numeric(char2end(age2$var,"_",noc = 2))

noage2 = noage[-which(grepl("u[0-9]",noage$var)),]
noage2 = noage2[-which(grepl("convert",noage2$var)),]

noage2$method = beg2char(noage2$var,"_")
noage2$order = rm_between(noage2$var,"_","_",extract = T)
noage2$pc = as.numeric(char2end(noage2$var,"_",noc = 2))


age2_12 = subset(age2,pc <= 12)
row.names(age2_12) = NULL


noage2_12 = subset(noage2,pc <= 12)
row.names(noage2_12) = NULL
