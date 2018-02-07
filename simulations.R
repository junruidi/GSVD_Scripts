#some simulation regarding the methdos

# 1. visual inspection of different direction of 2nd order and 3rd order basis

rm(list = ls())
library(sn)
library(ggplot2)
library(data.table)
library(MASS)

# 1.1 skew normal distribution
xi = c(3,5)
Omega = matrix(c(2,2,2,4),2,2)
alpha = c(120, 0)

rnd = rmsn(n = 1000, xi = xi, Omega = Omega, alpha = alpha)
DT = scale(rnd, center = T, scale = F)
DT = data.table(DT)
names(DT) = c("x","y")

##First PCA
PCA1=prcomp(DT,center = F)

##Dataframe containing the arrow with principal directions
data_arrow=data.table(x= c(0,0),y= c(0,0),xend=c(PCA1$rotation[1,]),yend=c(PCA1$rotation[2,]))
moment3 = MGT3(scale(rnd, center = T, scale = F))
V3 = hoevd(moment3,rank = 2)

#####################################################
###                Simulations Setup              ###
###                J Di. 02/07/2018               ###
#####################################################


# 1. simulation on the direction of 2nd and third order -------------------
rm(list = ls())
library(mgcv)
library(dplyr)
library(data.table)
library(ggplot2)

source("applyall.R")
creat_cov = function(p){
  A = matrix(runif(p^2)*2-1, ncol=p) 
  return(t(A) %*% A)
}

# 1.1 multivariate normal distribution
V_mvn = creat_cov(2)
V_mvn
mu_mvn = c(0,0)

MVN = rmvn(n = 1000, mu = mu_mvn, V = V_mvn)
MVN = as.data.frame(scale(MVN, center = TRUE, scale = FALSE))
names(MVN) = c("x","y")

PCA1=prcomp(MVN,center = F)
data_arrow=data.table(x= c(0,0),y= c(0,0),xend=c(PCA1$rotation[1,]),yend=c(PCA1$rotation[2,]))
moment3 = MGT3(scale(MVN, center = T, scale = F))
V3 = hoevd(moment3,rank = 2)$u
data_arrow2=data.table(x= c(0,0),y= c(0,0),xend=c(V3[1,]),yend=c(V3[2,]))


data_arrow
data_arrow2

#this line subject to change depending on actual arrow direction
data_arrow2[1,] = - data_arrow2[1,]

ggplot(MVN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order v.s. 3rd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(MVN$x),y=y,yend=yend*max(MVN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue") + 
  geom_segment(data=data_arrow2,aes(x = x,xend=xend*max(MVN$x),y=y,yend=yend*max(MVN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


score_2 = data.frame(PCA1$x)
moment3 = MGT3(PCA1$x)
V3_2 = hoevd(moment3,rank = 2)$u


ggplot(MVN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(MVN$x),y=y,yend=yend*max(MVN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue")
data_arrow2_2=data.table(x= c(0,0),y= c(0,0),xend=c(V3_2[1,]),yend=c(V3_2[2,]))
ggplot(score_2,aes(PC1,PC2))+geom_point(alpha=0.2) + coord_fixed()+ggtitle('3rd order')+
  geom_segment(data=data_arrow2_2,aes(x = x,xend=xend*max(score_2$PC1),y=y,yend=yend*max(score_2$PC1)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")

#2. two univariate normal data
x = rnorm(1000, sd = 2)
y = rnorm(1000, sd = 4)
UNM = as.data.frame(scale(cbind(x,y), center = T,scale = F))

PCA1=prcomp(UNM,center = F)
data_arrow=data.table(x= c(0,0),y= c(0,0),xend=c(PCA1$rotation[1,]),yend=c(PCA1$rotation[2,]))
moment3 = MGT3(scale(UNM, center = T, scale = F))
V3 = hoevd(moment3,rank = 2)$u
data_arrow2=data.table(x= c(0,0),y= c(0,0),xend=c(V3[1,]),yend=c(V3[2,]))


data_arrow
data_arrow2

#this line subject to change depending on actual arrow direction
data_arrow2[1,] = - data_arrow2[1,]

ggplot(UNM,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order v.s. 3rd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(UNM$x),y=y,yend=yend*max(UNM$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue") + 
  geom_segment(data=data_arrow2,aes(x = x,xend=xend*max(UNM$x),y=y,yend=yend*max(UNM$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


score_2 = data.frame(PCA1$x)
moment3 = MGT3(PCA1$x)
V3_2 = hoevd(moment3,rank = 2)$u


ggplot(UNM,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(UNM$x),y=y,yend=yend*max(UNM$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue")
data_arrow2_2=data.table(x= c(0,0),y= c(0,0),xend=c(V3_2[1,]),yend=c(V3_2[2,]))
ggplot(score_2,aes(PC1,PC2))+geom_point(alpha=0.2) + coord_fixed()+ggtitle('3rd order')+
  geom_segment(data=data_arrow2_2,aes(x = x,xend=xend*max(score_2$PC1),y=y,yend=yend*max(score_2$PC1)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


