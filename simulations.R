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
library(sn)

source("applyall.R")
creat_cov = function(p, min = 0, max = 1){
  A = matrix(runif(p^2,min = min , max = max)*2-1, ncol=p) 
  return(t(A) %*% A)
}

vv = creat_cov(2)
vv

# 1.1 multivariate normal distribution
V_mvn = vv
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
data_arrow2[2,] = - data_arrow2[2,]



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
y = rnorm(1000, sd = 2)
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


#3. multivariate skewed normal
V_mvsn = creat_cov(p = 2, min = -5, max = 5)
V_mvsn

MVSN = rmsn(n = 1000, xi = c(0,0), Omega = V_mvsn, alpha = c(-10,35))
MVSN = as.data.frame(scale(MVSN, center = TRUE, scale = FALSE))
names(MVSN) = c("x","y")

plot(MVSN$x, MVSN$y)


PCA1=prcomp(MVSN,center = F)
data_arrow=data.table(x= c(0,0),y= c(0,0),xend=c(PCA1$rotation[1,]),yend=c(PCA1$rotation[2,]))
moment3 = MGT3(scale(MVSN, center = T, scale = F))
V3 = hoevd(moment3,rank = 2)$u
data_arrow2=data.table(x= c(0,0),y= c(0,0),xend=c(V3[1,]),yend=c(V3[2,]))


data_arrow
data_arrow2

#this line subject to change depending on actual arrow direction
data_arrow2[1,] = - data_arrow2[1,]
data_arrow2[2,] = - data_arrow2[2,]

ggplot(MVSN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order v.s. 3rd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(MVSN$x),y=y,yend=yend*max(MVSN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue") + 
  geom_segment(data=data_arrow2,aes(x = x,xend=xend*max(MVSN$x),y=y,yend=yend*max(MVSN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


score_2 = data.frame(PCA1$x)
moment3 = MGT3(PCA1$x)
V3_2 = hoevd(moment3,rank = 2)$u


ggplot(MVSN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(MVSN$x),y=y,yend=yend*max(MVSN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue")
data_arrow2_2=data.table(x= c(0,0),y= c(0,0),xend=c(V3_2[1,]),yend=c(V3_2[2,]))
ggplot(score_2,aes(PC1,PC2))+geom_point(alpha=0.2) + coord_fixed()+ggtitle('3rd order')+
  geom_segment(data=data_arrow2_2,aes(x = x,xend=xend*max(score_2$PC1),y=y,yend=yend*max(score_2$PC1)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


#4. two univariate sn
V_mvsn = creat_cov(p = 2, min = -25, max = 25)
V_mvsn
x = rsn(1000, xi = 0, omega = V_mvsn[1,1], alpha = 0)
y = rsn(1000, xi = 0, omega = V_mvsn[2,2], alpha = 75)
USN = as.data.frame(scale(cbind(x,y), center = T,scale = F))
names(USN) = c("x","y")

plot(USN$x, USN$y)


PCA1=prcomp(USN,center = F)
data_arrow=data.table(x= c(0,0),y= c(0,0),xend=c(PCA1$rotation[1,]),yend=c(PCA1$rotation[2,]))
moment3 = MGT3(scale(USN, center = T, scale = F))
V3 = hoevd(moment3,rank = 2)$u
data_arrow2=data.table(x= c(0,0),y= c(0,0),xend=c(V3[1,]),yend=c(V3[2,]))


data_arrow
data_arrow2

#this line subject to change depending on actual arrow direction
data_arrow2[1,] = - data_arrow2[1,]
data_arrow2[2,] = - data_arrow2[2,]

ggplot(USN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order v.s. 3rd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(USN$x),y=y,yend=yend*max(USN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue") + 
  geom_segment(data=data_arrow2,aes(x = x,xend=xend*max(USN$x),y=y,yend=yend*max(USN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")


score_2 = data.frame(PCA1$x)
moment3 = MGT3(PCA1$x)
V3_2 = hoevd(moment3,rank = 2)$u


ggplot(USN,aes(x,y))+geom_point(alpha=0.2)+ coord_fixed()+ggtitle('2nd order')+
  geom_segment(data=data_arrow,aes(x = x,xend=xend*max(USN$x),y=y,yend=yend*max(USN$x)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="blue")
data_arrow2_2=data.table(x= c(0,0),y= c(0,0),xend=c(V3_2[1,]),yend=c(V3_2[2,]))
ggplot(score_2,aes(PC1,PC2))+geom_point(alpha=0.2) + coord_fixed()+ggtitle('3rd order')+
  geom_segment(data=data_arrow2_2,aes(x = x,xend=xend*max(score_2$PC1),y=y,yend=yend*max(score_2$PC1)),
               arrow = arrow(length = unit(0.03, "npc")),size=2,color="red")
