################################################################
##  Univariate Properties/Covariance, cumulant tensor slices  ##
##                    J Di 02/12/2018                         ##
################################################################

#1. univariate properties
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
library(qdap)
library(lubridate)
library(timeDate)
library(e1071)

TIME = char2end(as.character(timeSequence(from = hm("7:00"), to = hm("22:59"), by = "hour")),char = " ",noc=1)
TIME = beg2char(TIME,":",2)

Y = hr50[,-1]

varY = apply(Y,2,var)
skwY = apply(Y,2,skewness)
kurY = apply(Y,2,kurtosis)

pdf("Write Up/plots/univariate.pdf",width = 10,height = 10)
par(mfrow = c(3,1))

plot(varY,main = "Variance",type = "l",xaxt = "n",ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(skwY,main = "Skewness",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)
plot(kurY,main = "Kurtosis",type = "l",xaxt = "n", ylab = "")
axis(1, at = c(seq(1,32,2)),labels = TIME)

dev.off()

#2. Cumulant matrix and tensors...
rm(list = ls())
setwd("~/Dropbox/Junrui Di/tensor analysis/GSVD/")
load("Data/hr50.rda")
source("GSVDScripts/applyall.R")
library(corrplot)

Y = hr50[,-1]

#2.1 covariance
covY = cov(Y)
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
png(filename = "Write Up/cov2.png", width = 6, height = 5, units = "in",res = 300 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
corrplot(covY,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
         tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
dev.off()

#2.2. third order cumulant
Y = scale(hr50[,-1],center = T, scale = F)
cum3Y = third_cumulant_p(Y)

pdf(file = "Write Up/plots/c3_1.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
m1 = cum3Y[i,,]
corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
         tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()

pdf(file = "Write Up/plots/c3_2.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  m1 = cum3Y[,i,]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()

pdf(file = "Write Up/plots/c3_3.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  m1 = cum3Y[,,i]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()


#2.3. forth order cumulant
cum4Y = four_cumulants_direct(Y)
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
pdf(file = "Write Up/plots/c4_1_2.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  for(j in 1:32){
  m1 = cum4Y[i,j,,]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}}
dev.off()

pdf(file = "Write Up/plots/c4_2.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  m1 = cum3Y[,i,,]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()

pdf(file = "Write Up/plots/c4_3.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  m1 = cum3Y[,,i,]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()

pdf(file = "Write Up/plots/c4_4.pdf", width = 8,height = 5 )
par(mar = c(4,5,9,6))
par(oma = c(1,0,1,0))
for(i in 1:32){
  m1 = cum3Y[,,,i]
  corrplot(m1,method = "color",is.corr = F,cl.pos = "r",col = col4(200), cl.length = 21,cl.cex = 0.6,
           tl.cex = 0.8, tl.pos = "n",cl.align.text = "l",cl.ratio = 0.1, number.digits = 1)
}
dev.off()
