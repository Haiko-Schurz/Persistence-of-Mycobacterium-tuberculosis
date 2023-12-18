#Install packages
install.packages("car")
install.packages("nlshelper")
install.packages("scales")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("ggplot2")
install.packages("dplyr") 

#Load libraries
library(car)
library(nlshelper)
library(scales)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr) 

#Function
standard_error<- function(x) sd(x) / sqrt(length(x))

########### Metabolic activity NEW ###########################################################################
meta.c<-read.table("sorted_counts_metabolic.txt",header=T)

#Set missing values to NA
meta.c$D0_Bac[c(29,30)]<-NA
meta.c$D3_MP[c(36:39)]<-NA
meta.c$D3_m_MP_High_red[c(36:39)]<-NA
meta.c$D3_m_MP_Low_red[c(36:39)]<-NA
meta.c$D5_MP[c(31:39)]<-NA
meta.c$D5_m_MP_High_red[c(31:39)]<-NA
meta.c$D5_m_MP_Low_red[c(31:39)]<-NA

#Correlation analysis
#Day 0 Bacteria vs. MFI (metabolic activity) D3 and D5 for active (low red) and persister (High red)

#Day 0 bacteria vs. D3 persister MFI
l <- loess(D3_m_MP_High_red ~ D0_Bac, data = meta.c)
with(meta.c, plot(D0_Bac, D3_m_MP_High_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D3 persister MFI", xlab="D0 Bacteria/ml"),ylab="D3 persister MFI", xlab="D0 Bacteria/ml"))
myline.fit <- lm( meta.c$D3_m_MP_High_red ~meta.c$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D3_m_MP_High_red , meta.c$D0_Bac)

#Day 0 bacteria vs. D3 active MFI
l <- loess(D3_m_MP_Low_red ~ D0_Bac, data = meta.c)
with(meta.c, plot(D0_Bac, D3_m_MP_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D3 active MFI", xlab="D0 Bacteria/ml"),ylab="D3 active MFI", xlab="D0 Bacteria/ml"))
myline.fit <- lm( meta.c$D3_m_MP_Low_red ~meta.c$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D3_m_MP_Low_red , meta.c$D0_Bac)

#Day 0 bacteria vs. D5 persister MFI
l <- loess(D5_m_MP_High_red ~ D0_Bac, data = meta.c)
with(meta.c, plot(D0_Bac, D5_m_MP_High_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D5 persister MFI", xlab="D0 Bacteria/ml"),ylab="D5 persister MFI", xlab="D0 Bacteria/ml"))
myline.fit <- lm( meta.c$D5_m_MP_High_red ~meta.c$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D5_m_MP_High_red , meta.c$D0_Bac)

#Day 0 bacteria vs. D5 active MFI
l <- loess(D5_m_MP_Low_red ~ D0_Bac, data = meta.c)
with(meta.c, plot(D0_Bac, D5_m_MP_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D5 active MFI", xlab="D0 Bacteria/ml"),ylab="D5 active MFI", xlab="D0 Bacteria/ml"))
myline.fit <- lm( meta.c$D5_m_MP_Low_red ~meta.c$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D5_m_MP_Low_red , meta.c$D0_Bac)

#Day 3 and 5 macrophage bacteria MIF vs. MIFcfor persistrs adn active #####################################################

#Day 3 macrophge bacteria vs. Day 3 prsisters
l <- loess(D3_m_MP_High_red ~ D3_MP, data = meta.c)
with(meta.c, plot(D3_MP, D3_m_MP_High_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D3 persister MFI", xlab="D3 MP Bacteria/ml"),ylab="D3 persister MFI", xlab="D3 MP Bacteria/ml"))
myline.fit <- lm( meta.c$D3_m_MP_High_red ~meta.c$D3_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D3_m_MP_High_red , meta.c$D3_MP)

#Day 3 macrophge bacteria vs. Day 3 active
l <- loess(D3_m_MP_Low_red ~ D3_MP, data = meta.c)
with(meta.c, plot(D3_MP, D3_m_MP_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D3 active MFI", xlab="D3 MP Bacteria/ml"),ylab="D3 active MFI", xlab="D3 MP Bacteria/ml"))
myline.fit <- lm( meta.c$D3_m_MP_Low_red ~meta.c$D3_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D3_m_MP_Low_red , meta.c$D3_MP)

#Day 5 macrophge bacteria vs. Day 5 prsisters
l <- loess(D5_m_MP_High_red ~ D5_MP, data = meta.c)
with(meta.c, plot(D5_MP, D5_m_MP_High_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D5 persister MFI", xlab="D5 MP Bacteria/ml"),ylab="D5 persister MFI", xlab="D5 MP Bacteria/ml"))
myline.fit <- lm( meta.c$D5_m_MP_High_red ~meta.c$D5_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D5_m_MP_High_red , meta.c$D5_MP)

#Day 5 macrophge bacteria vs. Day 5 active
l <- loess(D5_m_MP_Low_red ~ D5_MP, data = meta.c)
with(meta.c, plot(D5_MP, D5_m_MP_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE,ylab="D5 active MFI", xlab="D5 MP Bacteria/ml"),ylab="D5 active MFI", xlab="D5 MP Bacteria/ml"))
myline.fit <- lm( meta.c$D5_m_MP_Low_red ~meta.c$D5_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(meta.c$D5_m_MP_Low_red , meta.c$D5_MP)

#Day 3 and Day 5 active vs. persisters correlation analysis and plots

#Day 3 active vs. persisters
png(file="D3persister MFI vs D3active MFI.png",
    width=1500, height=1500, res=300)
l <- loess(D3_m_MP_High_red ~ D3_m_MP_Low_red, data = meta.c)
with(meta.c, plot(D3_m_MP_Low_red, D3_m_MP_High_red, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,8.5e3),xlim=c(0,7.2e3),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,8.5e3),xlim=c(0,7.2e3))))
axis(1, at = c( 1e3, 3e3, 5e3, 7e3), labels = c("1e+03", "3e+03", "5e+03", "7e+03"))
axis(2, at = c(0e00, 2e3, 4e3, 6e3,8e3), labels = c( "0e+00","2e+03", "4e+03", "6e+03", "8e+03"))
myline.fit <- lm( meta.c$D3_m_MP_High_red ~meta.c$D3_m_MP_Low_red)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(meta.c$D3_m_MP_High_red , meta.c$D3_m_MP_Low_red)

#Day 5 active vs. persisters
png(file="D5persister MFI vs D5active MFI.png",
    width=1500, height=1500, res=300)
l <- loess(D5_m_MP_High_red ~ D5_m_MP_Low_red, data = meta.c)
with(meta.c, plot(D5_m_MP_Low_red, D5_m_MP_High_red, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,8.5e3),xlim=c(0,7.2e3),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,8.5e3),xlim=c(0,7.2e3))))
axis(1, at = c( 1e3, 3e3, 5e3, 7e3), labels = c("1e+03", "3e+03", "5e+03", "7e+03"))
axis(2, at = c(0e00, 2e3, 4e3, 6e3,8e3), labels = c( "0e+00","2e+03", "4e+03", "6e+03", "8e+03"))
myline.fit <- lm( meta.c$D5_m_MP_High_red ~meta.c$D5_m_MP_Low_red)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(meta.c$D5_m_MP_High_red , meta.c$D5_m_MP_Low_red)

#Boxplot of MFi distribution for Day 3 and 5 active and persisters
png(file="D3 and D5 active and persistr MFI.png",
    width=1000, height=1000)

png(file="Box Plot D3 and D5 persister MFI (grey) and active MFI (orange).png",
    width=1500, height=1500, res=300)
boxplot(meta.c$D3_m_MP_High_red, meta.c$D3_m_MP_Low_red,
        meta.c$D5_m_MP_High_red, meta.c$D5_m_MP_Low_red,
        ylab="MFI",yaxt='n',
        names=c("Day3","Day3","Day5","Day5"),at = c(1,2,3,4) ,col=c("grey", "orange","grey", "orange"))
        axis(2, at = c(0e00, 2e3, 4e3, 6e3,8e3), labels = c( "0e+00","2e+03", "4e+03", "6e+03", "8e+03"))
#legend("topleft", c("High Red" ,"Low Red"), col=c("red","orange"), pch=c(15,15))
dev.off()


#Line plot with error bars for the active and persisters MFI at day 3 and Day 5
#Clculate required values
day<-c(3,5,3,5)
high.low<-c("Low red","Low red","High red","High red")
med<-c(median(na.omit(meta.c$D3_m_MP_Low_red)),median(na.omit(meta.c$D5_m_MP_Low_red)),median(na.omit(meta.c$D3_m_MP_High_red)),median(na.omit(meta.c$D5_m_MP_High_red)))
se<-c(standard_error(na.omit(meta.c$D3_m_MP_Low_red)),standard_error(na.omit(meta.c$D5_m_MP_Low_red)),standard_error(na.omit(meta.c$D3_m_MP_High_red)),standard_error(na.omit(meta.c$D5_m_MP_High_red)))
actrep.m<-as.data.frame(cbind(day,high.low,med,se),stringsAsfactors=F)
actrep.m$med<-round(as.numeric(actrep.m$med),0)
actrep.m$se<-as.numeric(actrep.m$se)

png(file="~/Schlumpf flow/20220415/publish/FgureS3_MFIxTimeSE.png",
    width=1500, height=1500, res=300)
ggplot(actrep.m, aes(day, med)) +
  geom_line(aes(color = high.low, group = high.low))+
  geom_point() +
  geom_errorbar(
   aes(ymin = med-se, ymax = med+se, group = high.low,color = high.low),
   width = 0.2) +
  labs(y= "MFI", x = "Day")  + scale_color_manual(values=c("grey","orange"))+ theme_bw() + theme(axis.text.y=element_text(angle=90)) +  scale_y_continuous(labels = scientific)
dev.off()

#Stats test for group comparison, ither t.test or wilcox test depending on distribution
#Day 3 active vs. persisters
t.test(meta.c$D3_m_MP_High_red, meta.c$D3_m_MP_Low_red,alternative = c("two.sided"),paired=F)
wilcox.test(meta.c$D3_m_MP_High_red, meta.c$D3_m_MP_Low_red,alternative = c("two.sided"),paired=F)

#Day 5 active vs. persisters 
t.test(meta.c$D5_m_MP_High_red, meta.c$D5_m_MP_Low_red,alternative = c("two.sided"),paired=F)
wilcox.test(meta.c$D5_m_MP_High_red, meta.c$D5_m_MP_Low_red,alternative = c("two.sided"),paired=F)

#Day 3 active vs. Day 5 active
t.test(meta.c$D3_Low_red, meta.c$D5_Low_red,alternative = c("two.sided"),paired=F)
wilcox.test(meta.c$D3_Low_red, meta.c$D5_Low_red,alternative = c("two.sided"),paired=F)

#Day 3 persisters vs. Day 5 persisters
t.test(meta.c$D3_m_MP_High_red, meta.c$D5_m_MP_High_red,alternative = c("two.sided"),paired=F)
wilcox.test(meta.c$D3_m_MP_High_red, meta.c$D5_m_MP_High_red,alternative = c("two.sided"),paired=F)

#Boxplot
png(file="D3 and D5 active and persistr bacteria per ml.png",
    width=1000, height=1000)

boxplot(meta.c$D3_Persisters_High_Red, meta.c$D3_Low_red,
        meta.c$D5_Persisters_High_Red, meta.c$D5_Low_red,
        xlab="Day", ylab="Bacteria/ml",ylim=c(0,1e6),
        names=c("D3","D3","D5","D5"),at = c(1,2,3,4) ,col=c("red", "orange","red", "orange"))
legend("topleft", c("High Red" ,"Low Red"), col=c("red","orange"), pch=c(15,15))
dev.off()
