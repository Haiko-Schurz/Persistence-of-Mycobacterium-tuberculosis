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

#########################################################################################################################################
###  Filter outliers based on D0%_uptake (chose cutoff of 30%)
#########################################################################################################################################
uptake<-read.table("Counts20220415.txt", header = T)
uptake$D0_uptake<-uptake$D0_uptake*100
uptake$D1_uptake<-uptake$D1_uptake*100

#Plot raw D0%_uptake for different MOI

boxplot((uptake$D0_uptake[uptake$MOI==5]),(uptake$D0_uptake[uptake$MOI==10]), (uptake$D0_uptake[uptake$MOI==20]), main = "D0 MP uptake", ylab="% uptake",
        names=c("MOI 5","MOI 10","MOI 20"),at = c(1,2,3) )

boxplot((uptake$D1_uptake[uptake$MOI==5]),(uptake$D1_uptake[uptake$MOI==10]), (uptake$D1_uptake[uptake$MOI==20]), main = "D1 MP uptake", ylab="% uptake",
        names=c("MOI 5","MOI 10","MOI 20"),at = c(1,2,3) )

boxplot(uptake$D0_Bac[uptake$MOI==5],uptake$D0_Bac[uptake$MOI==10], uptake$D0_Bac[uptake$MOI==20], main = "D0 in vitro Bacteria", ylab="Bacteria/ml",
        names=c("MOI 5","MOI 10","MOI 20"),at = c(1,2,3) )

mean(na.omit(uptake$D0_uptake[uptake$MOI==5]))
mean(na.omit(uptake$D0_uptake[uptake$MOI==10]))
mean(na.omit(uptake$D0_uptake[uptake$MOI==20]))

#Filter out outliers at 30% uptake
u30<-uptake[uptake$D0_uptake<=30,]

b<-which(u30$D0_Bac >=500000)
u30<-u30[-b,]

#Plot Day 0 uptake %
boxplot(u30$D0_uptake[u30$MOI==5],u30$D0_uptake[u30$MOI==10], u30$D0_uptake[u30$MOI==20], main = "D0 %uptake QC", ylab="% uptake",
        names=c("MOI 5","MOI 10","MOI 20"),at = c(1,2,3) )

mean(na.omit(u30$D0_uptake[u30$MOI==5]))
mean(na.omit(u30$D0_uptake[u30$MOI==10]))
mean(na.omit(u30$D0_uptake[u30$MOI==20]))

#plot day 0 invitro bacterial numbers
boxplot(u30$D0_Bac[u30$MOI==5],u30$D0_Bac[u30$MOI==10], u30$D0_Bac[u30$MOI==20], main = "D0 in vitro Bacteria", ylab="Bacteria/ml",
        names=c("MOI 5","MOI 10","MOI 20"),at = c(1,2,3) )

###################################################################################################################
#Make line plots
##################################################################################################################
#Calculate point and error bars for plot
bact_d0_5<-as.numeric(u30$D0_Bac[u30$MOI==5])
bact_d0_5m<-mean(na.omit(bact_d0_5))
bact_d0_5sd<-sd(na.omit(bact_d0_5))

bact_d0_10<-as.numeric(u30$D0_Bac[u30$MOI==10])
bact_d0_10m<-mean(na.omit(bact_d0_10))
bact_d0_10sd<-sd(na.omit(bact_d0_10))

bact_d0_20<-as.numeric(u30$D0_Bac[u30$MOI==20])
bact_d0_20m<-mean(na.omit(bact_d0_20))
bact_d0_20sd<-sd(na.omit(bact_d0_20))

ym<-c(bact_d0_5m,bact_d0_10m,bact_d0_20m)
ysd<-c(bact_d0_5sd,bact_d0_10sd,bact_d0_20sd)
x<-c(5,10,20)

# Plot Day 0 bacterial numbers for each MOI
plot(x, ym, xlab="MOI", ylab="Bacteria/ml", pch=16, cex=2, type='b', ylim=c(0,500000), main="D0 Bact")
# Add error bars
arrows(x0=x, y0=ym-ysd, x1=x, y1=ym+ysd, code=3, angle=90, length=0.1)

################# Sorted numbers correlation analysis  ##############################################################################
#Read in and QC sorted data 
sort<-read.table("counts_sorted.txt", header = T)
sortqc<-sort[-c(28:33),]  #remove outliers

#Plot correlation between different groups and calculte correaltion statistics
#Seeded bacteria numbers vs macrophage ingested bacteria at day 0
plot(sortqc$D0_Bac,sortqc$D0_MP)
cor.test(sortqc$D0_MP,sortqc$D0_Bac)
myline.fit <- lm( sortqc$D0_MP ~sortqc$D0_Bac)
abline(myline.fit)
summary(myline.fit)
scatterplot(D0_MP ~ D0_Bac, data = sortqc)

#Seeded bacteria numbers vs persisters at day 3
plot(sortqc$D0_Bac,sortqc$D3_Persisters_HighRed, main = "D0 bact vs. D3 Persister")
cor.test(sortqc$D3_Persisters_HighRed,sortqc$D0_Bac)
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D0_Bac)
abline(myline.fit)
summary(myline.fit)
scatterplot(D3_Persisters_HighRed ~ D0_Bac, data = sortqc)

#Seeded bacteria numbers vs persisters at day 5
plot(sortqc$D0_Bac,sortqc$D5_Persisters_HighRed, main = "D0 bact vs. D5 Persister")
cor.test(sortqc$D5_Persisters_HighRed,sortqc$D0_Bac)
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D0_Bac)
abline(myline.fit)
summary(myline.fit)
scatterplot(D5_Persisters_HighRed ~ D0_Bac, data = sortqc)

# non-parametric regression of the bacterial numbers 
# To plot behind the data:

#Day 0 Bacteria numbers vs, Day 0,1,3,5 macrophage internal bacteria

#Day 0 bacteria vs Day 0 Macrophage bacteria
png(file="D0Bact vs D0MP.png",
    width=1500, height=1500, res=300)
l <- loess(D0_MP ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D0_MP, pch=19,xlab = "",ylab = "",xaxt='n',
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",xaxt='n')))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
myline.fit <- lm( sortqc$D0_MP ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()

summary(myline.fit)
cor.test(sortqc$D0_MP , sortqc$D0_Bac)

#Day 0 bacteria vs Day 1 Macrophage bacteria
l <- loess(D1_MP ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D1_MP, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D1_MP ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D1_MP , sortqc$D0_Bac)

#Day 0 bacteria vs Day 3 Macrophage bacteria
png(file="D0Bact vs D3MP.png",
    width=1500, height=1500, res=300)
l <- loess(D3_MP ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D3_MP, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n', 
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n')))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
myline.fit <- lm( sortqc$D3_MP ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_MP , sortqc$D0_Bac)
#axis(1, at = seq(2E6, 10E6, 2E6), labels = c("2^6", "4^6", "6^6", "8^6", "10^6"))

#Day 0 bacteria vs Day 5 Macrophage bacteria
png(file="D0Bact vs D5MP.png",
    width=1500, height=1500, res=300)
l <- loess(D5_MP ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D5_MP, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n')))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
myline.fit <- lm( sortqc$D5_MP ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_MP , sortqc$D0_Bac)
method = c("pearson", "kendall", "spearman")
cor.test(sortqc$D5_MP , sortqc$D0_Bac,method = c( "spearman"))
cor.test(sortqc$D5_MP , sortqc$D0_Bac,method = c( "kendall"))


#Day 0 Bacteria numbers vs. Day 0 and Day 1 % Uptake
#Calculate % uptake and set missing to NA
sortqc$D0_uptake<-(sortqc$D0_MP/sortqc$D0_Bac)*100
sortqc$D1_uptake<-(sortqc$D1_MP/sortqc$D0_Bac)*100
sortqc$D0_uptake[c(1,2)]<-NA
sortqc$D1_uptake[c(1,2)]<-NA

#Day 0 bacteria vs Day 0 % uptake
l <- loess(D0_uptake ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D0_uptake, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D0_uptake ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D0_uptake , sortqc$D0_Bac)

#Day 0 bacteria vs Day 1 % uptake
l <- loess(D1_uptake ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D1_uptake, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D1_uptake ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D1_uptake , sortqc$D0_Bac)

boxplot(sortqc$D0_uptake,sortqc$D1_uptake, ylab="Uptake %",
        names=c("Day 0","Day 1"),at = c(1,2) )

mean(as.numeric(na.omit(sortqc$D0_uptake)))
sd(as.numeric(na.omit(sortqc$D0_uptake)))
mean(as.numeric(na.omit(sortqc$D1_uptake)))
sd(as.numeric(na.omit(sortqc$D1_uptake)))


###  Day 0 Bacteria vs. Day 3 and 5 persister numbers (High red) correlation analysis

#Day 0 Bacteria vs. Day 3 persister numbers
png(file="FIG 1C D0Bact vs D3HighRed.png",
    width=1500, height=1500, res=300)
l <- loess(D3_Persisters_HighRed ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D3_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE, xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_Persisters_HighRed , sortqc$D0_Bac)

#Day 0 Bacteria vs. Day 5 persister numbers
png(file="FIG 1C D0Bact vs D5HighRed.png",
    width=1500, height=1500, res=300)
l <- loess(D5_Persisters_HighRed ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D5_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_Persisters_HighRed , sortqc$D0_Bac)

###  Day 0 Bacteria vs. Day 3 and 5 replicating numbers (low red) correlation analysis

# Day 0 Bacteria vs. Day 3 replicating
png(file="FIG 1D D0Bact vs D3LowRed.png",
    width=1500, height=1500, res=300)
l <- loess(D3_Low_red ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D3_Low_red, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D3_Low_red ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_Low_red , sortqc$D0_Bac)

# Day 0 Bacteria vs. Day 5 replicating
png(file="FIG 1D D0Bact vs D5LowRed same scale.png",
    width=1500, height=1500, res=300)
l <- loess(D5_Low_red ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D5_Low_red, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e5, 1e6, 1.5e6), labels = c( "0e+00","5e+05", "1e+06", "1.5e+06"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D5_Low_red ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_Low_red , sortqc$D0_Bac)
cor.test(sortqc$D5_Low_red , sortqc$D0_Bac,method = c( "kendall"))

# Plot distribution between day 3 and 5 persisters (high red) and replicating (low red)
png(file="Box plot D3 and D5 active and persister numbers.png",
    width=1500, height=1500, res=300)
boxplot(sortqc$D3_Persisters_HighRed, sortqc$D3_Low_red, sortqc$D5_Persisters_HighRed,sortqc$D5_Low_red, yaxt='n', ylab="Bacteria/ml",
        names=c("Day 3","Day 3","Day 5","Day 5"),at = c(1,2,3,4) ,col=c("grey","orange","grey","orange"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
dev.off()


#########################################
days<-as.numeric(c(0,3,5,0,3,5,0,3,5))
group<-c(rep("Un",3),rep("CytD",3),rep("BafA1",3))
un.m<-c(median(inhibitor.c$D0_MP_Un),median(inhibitor.c$D3_MP_Un),median(inhibitor.c$D5_MP_Un))
cy.m<-c(median(na.omit(inhibitor.c$D0_MP_Cyt)),median(na.omit(inhibitor.c$D3_MP_Cyt)),median(na.omit(inhibitor.c$D5_MP_Cyt)))
baf.m<-c(median(inhibitor.c$D0_MP_Baf),median(inhibitor.c$D3_MP_Baf),median(inhibitor.c$D5_MP_Baf))
un.se<-c(standard_error(inhibitor.c$D0_MP_Un),standard_error(inhibitor.c$D3_MP_Un),standard_error(inhibitor.c$D5_MP_Un))
cy.se<-c(standard_error(na.omit(inhibitor.c$D0_MP_Cyt)),standard_error(na.omit(inhibitor.c$D3_MP_Cyt)),standard_error(na.omit(inhibitor.c$D5_MP_Cyt)))
baf.se<-c(standard_error(inhibitor.c$D0_MP_Baf),standard_error(inhibitor.c$D3_MP_Baf),standard_error(inhibitor.c$D5_MP_Baf))
med<-c(un.m,cy.m,baf.m)
se<-c(un.se,cy.se,baf.se)

day<-c(3,5,3,5)
high.low<-c("Low red","Low red","High red","High red")
med<-c(median(sortqc$D3_Low_red),median(sortqc$D5_Low_red),median(sortqc$D3_Persisters_HighRed),median(sortqc$D5_Persisters_HighRed))
se<-c(standard_error(sortqc$D3_Low_red),standard_error(sortqc$D5_Low_red),standard_error(sortqc$D3_Persisters_HighRed),standard_error(sortqc$D5_Persisters_HighRed))
actrep.c<-as.data.frame(cbind(day,high.low,med,se),stringsAsfactors=F)
actrep.c$med<-round(as.numeric(actrep.c$med),0)
actrep.c$se<-as.numeric(actrep.c$se)



png(file="~/Schlumpf flow/20220415/publish/Fgure3_CFUxTimeSE.png",
    width=1500, height=1500, res=300)
ggplot(actrep.c, aes(day, med)) +
  geom_line(aes(color = high.low, group = high.low))+
  geom_point() +
  geom_errorbar(
      aes(ymin = med-se, ymax = med+se, group = high.low,color = high.low),
      width = 0.2) +
  labs(y= "Bacteria/ml", x = "Day")  + scale_color_manual(values=c("grey","orange"))+ theme_bw() + theme(axis.text.y=element_text(angle=90)) + scale_y_continuous(labels = scientific)
dev.off()

+
  scale_y_continuous(name="Bacteria/ml", limits=c(1500, 3000),labels=c("1500" = "1.5e+03", "2000" = "2.0e+03", "2500" = "2.5e+03", "3000" = "3.0e+03"))
ggplot(data.inh.c, aes(days, med)) +
  geom_line(aes(color = group, group = group))+
  geom_point() + theme_bw() + theme(axis.text.y=element_text(angle=90)) +
  geom_errorbar(
    aes(ymin = med-se, ymax = med+se, group = group,color = group),
    width = 0.2) +
  labs(y= "Bacteria/ml", x = "Days")  + scale_color_manual(values=c("orange","red","grey"))

##########################################################################
### Intracellular (MP) bacteria vs. Persisters and active numbers at different days
# Correlation analysis

#   Day 0 Intracellular (MP) bacteria vs. Day 3 Persisters 
l <- loess(D3_Persisters_HighRed ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D3_Persisters_HighRed, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_Persisters_HighRed , sortqc$D0_MP)

#   Day 1 Intracellular (MP) bacteria vs. Day 3 Persisters 
l <- loess(D3_Persisters_HighRed ~ D1_MP, data = sortqc)
with(sortqc, plot(D1_MP, D3_Persisters_HighRed, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D1_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_Persisters_HighRed , sortqc$D1_MP)

#   Day 0 Intracellular (MP) bacteria vs. Day 3 active
l <- loess(D3_Low_red ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D3_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_Low_red ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_Low_red , sortqc$D0_MP)

#   Day 1 Intracellular (MP) bacteria vs. Day 3 active
l <- loess(D3_Low_red ~ D1_MP, data = sortqc)
with(sortqc, plot(D1_MP, D3_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_Low_red ~sortqc$D1_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_Low_red , sortqc$D1_MP)

#   Day 0 Intracellular (MP) bacteria vs. Day 5 Persisters
l <- loess(D5_Persisters_HighRed ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D5_Persisters_HighRed, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_Persisters_HighRed , sortqc$D0_MP)

#   Day 1 Intracellular (MP) bacteria vs. Day 5 Persisters
l <- loess(D5_Persisters_HighRed ~ D1_MP, data = sortqc)
with(sortqc, plot(D1_MP, D5_Persisters_HighRed, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D1_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_Persisters_HighRed , sortqc$D1_MP)

#   Day 0 Intracellular (MP) bacteria vs. Day 5 active
l <- loess(D5_Low_red ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D5_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_Low_red ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_Low_red , sortqc$D0_MP)

#   Day 1 Intracellular (MP) bacteria vs. Day 5 active
l <- loess(D5_Low_red ~ D1_MP, data = sortqc)
with(sortqc, plot(D1_MP, D5_Low_red, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_Low_red ~sortqc$D1_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_Low_red , sortqc$D1_MP)

#Plot distribution of different bacterial numbers
boxplot(sortqc$D0_MP, sortqc$D1_MP, sortqc$D3_Persisters_HighRed, sortqc$D3_Low_red, sortqc$D5_Persisters_HighRed,sortqc$D5_Low_red, ylab="Bacteria/ml",
        names=c("D0 MP", "D1 MP", "Day 3","Day 3","Day 5","Day 5"),at = c(1,2,3,4,5,6) ,col=c("green","green", "grey","orange","grey","orange"), ylim=c(0,2e5))
legend("topright", c("MP Bacteria" ,"Persisters","Active"), col=c("green", "gray","orange"), pch=c(15,15,15))

#D3 MP vs. D3 persistr/active and D5MP vs. D5 persister/active
#PLot figure 2 for publication 
png(file="FIG 2A replicating D3MP vs D3HighRed.png",
    width=1500, height=1500, res=300)
l <- loess(D3_Persisters_HighRed ~ D3_MP, data = sortqc)
with(sortqc, plot(D3_MP, D3_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D3_MP)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_Persisters_HighRed , sortqc$D3_MP)

png(file="FIG 2A persister D3MP vs D3LowRed.png",
    width=1500, height=1500, res=300)
l <- loess(D3_Low_red  ~ D3_MP, data = sortqc)
with(sortqc, plot(D3_MP, D3_Low_red , pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D3_Low_red  ~sortqc$D3_MP)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_Low_red , sortqc$D3_MP)

png(file="FIG 2B active D5MP vs D5HighRed.png",
    width=1500, height=1500, res=300)
l <- loess(D5_Persisters_HighRed ~ D5_MP, data = sortqc)
with(sortqc, plot(D5_MP, D5_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D5_MP)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_Persisters_HighRed , sortqc$D5_MP)

png(file="FIG 2B persister D5MP vs D5LowRed.png",
    width=1500, height=1500, res=300)
l <- loess(D5_Low_red  ~ D5_MP, data = sortqc)
with(sortqc, plot(D5_MP, D5_Low_red , pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7.2e4))))
axis(1, at = c(0e00, 5e4, 1e5, 1.5e5), labels = c( "0e+00","5e+04", "1e+05", "1.5e+05"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D5_Low_red  ~sortqc$D5_MP)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_Low_red , sortqc$D5_MP)

boxplot(sortqc$D3_MP, sortqc$D3_Persisters_HighRed, sortqc$D3_Low_red, sortqc$D5_MP, sortqc$D5_Persisters_HighRed,sortqc$D5_Low_red, ylab="Bacteria/ml",
        names=c("D3 MP", "Day 3","Day 3", "D5 MP","Day 5","Day 5"),at = c(1,2,3,4,5,6) ,col=c("red", "grey","orange","red","grey","orange"), ylim=c(0,2e5))
legend("topright", c("MP Bacteria" ,"Persisters","Active"), col=c("red", "gray","orange"), pch=c(15,15,15))


#Plot and correlation analysis of persisters and active for Day 3 and 5

#Day 3 persister vs. Day 3 active
png(file="D3persister vs D3active.png",
    width=1500, height=1500, res=300)
l <- loess(D3_Persisters_HighRed ~ D3_Low_red, data = sortqc)
with(sortqc, plot(D3_Low_red, D3_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7e4),
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n',ylim=c(0,7e4))))
axis(1, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D3_Persisters_HighRed ~sortqc$D3_Low_red)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D3_Persisters_HighRed , sortqc$D3_Low_red)

#Day 5 persister vs. Day 5 active
png(file="D5persister vs D5active.png",
    width=1500, height=1500, res=300)
l <- loess(D5_Persisters_HighRed ~ D5_Low_red, data = sortqc)
with(sortqc, plot(D5_Low_red, D5_Persisters_HighRed, pch=19,xlab = "",ylab = "",yaxt='n',xaxt='n',
                  panel.first=plot_loess(l, plotdata=FALSE,xlab = "",ylab = "",yaxt='n',xaxt='n')))
axis(1, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
axis(2, at = c(0e00, 2e4, 4e4, 6e4), labels = c( "0e+00","2e+04", "4e+04", "6e+04"))
myline.fit <- lm( sortqc$D5_Persisters_HighRed ~sortqc$D5_Low_red)
abline(myline.fit, col="blue")
dev.off()
summary(myline.fit)
cor.test(sortqc$D5_Persisters_HighRed , sortqc$D5_Low_red)

boxplot(sortqc$D0_Bac,sortqc$D0_MP, sortqc$D1_MP, sortqc$D3_MP, sortqc$D3_Persisters_HighRed, sortqc$D3_Low_red, sortqc$D5_MP, sortqc$D5_Persisters_HighRed,sortqc$D5_Low_red, ylab="Bacteria/ml",
        names=c("D0 Bact","D0 MP","D1 MP", "D3 MP", "Day 3","Day 3", "D5 MP","Day 5","Day 5"),at = c(1,2,3,4,5,6,7,8,9) ,col=c("green","red","red", "red", "grey","orange","red","grey","orange"), ylim=c(0,7e5))
legend("topright", c("Bacteria", "MP Bacteria" ,"Persisters","Active"), col=c("green", "red", "gray","orange"), pch=c(15,15,15,15))

###################### percent persister  #######################################
#Instead of numbers do correlation analysis using proportion of persisters relative to active numbers

#Day 0 bacteria vs. proportion of persisters at Day 3
l <- loess(D3_p_Persisters ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D3_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_p_Persisters ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_p_Persisters , sortqc$D0_Bac)

#Day 0 bacteria vs. proportion of persisters at Day 5
l <- loess(D5_p_Persisters ~ D0_Bac, data = sortqc)
with(sortqc, plot(D0_Bac, D5_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_p_Persisters ~sortqc$D0_Bac)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_p_Persisters , sortqc$D0_Bac)

#Day 0 macrophage bacteria vs. proportion of persisters at Day 3
l <- loess(D3_p_Persisters ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D3_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_p_Persisters ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_p_Persisters , sortqc$D0_MP)

#Day 0 macrophage bacteria vs. proportion of persisters at Day 5
l <- loess(D5_p_Persisters ~ D0_MP, data = sortqc)
with(sortqc, plot(D0_MP, D5_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D5_p_Persisters ~sortqc$D0_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_p_Persisters , sortqc$D0_MP)

#Day 3 macrophage bacteria vs. proportion of persisters at Day 3
l <- loess(D3_p_Persisters ~ D3_MP, data = sortqc)
with(sortqc, plot(D3_MP, D3_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
myline.fit <- lm( sortqc$D3_p_Persisters ~sortqc$D3_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D3_p_Persisters , sortqc$D3_MP)

#Day 5 macrophage bacteria vs. proportion of persisters at Day 5
l <- loess(D5_p_Persisters ~ D5_MP, data = sortqc)
with(sortqc, plot(D5_MP, D5_p_Persisters, pch=19,
                  panel.first=plot_loess(l, plotdata=FALSE)))
boxplot(sortqc$D0_Bac,sortqc$D0_MP, sortqc$D1_MP, sortqc$D3_MP, sortqc$D3_Persisters_HighRed, sortqc$D3_Low_red, sortqc$D5_MP, sortqc$D5_Persisters_HighRed,sortqc$D5_Low_red, ylab="Bacteria/ml",
        names=c("D0 Bact","D0 MP","D1 MP", "D3 MP", "Day 3","Day 3", "D5 MP","Day 5","Day 5"),at = c(1,2,3,4,5,6,7,8,9) ,col=c("green","red","red", "red", "grey","orange","red","grey","orange"), ylim=c(0,7e5))
legend("topright", c("Bacteria", "MP Bacteria" ,"Persisters","Active"), col=c("green", "red", "gray","orange"), pch=c(15,15,15,15))
myline.fit <- lm( sortqc$D5_p_Persisters ~sortqc$D5_MP)
abline(myline.fit, col="blue")
summary(myline.fit)
cor.test(sortqc$D5_p_Persisters , sortqc$D5_MP)

#Day 3 vs. Day 5 proportion of persisters 
png(file="D3 percent persister and D5 percent persister.png",
    width=1500, height=1500, res=300)
boxplot((sortqc$D3_p_Persisters)*100, (sortqc$D5_p_Persisters)*100,
        names=c("Day 3","Day 5"),at = c(1,2) ,col=c("grey","grey"))
#legend("topright", c("Persisters","Active"), col=c("gray","orange"), pch=c(15,15))
dev.off()
median((sortqc$D3_p_Persisters)*100)
median((sortqc$D5_p_Persisters)*100)
wilcox.test((sortqc$D3_p_Persisters)*100,(sortqc$D5_p_Persisters)*100)

