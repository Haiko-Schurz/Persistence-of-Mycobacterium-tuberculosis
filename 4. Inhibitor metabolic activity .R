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

##########################################################################################################
######### Inhibitor metabolic
##########################################################################################################
inhibitor.m<-read.table("Inhibitor_metabolic_sorted.txt",header=T)

#Plot distribution of Day 3 nd 5 active and persister metabolic activity (MFI) for all inhibitors
png(file="~/Schlumpf flow/20220415/publish/Box plot D3 and D5 inhibitor MFI persister (grey) and active (orange).png",
    width=7500, height=5000, res=300)
par(mar=c(5,5,4,2))
boxplot(inhibitor.m$D3_Un_Persisters,inhibitor.m$D5_Un_Persisters,inhibitor.m$D3_Un_Replicating,inhibitor.m$D5_Un_Replicating,
        inhibitor.m$D3_Cyt_Persisters,inhibitor.m$D5_Cyt_Persisters,inhibitor.m$D3_Cyt_Replicating,inhibitor.m$D5_Cyt_Replicating,
        inhibitor.m$D3_Baf_Persisters,inhibitor.m$D5_Baf_Persisters,inhibitor.m$D3_Baf_Replicating,inhibitor.m$D5_Baf_Replicating,
        ylab="MFI", cex.axis = 2, cex.lab = 2,
        names=c("D3_Un","D5_Un","D3_Un","D5_Un","D3_CytD","D5_CytD","D3_CytD","D5_CytD","D3_BafA1","D5_BafA1","D3_BafA1","D5_BafA1"),at = c(1,2,3,4,5,6,7,8,9,10,11,12),
        col=c("grey","grey","orange","orange","grey","grey","orange","orange","grey","grey","orange","orange"),
        horizontal = F)
dev.off()

#Plot distribution of Day 3 nd 5 active and persister metabolic activity (MFI) for all inhibitors. Alternative layout
png(file="~/Schlumpf flow/20220415/publish/Box plot D3 and D5 inhibitor MFI persister (grey) and active (orange) Alternative.png",
    width=7500, height=5000, res=300)
par(mar=c(5,5,4,2))
boxplot(inhibitor.m$D3_Un_Persisters,inhibitor.m$D3_Un_Replicating,inhibitor.m$D5_Un_Persisters,inhibitor.m$D5_Un_Replicating,
        inhibitor.m$D3_Cyt_Persisters,inhibitor.m$D3_Cyt_Replicating,inhibitor.m$D5_Cyt_Persisters,inhibitor.m$D5_Cyt_Replicating,
        inhibitor.m$D3_Baf_Persisters,inhibitor.m$D3_Baf_Replicating,inhibitor.m$D5_Baf_Persisters,inhibitor.m$D5_Baf_Replicating,
        ylab="MFI", cex.axis = 2, cex.lab = 2,
        names=c("D3_Un","D3_Un","D5_Un","D5_Un","D3_CytD","D3_CytD","D5_CytD","D5_CytD","D3_BafA1","D3_BafA1","D5_BafA1","D5_BafA1"),at = c(1,2,3,4,5,6,7,8,9,10,11,12),
        col=c("grey","orange","grey","orange","grey","orange","grey","orange","grey","orange","grey","orange"),
        horizontal = F)
dev.off()

#Calculate values for line plot with error bars
day<-c(3,3,5,5,3,3,5,5,3,3,5,5)
group<-c("Un Persisters","Un Replicating","Un Persisters","Un Replicating","CytD Persisters","CytD Replicating","CytD Persisters","CytD Replicating","BafA1 Persisters","BafA1 Replicating","BafA1 Persisters","BafA1 Replicating")
med<-c(median(inhibitor.m$D3_Un_Persisters),median(inhibitor.m$D3_Un_Replicating),median(inhibitor.m$D5_Un_Persisters),median(inhibitor.m$D5_Un_Replicating),
       median(na.omit(inhibitor.m$D3_Cyt_Persisters)),median(na.omit(inhibitor.m$D3_Cyt_Replicating)),median(na.omit(inhibitor.m$D5_Cyt_Persisters)),median(na.omit(inhibitor.m$D5_Cyt_Replicating)),
       median(inhibitor.m$D3_Baf_Persisters),median(inhibitor.m$D3_Baf_Replicating),median(inhibitor.m$D5_Baf_Persisters),median(inhibitor.m$D5_Baf_Replicating))
se<-c(standard_error(inhibitor.m$D3_Un_Persisters),standard_error(inhibitor.m$D3_Un_Replicating),standard_error(inhibitor.m$D5_Un_Persisters),standard_error(inhibitor.m$D5_Un_Replicating),
      standard_error(na.omit(inhibitor.m$D3_Cyt_Persisters)),standard_error(na.omit(inhibitor.m$D3_Cyt_Replicating)),standard_error(na.omit(inhibitor.m$D5_Cyt_Persisters)),standard_error(na.omit(inhibitor.m$D5_Cyt_Replicating)),
      standard_error(inhibitor.m$D3_Baf_Persisters),standard_error(inhibitor.m$D3_Baf_Replicating),standard_error(inhibitor.m$D5_Baf_Persisters),standard_error(inhibitor.m$D5_Baf_Replicating))

line<-c("solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed")
line<-c("Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating")
inhib.pr<-as.data.frame(cbind(day,group,med,line,se))
inhib.pr$med<-as.numeric(inhib.pr$med)
inhib.pr$se<-as.numeric(inhib.pr$se)

#Figure 7A of publication. Line plot with error bars for mtabolic activity of all groups at day 3 and 5
png(file="~/Schlumpf flow/20220415/publish/Figure7_CFUxTimeSE.png",
    width=1500, height=1500, res=300)
ggplot(inhib.pr, aes(day, med)) +
  geom_line(aes(color = group, group = group,linetype=line))+
  geom_point() +
  geom_errorbar(
    aes(ymin = med-se, ymax = med+se, group = group,color = group,linetype=line),
    width = 0.2) +
  labs(y= "Bacteria/ml", x = "Day") + scale_color_manual(values=c("orange","orange","red","red","grey","grey")) + theme_bw() + theme(axis.text.y=element_text(angle=90)) +  scale_y_continuous(labels = scientific)

dev.off()

#Figure 7A of publication. distribution mtabolic activity of all groups at day 3
png(file="~/Schlumpf flow/20220415/publish/updated/FIG 7A Box plot D3 inhibitor MFI persister (grey) and active (orange).png",
    width=3500, height=2000, res=300)
par(mar=c(5,5,4,2))
boxplot(inhibitor.m$D3_Un_Persisters,inhibitor.m$D3_Un_Replicating,
        inhibitor.m$D3_Cyt_Persisters,inhibitor.m$D3_Cyt_Replicating,
        inhibitor.m$D3_Baf_Persisters,inhibitor.m$D3_Baf_Replicating,
        ylab="MFI", ylim = c(0,550),
        names=c("D3_Un","D3_Un","D3_CytD","D3_CytD","D3_BafA1","D3_BafA1"),at = c(1,2,3,4,5,6),
        col=c("grey","orange","grey","orange","grey","orange"),
        horizontal = F)
dev.off()

#Figure 7A of publication. distribution mtabolic activity of all groups at day 5
png(file="~/Schlumpf flow/20220415/publish/updated/FIG 7B Box plot D5 inhibitor MFI persister (grey) and active (orange).png",
    width=3500, height=2000, res=300)
par(mar=c(5,5,4,2))
boxplot(inhibitor.m$D5_Un_Persisters,inhibitor.m$D5_Un_Replicating,
        inhibitor.m$D5_Cyt_Persisters,inhibitor.m$D5_Cyt_Replicating,
        inhibitor.m$D5_Baf_Persisters,inhibitor.m$D5_Baf_Replicating,
        ylab="MFI", ylim = c(0,550),
        names=c("D5_Un","D5_Un","D5_CytD","D5_CytD","D5_BafA1","D5_BafA1"),at = c(1,2,3,4,5,6),
        col=c("grey","orange","grey","orange","grey","orange"),
        horizontal = F)
#legend("topright", c("Bac" ,"Persister","Active"), col=c("grey","red","orange"), pch=c(15,15,15))
dev.off()

#Statistical analysis for all groups at Day 3
d3anova<-inhibitor.m[,c(1,4,5,6,7,8,9)]
d3anova <- d3anova %>%
  gather(key = "Exp", value = "MFI", D3_Un_Persisters, D3_Un_Replicating,D3_Cyt_Persisters,D3_Cyt_Replicating,D3_Baf_Persisters,D3_Baf_Replicating) 
d3anova$rep<-c(seq(1,12),seq(1,12),seq(1,12),seq(1,12),seq(1,12),seq(1,12))
colnames(d3anova)<-c("Group","count","rep")
d3anova$rep<-as.factor(d3anova$rep)
d3anova$Group<-as.factor(d3anova$Group)
d3anova<-d3anova[-c(31:36,43:48),]

res.aov <- anova_test(data = d3anova, dv = count, wid = rep, within = Group)
get_anova_table(res.aov)

pwc <- d3anova %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc

write.table(pwc, file = "D3 Inhibitor metabolic ttest.txt", col.names = T, row.names = F, quote = F)

#Statistical analysis for all groups at Day 5
d5anova<-inhibitor.m[,c(1,10,11,12,13,14,15)]
d5anova <- d5anova %>%
  gather(key = "Exp", value = "MFI", D5_Un_Persisters, D5_Un_Replicating,D5_Cyt_Persisters,D5_Cyt_Replicating,D5_Baf_Persisters,D5_Baf_Replicating) 
d5anova$rep<-c(seq(1,12),seq(1,12),seq(1,12),seq(1,12),seq(1,12),seq(1,12))
colnames(d5anova)<-c("Group","count","rep")
d5anova$rep<-as.factor(d5anova$rep)
d5anova$Group<-as.factor(d5anova$Group)
d5anova<-d5anova[-c(31:36,43:48),]

res.aov <- anova_test(data = d5anova, dv = count, wid = rep, within = Group)
get_anova_table(res.aov)

pwc <- d5anova %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc

write.table(pwc, file = "D5 Inhibitor metabolic ttest.txt", col.names = T, row.names = F, quote = F)

#Statistical analysis for all groups between day 3 and day 5
d3anova$Day<-"Day3"
#colnames(d3anova)<-
d5anova$Day<-"Day5"
d3vd5<-as.data.frame(rbind(d3anova,d5anova))
d3vd5$Day<-as.factor(d3vd5$Day)

res.aov <- anova_test(
  data = d3vd5, dv = count, wid = rep,
  within = c(Day, Group)
)
get_anova_table(res.aov)
library(lme4)
model = lm(count ~Group*Day , data = d3vd5)
summary(model)

pwc <- d3vd5 %>%
  group_by(Group) %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc
print(pwc)
write.table(pwc, file = "d3vd5 Inhibitor meta ttest.txt", col.names = T, row.names = F, quote = F)

pwc<-pairwise_t_test(count ~ Group, paired = F,   p.adjust.method = "bonferroni",data=d3vd5)

