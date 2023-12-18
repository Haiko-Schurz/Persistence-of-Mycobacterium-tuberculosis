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

#################################################################################################
#Inhibitors
#################################################################################################
inhibitor.c<-read.table("Inhibitor_counts_sortd.txt",header=T)
#Set missing values to NA
inhibitor.c$D0_Bac[12]<-NA
inhibitor.c$D3pPersisters.Un.[c(8:12)]<-NA

#Bacteria and macrophage bacteria counts
ggline(ToothGrowth, x = "dose", y = "len", 
       add = c("mean_sd", "jitter"),
       color = "supp", palette = c("#00AFBB", "#E7B800"))

#Get values for detailed plot
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

data.inh.c<-as.data.frame(cbind(days,group,med,se),stringsAsFactors = F)
data.inh.c$days<-as.numeric(data.inh.c$days)
data.inh.c$se<-as.numeric(data.inh.c$se)
data.inh.c$med<-as.numeric(data.inh.c$med)

#Plot macrophge bacterial numbers at the different days for the different treatments 
ggline(data.inh.c, x = "days", y = "med", 
       add = c("mean_sd", "jitter"),
       color = "group", palette = c("grey", "red","orange"))
c("grey","grey","grey", "red","red","red","orange","orange","orange")

#Plot firgure 4 for publication
png(file="Fgure4_CFUxTimeSE.png",
    width=1500, height=1500, res=300)
ggplot(data.inh.c, aes(days, med)) +
  geom_line(aes(color = group, group = group))+
  geom_point() + theme_bw() + theme(axis.text.y=element_text(angle=90)) +
  geom_errorbar(
    aes(ymin = med-se, ymax = med+se, group = group,color = group),
    width = 0.2) +
  labs(y= "Bacteria/ml", x = "Days")  + scale_color_manual(values=c("orange","red","grey"))
dev.off()

#Plot persister and active numbers for Day 3 and 5
ggplot(actrep.c, aes(day, med)) +
  geom_line(aes(color = high.low, group = high.low))+
  geom_point() +
  #  geom_errorbar(
  #    aes(ymin = med-se, ymax = med+se, group = group,color = group),
  #    width = 0.2) +
  labs(y= "Bacteria/ml", x = "Day")  + scale_color_manual(values=c("grey","orange"))+ theme_bw() + theme(axis.text.y=element_text(angle=90))+
  scale_y_continuous(name="Bacteria/ml", limits=c(1500, 3000),labels=c("1500" = "1.5e+03", "2000" = "2.0e+03", "2500" = "2.5e+03", "3000" = "3.0e+03"))
dev.off()

#plot distribution of bacterial macropahge numbers at Day 0
png(file="Box plot D0 inhibitor counts Un (grey) cyt (red) baf (orange) Un_SN (blue) Cyt_SN (maroon).png",
    width=1500, height=1500, res=300)
boxplot(inhibitor.c$D0_MP_Un,inhibitor.c$D0_MP_Cyt,inhibitor.c$D0_MP_Baf,inhibitor.c$D0_Un_SN,inhibitor.c$D0_Cyt_SN,
        ylab="Bacteria/ml", xlab = "Day 0",
        at = c(1,2,3,4,5) ,col=c("grey","red","orange","lightblue","brown"),
        horizontal = F)
dev.off()
#legend("topright", c("D0_Bac" ,"D0_MP_Un","D0_MP_Cyt","D0_MP_Baf","D0_Un_SN","D0_Cyt_SN"), col=c("grey","green","red","orange","blue","brown"), pch=c(15,15,15,15,15,15))

#plot distribution of bacterial macropahge numbers at Day 1
boxplot(inhibitor.c$D1_Bac,inhibitor.c$D1_MP_Un,inhibitor.c$D1_MP_Cyt,inhibitor.c$D1_MP_Baf,
        ylab="Bacteria/ml",
        names=c("D1","D1","D1","D1"),at = c(1,2,3,4) ,col=c("grey","green","red","orange"),
        horizontal = F,ylim=c(0,8e7))
legend("topright", c("D1_Bac" ,"D1_MP_Un","D1_MP_Cyt","D1_MP_Baf"), col=c("grey","green","red","orange"), pch=c(15,15,15,15))

#plot distribution of bacterial macropahge numbers at Day 3, for figure 4b in the publication
png(file="FIG 4B Box plot D3 inhibitor counts Un (grey) cyt (red) baf (orange).png",
    width=1500, height=1500, res=300)
boxplot(inhibitor.c$D3_MP_Un,inhibitor.c$D3_MP_Cyt,inhibitor.c$D3_MP_Baf,
        ylab="Bacteria/ml",xlab="Day 3",
        at = c(1,2,3) ,col=c("grey","red","orange"), ylim = c(0,8e7),
        horizontal = F)
dev.off()

#plot distribution of bacterial macropahge numbers at Day 5, for figure 4c in the publication
png(file="~/Schlumpf flow/20220415/publish/updated/FIG 4C Box plot D5 inhibitor counts Un (grey) cyt (red) baf (orange).png",
    width=1500, height=1500, res=300)
boxplot(inhibitor.c$D5_MP_Un,inhibitor.c$D5_MP_Cyt,inhibitor.c$D5_MP_Baf,
        ylab="Bacteria/ml",xlab="Day 5",
        at = c(1,2,3) ,col=c("grey","red","orange"), ylim = c(0,8e7),
        horizontal = F)
#legend("topright", c("D5_Bac" ,"D5_MP_Un","D5_MP_Cyt","D5_MP_Baf"), col=c("grey","green","red","orange"), pch=c(15,15,15,15))
dev.off()

#Plot distribution of persister and active numbers for all inhibitors at Day 3, for figure 5A in the publicaiton
png(file="FIG 5A Box plot D3 inhibitor counts for active (orange) and persisiter (grey).png",
    width=1600, height=1600, res=300)
boxplot(inhibitor.c$D3_Un_Persisters,inhibitor.c$D3_Un_Replicating,
        inhibitor.c$D3_Cyt_Persisters,inhibitor.c$D3_Cyt_Replicating,
        inhibitor.c$D3_Baf_Persisters,inhibitor.c$D3_Baf_Replicating,
        ylab="Bacteria/ml", xlab="Day 3", ylim = c(0,5e7),
        names=c("Un","Un","CytD","CytD","BafA1","BafA1"),at = c(1,2,3,4,5,6) ,col=c("grey","orange","grey","orange","grey","orange"),
        horizontal = F)
dev.off()

#Plot distribution of persister and active numbers for all inhibitors at Day 5, for figure 5A in the publication
png(file="FIG 5B Box plot D5 inhibitor counts for active (orange) and persisiter (grey).png",
    width=1600, height=1600, res=300)
boxplot(inhibitor.c$D5_Un_Persisters,inhibitor.c$D5_Un_Replicating,
        inhibitor.c$D5_Cyt_Persisters,inhibitor.c$D5_Cyt_Replicating,
        inhibitor.c$D5_Baf_Persisters,inhibitor.c$D5_Baf_Replicating,
        ylab="Bacteria/ml", xlab="Day 5",  ylim = c(0,5e7),
        names=c("Un","Un","CytD","CytD","BafA1","BafA1"),at = c(1,2,3,4,5,6) ,col=c("grey","orange","grey","orange","grey","orange"),
        horizontal = F)
dev.off()

#Plot distribution of active and persister for both Day 3 nd 5 across all treatments
png(file="Box plot D3 and D5 inhibitor counts for active (orange) and persisiter (grey).png",
    width=5800, height=3500, res=300)
boxplot(inhibitor.c$D3_Un_Persisters,inhibitor.c$D5_Un_Persisters,inhibitor.c$D3_Un_Replicating,inhibitor.c$D5_Un_Replicating,
        inhibitor.c$D3_Cyt_Persisters,inhibitor.c$D5_Cyt_Persisters,inhibitor.c$D3_Cyt_Replicating,inhibitor.c$D5_Cyt_Replicating,
        inhibitor.c$D3_Baf_Persisters,inhibitor.c$D5_Baf_Persisters,inhibitor.c$D3_Baf_Replicating,inhibitor.c$D5_Baf_Replicating,
        ylab="Bacteria/ml", cex.lab = 1.5, cex.axis = 1.5,
        names=c("D3_Un","D5_Un","D3_Un","D5_Un","D3_CytD","D5_CytD","D3_CytD","D5_CytD","D3_BafA1","D5_BafA1","D3_BafA1","D5_BafA1"),at = c(1,2,3,4,5,6,7,8,9,10,11,12),
        col=c("grey","grey","orange","orange","grey","grey","orange","orange","grey","grey","orange","orange"),
        horizontal = F)
#legend("topright", c("Bac" ,"MP","Persister","Active"), col=c("grey","green","red","orange"), pch=c(15,15,15,15))
dev.off()

#Plot distribution of active and persister for both Day 3 nd 5 across all treatments, alternative plot layout
png(file="Box plot D3 and D5 inhibitor counts for active (orange) and persisiter (grey) Alterntive.png",
    width=5800, height=3500, res=300)
boxplot(inhibitor.c$D3_Un_Persisters,inhibitor.c$D3_Un_Replicating,inhibitor.c$D5_Un_Persisters,inhibitor.c$D5_Un_Replicating,
        inhibitor.c$D3_Cyt_Persisters,inhibitor.c$D3_Cyt_Replicating,inhibitor.c$D5_Cyt_Persisters,inhibitor.c$D5_Cyt_Replicating,
        inhibitor.c$D3_Baf_Persisters,inhibitor.c$D3_Baf_Replicating,inhibitor.c$D5_Baf_Persisters,inhibitor.c$D5_Baf_Replicating,
        ylab="Bacteria/ml", cex.lab = 1.5, cex.axis = 1.5,
        names=c("D3_Un","D3_Un","D5_Un","D5_Un","D3_CytD","D3_CytD","D5_CytD","D5_CytD","D3_BafA1","D3_BafA1","D5_BafA1","D5_BafA1"),at = c(1,2,3,4,5,6,7,8,9,10,11,12),
        col=c("grey","orange","grey","orange","grey","orange","grey","orange","grey","orange","grey","orange"),
        horizontal = F)
dev.off()

#Calcluate values for the line plot (Figure 5 in publications)
day<-c(3,3,5,5,3,3,5,5,3,3,5,5)
group<-c("Un Persisters","Un Replicating","Un Persisters","Un Replicating","CytD Persisters","CytD Replicating","CytD Persisters","CytD Replicating","BafA1 Persisters","BafA1 Replicating","BafA1 Persisters","BafA1 Replicating")
med<-c(median(inhibitor.c$D3_Un_Persisters),median(inhibitor.c$D3_Un_Replicating),median(inhibitor.c$D5_Un_Persisters),median(inhibitor.c$D5_Un_Replicating),
       median(na.omit(inhibitor.c$D3_Cyt_Persisters)),median(na.omit(inhibitor.c$D3_Cyt_Replicating)),median(na.omit(inhibitor.c$D5_Cyt_Persisters)),median(na.omit(inhibitor.c$D5_Cyt_Replicating)),
       median(inhibitor.c$D3_Baf_Persisters),median(inhibitor.c$D3_Baf_Replicating),median(inhibitor.c$D5_Baf_Persisters),median(inhibitor.c$D5_Baf_Replicating))
se<-c(standard_error(inhibitor.c$D3_Un_Persisters),standard_error(inhibitor.c$D3_Un_Replicating),standard_error(inhibitor.c$D5_Un_Persisters),standard_error(inhibitor.c$D5_Un_Replicating),
      standard_error(na.omit(inhibitor.c$D3_Cyt_Persisters)),standard_error(na.omit(inhibitor.c$D3_Cyt_Replicating)),standard_error(na.omit(inhibitor.c$D5_Cyt_Persisters)),standard_error(na.omit(inhibitor.c$D5_Cyt_Replicating)),
      standard_error(inhibitor.c$D3_Baf_Persisters),standard_error(inhibitor.c$D3_Baf_Replicating),standard_error(inhibitor.c$D5_Baf_Persisters),standard_error(inhibitor.c$D5_Baf_Replicating))

line<-c("solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed")
line<-c("Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating","Persisters","replicating")
inhib.pr<-as.data.frame(cbind(day,group,med,line,se))
inhib.pr$med<-as.numeric(inhib.pr$med)
inhib.pr$se<-as.numeric(inhib.pr$se)

png(file="Figure5_CFUxTimeSE.png",
    width=1500, height=1500, res=300)
ggplot(inhib.pr, aes(day, med)) +
  geom_line(aes(color = group, group = group,linetype=line))+
  geom_point() +
  geom_errorbar(
     aes(ymin = med-se, ymax = med+se, group = group,color = group,linetype=line),
     width = 0.2) +
  labs(y= "Bacteria/ml", x = "Day") + scale_color_manual(values=c("orange","orange","red","red","grey","grey")) + theme_bw() + theme(axis.text.y=element_text(angle=90)) +  scale_y_continuous(labels = scientific)
  
dev.off()

#Plot the proportion of persisters instead of the actual numbers for all treatments at Day 3
png(file="Box plot D3 inhibitor percent persister.png",
    width=1600, height=1500, res=300)
boxplot(inhibitor.c$D3pPersisters.Un.,inhibitor.c$D3pPersisters.Cyt.,inhibitor.c$D3pPersisters.Baf.,
        ylab = "% Persister", names = c("Untreated", "CytD", "BafA1"), at = c(1,2,3), cex.lab = 1.2, cex.axis = 1.2, ylim=c(0,39))
dev.off()

#Plot the proportion of persisters instead of the actual numbers for all treatments at Day 3, Alternative layout
png(file="Box plot D3 inhibitor percent persister Alternative axes.png",
    width=1600, height=1500, res=300)
boxplot(inhibitor.c$D3pPersisters.Un.,inhibitor.c$D3pPersisters.Cyt.,inhibitor.c$D3pPersisters.Baf.,
        ylab = "% Persister", names = c("Untreated", "CytD", "BafA1"), at = c(1,2,3), cex.lab = 1.2, cex.axis = 1.2, ylim=c(10,40))
dev.off()

#Plot the proportion of persisters instead of the actual numbers for all treatments at Day 5
png(file="Box plot D5 inhibitor percent persister same scale.png",
    width=1600, height=1500, res=300)
boxplot(inhibitor.c$D5pPersisters.Un.,inhibitor.c$D5pPersisters.Cyt.,inhibitor.c$D5pPersisters.Baf.,
        ylab = "% Persister", names = c("Untreated", "CytD", "BafA1"), at = c(1,2,3),
        cex.lab = 1.2, cex.axis = 1.2, ylim=c(10,40))
dev.off()

#Stats test
#Anova btw Day 0/3/5 MP groups
#Day 0
d0anova<-inhibitor.c[,c(1,6,7,8)]
d0anova$rep<-d0anova$Exp
d0anova <- d0anova %>%
  gather(key = "Exp", value = "MFI", D0_MP_Un, D0_MP_Cyt, D0_MP_Baf) 
d0anova$rep<-c(seq(1,12),seq(1,12),seq(1,12))
colnames(d0anova)<-c("rep","Group","count")
d0anova$rep<-as.factor(d0anova$rep)
d0anova$Group<-as.factor(d0anova$Group)
d0anova<-d0anova[-c(19:24),]

res.aov <- anova_test(data = d0anova, dv = count, wid = rep, within = Group)
get_anova_table(res.aov)

pwc <- d0anova %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc
t.test(inhibitor.c$D0_Un_SN, inhibitor.c$D0_Cyt_SN , alternative = c("two.sided"),paired=F)
wilcox.test(inhibitor.c$D0_Un_SN, meta.c$D0_Cyt_SN, alternative = c("two.sided"),paired=F)

#Day 3
d3anova<-inhibitor.c[,c(1,14,15,16)]
d3anova$rep<-d3anova$Exp
d3anova <- d3anova %>%
  gather(key = "Exp", value = "MFI", D3_MP_Un, D3_MP_Cyt, D3_MP_Baf) 
d3anova$rep<-c(seq(1,12),seq(1,12),seq(1,12))
colnames(d3anova)<-c("rep","Group","count")
d3anova$rep<-as.factor(d3anova$rep)
d3anova$Group<-as.factor(d3anova$Group)
d3anova<-d3anova[-c(19:24),]

res.aov <- anova_test(data = d3anova, dv = count, wid = rep, within = Group)
get_anova_table(res.aov)

pwc <- d3anova %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc

#Day 5
d5anova<-inhibitor.c[,c(1,17,18,19)]
d5anova$rep<-d5anova$Exp
d5anova <- d5anova %>%
  gather(key = "Exp", value = "MFI", D5_MP_Un, D5_MP_Cyt, D5_MP_Baf) 
d5anova$rep<-c(seq(1,12),seq(1,12),seq(1,12))
colnames(d5anova)<-c("rep","Group","count")
d5anova$rep<-as.factor(d5anova$rep)
d5anova$Group<-as.factor(d5anova$Group)
d5anova<-d5anova[-c(19:24),]

res.aov <- anova_test(data = d5anova, dv = count, wid = rep, within = Group)
get_anova_table(res.aov)

pwc <- d5anova %>%
  pairwise_t_test(
    count ~ Group, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc

# Day 3/5 persisters and Day 3/5 active
#Day 3
d3anova<-inhibitor.c[,c(1,20,21,22,23,24,25)]
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

#Day 5
d5anova<-inhibitor.c[,c(1,26,27,28,29,30,31)]
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

# D3 persisters/active vs D5 persisters/active
d3anova$Day<-"Day3"
d5anova$Day<-"Day5"
d3vd5<-as.data.frame(rbind(d3anova,d5anova))

res.aov <- anova_test(
  data = d3vd5, dv = count, wid = rep,
  within = c(Day, Group)
)
get_anova_table(res.aov)
library(lme4)
model = lm(count ~Group*Day , data = d3vd5)
summary(model)

pwc <- d3vd5 %>%
  group_by(Day) %>%
  pairwise_t_test(
    count ~ Group*Day, paired = F,
    p.adjust.method = "bonferroni"
  )
pwc
print(pwc)
write.table(pwc, file = "d3vd5 Inhibitor count ttest.txt", col.names = T, row.names = F, quote = F)
pairwise_t_test(count ~ Group, paired = F,   p.adjust.method = "bonferroni",data=d3vd5)
