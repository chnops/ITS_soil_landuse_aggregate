#Elizabeth Bach
#Diversity stats for Fan curated ITS data
#24 December 2014

#removes all objects previously assigned in R session, ensures no conflicting objects from another analysis
rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(lme4)
library(lmerTest)
library(bbmle)
library(reshape)

library(vegan)
library(ggplot2)
library(plyr)

#use "COBS_ITS_data_rar.csv", generated in COBS_ITS_rarifying.R
data.nosing.rar<-read.csv(file.choose())
head(data.nosing.rar[,1:10])
str(data.nosing.rar)

#calculating richness, shannons, and evenness
#for diversity measures, we will use Fisher's alpha for "richness", which is not the same as the sum of OTUs, see taxonomy_merge.R for sum of OTUs per sample
richness<-fisher.alpha(data.nosing.rar[,-c(1:12)],1)
head(richness)
shannons<-diversity(data.nosing.rar[,-c(1:12)])
evenness<-shannons/log(richness)
hist(richness)
div_stats<-data.frame(data.nosing.rar[,1:12],richness,shannons,evenness)

head(div_stats)
str(div_stats)

#looking at data distribution
ggplot(div_stats)+geom_histogram(aes(shannons))
ggplot(div_stats)+geom_histogram(aes(richness))
#skewed left
ggplot(div_stats)+geom_histogram(aes(evenness))

#richness is skewed heavily to left
ggplot(div_stats)+geom_histogram(aes(log(richness)))
#log transformation helps normalize data, still a tail at high end

#Testing main effects of date, crop, soil frac on diversity measures, linear model
summary(test<-aov(richness~Date+Crop+SoilFrac, data=div_stats))
TukeyHSD(test)
#looking for significant interactions
summary(test<-aov(richness~Date*Crop*SoilFrac, data=div_stats))
#Date*Crop and Date*SoilFrac only significant
summary(test<-aov(richness~Date+Crop+SoilFrac+Date*Crop+Date*SoilFrac, data=div_stats))

#Using mixed model to look with block effect, no diff
summary(test4<-lmer(richness~Date+Crop+SoilFrac+(1|block), data=div_stats, REML=FALSE))

#Looking with nested agg + random block effect, all interactions:
test2<-lmer(richness~Date*Crop*SoilFrac+(1|block)+(1|Date/Crop/SoilFrac), data=div_stats, REML=FALSE)
summary(test2)
#Crop and Date significance drops out

#Looking with nested agg, just Date*Crop interaction
test3<-lmer(richness~Date+Crop+SoilFrac+Date:Crop+(1|block)+(1|Date/Crop/SoilFrac), data=div_stats, REML=FALSE)
summary(test3)
#Crop and agg signifcant, date NS

anova(test2, test3, test4)
#AIC values all within a few digits of each other (and fairly high)
#Proceed with basic ANOVA with Date*Crop and Date*SoilFrac interactions, simplest model
#No radical difference in significant results between models

summary(test)
TukeyHSD(test)
#Date, Crop, SoilFrac and interactions are significant, look at data to understand these differences

#Making quick graphs to look at differences among factors
#Requires R. Williams bootstrap function
#Bootstrap functions from R. Williams to generate confidence intervals (not standard error)
boot.high<-function(XX){
boot.mean<-numeric(1000)
for (i in 1:1000){
 boot.mean[i]<-mean(sample(XX,replace=T))
}
return(quantile(boot.mean,(0.975)))
}

boot.low<-function(XX){
boot.mean<-numeric(1000)
for (i in 1:1000){
 boot.mean[i]<-mean(sample(XX,replace=T))
}
return(quantile(boot.mean,(0.025)))
}

#Inital look at management system differences in richness
Rich<-ddply(div_stats, .(Crop), summarise,.progress="text",
mean=mean(richness),
high95=boot.high(richness),
low95=boot.low(richness)
)
Rich
ggplot(Rich, aes(Crop, mean))+geom_pointrange(aes(ymax=high95, ymin=low95))

#Now look at SoilFrac differences in richness
Rich2<-ddply(div_stats, .(SoilFrac), summarise,.progress="text",
mean=mean(richness),
high95=boot.high(richness),
low95=boot.low(richness)
)
Rich2
ggplot(Rich2, aes(SoilFrac, mean))+geom_pointrange(aes(ymax=high95, ymin=low95))


#Evenness
summary(test<-aov(evenness~Date+Crop+SoilFrac, data=div_stats))
summary(test2<-aov(evenness~Date*Crop*SoilFrac, data=div_stats))
#Date*SoilFrac significant
summary(test3<-aov(evenness~Date+Crop+SoilFrac+Date*SoilFrac, data=div_stats))
TukeyHSD(test3)
#Date, Crop, and Date*SoilFrac significant

Even<-ddply(div_stats, .(Crop), summarise,.progress="text",
mean=mean(evenness),
high95=boot.high(evenness),
low95=boot.low(evenness)
)
Even
ggplot(Even, aes(Crop, mean))+geom_pointrange(aes(ymax=high95, ymin=low95))
#Date*SoilFrac interaction
Even2<-ddply(div_stats, .(Date,SoilFrac), summarise,.progress="text",
mean=mean(evenness),
high95=boot.high(evenness),
low95=boot.low(evenness)
)
Even2
ggplot(Even2, aes(Date, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=SoilFrac, color=SoilFrac))

Even3<-ddply(div_stats, .(SoilFrac), summarise,.progress="text",
mean=mean(evenness),
high95=boot.high(evenness),
low95=boot.low(evenness)
)
Even3
ggplot(Even3, aes(SoilFrac, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=SoilFrac, color=SoilFrac))

#Shannon's diversity
summary(test<-aov(shannons~Date+Crop+SoilFrac, data=div_stats))
summary(test2<-aov(shannons~Date*Crop*SoilFrac, data=div_stats))
#Date*Crop significant interaction
summary(test3<-aov(shannons~Date+Crop+SoilFrac+Date*Crop, data=div_stats))
TukeyHSD(test3)
#Crop, SoilFrac, Date*Crop significant
#Date*Crop
Shan<-ddply(div_stats, .(Date,Crop), summarise,.progress="text",
mean=mean(shannons),
high95=boot.high(shannons),
low95=boot.low(shannons)
)
Shan
ggplot(Shan, aes(Date, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=Crop, color=Crop))

#Crop
Shan2<-ddply(div_stats, .(Crop), summarise,.progress="text",
mean=mean(shannons),
high95=boot.high(shannons),
low95=boot.low(shannons)
)
Shan2
ggplot(Shan2, aes(Crop, mean))+geom_pointrange(aes(ymax=high95, ymin=low95))

#SoilFrac
Shan3<-ddply(div_stats, .(SoilFrac), summarise,.progress="text",
mean=mean(shannons),
high95=boot.high(shannons),
low95=boot.low(shannons)
)
Shan3
ggplot(Shan3, aes(SoilFrac, mean))+geom_pointrange(aes(ymax=high95, ymin=low95))
