#Elizabeth Bach
#COBS ITS
#Glomeromycota taxa analysis
#26 Jan 2015


rm(list=ls())
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(ggplot2)
library(plyr)

#Use "COBS_ITS_data_taxa.csv" generated in "COBS_ITS_taxonomy_merge.R" file
data_taxa2<-read.csv(file.choose())
head(data_taxa2)
str(data_taxa2)

#Bootstrap function from R. Williams
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

Glom.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Glomeromycota"))
str(Glom.data)
#Date*SoilFrac and Crop*SoilFrac interactions significant
Glom.inter<-ddply(Glom.data, .(Crop,Date,SoilFrac, class), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Glom.inter) 
str(Glom.inter)
ggplot(Glom.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
levels(Glom.data$genus)
#looking at ranked abundance of genera
Glom.genera<-ddply(Glom.data2, .(genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Glom.genera[order(-Glom.genera$mean),]

#Diversispora
Diver.data<-droplevels(subset(Glom.data, Glom.data$genus=="Diversispora"))
str(Diver.data)
#26 observations, 3 species
Diver.inter<-ddply(Diver.data2, .(Crop,Date,SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Diver.inter) 
str(Diver.inter)
ggplot(Diver.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
max(Diver.data$value)
#may have outlier
Glom.data2<-subset(Glom.data, Glom.data$value<385)
#Removing outlier brings Diverspora to lower point

#unclassified Glomeraceae
unGlom.data<-droplevels(subset(Glom.data, Glom.data$genus=="unclassified_Glomeraceae"))
str(unGlom.data)
#129 observations
unGlom.inter<-ddply(unGlom.data, .(Crop,Date,SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(unGlom.inter) 
str(unGlom.inter)
ggplot(unGlom.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

unGlom.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=unGlom.data, REML=FALSE)
unGlom.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=unGlom.data, REML=FALSE)
AICtab(unGlom.model.full,unGlom.model.main)
#Main model best fit, not significant interactions
anova(unGlom.model.main)

unGlom.high<-subset(unGlom.data, unGlom.data$value>100)
unGlom.high
unGlom.data2<-subset(unGlom.data, unGlom.data$value<300)
#Not really an outlier, doesn't change story, so leave all values in

#Glomus
Glomus.data<-droplevels(subset(Glom.data, Glom.data$genus=="Glomus"))
str(Glomus.data)
#2098 observations of 10 species
Glomus.inter<-ddply(Glomus.data, .(Crop,Date,SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Glomus.inter) 
str(Glomus.inter)
ggplot(Glomus.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

Glomus.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Glomus.data, REML=FALSE)
Glomus.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Glomus.data, REML=FALSE)
AICtab(Glomus.model.full,Glomus.model.main)
anova(Glomus.model.full)
Glomus.model<-lmer(value~Date+Crop+SoilFrac+Crop*SoilFrac+(1|block), data=Glomus.data, REML=FALSE)
AICtab(Glomus.model.full,Glomus.model.main,Glomus.model)
anova(Glomus.model)

Glomus.species<-ddply(Glomus.data, .(species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Glomus.species[order(-Glomus.species$mean),]
#look into these species

#Paraglomus
Para.data<-droplevels(subset(Glom.data, Glom.data$genus=="Paraglomus"))
str(Para.data)
#97 observations, 3 species
Para.inter<-ddply(Para.data, .(Crop,Date,SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Para.inter) 
str(Para.inter)
ggplot(Para.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

Para.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Para.data, REML=FALSE)
Para.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Para.data, REML=FALSE)
AICtab(Para.model.full,Para.model.main)
anova(Para.model.main)
#Main model best fit, no interactions
difflsmeans(Para.model.main)
Para.Date<-ddply(Para.data, .(Date, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Para.Date

Para.SoilFrac<-ddply(Para.data, .(SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Para.SoilFrac
