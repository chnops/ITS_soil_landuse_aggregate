#Elizabeth Bach
#COBS ITS
#Ascoomycota taxa analysis
#24 Jan 2015

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

#Subsetting data to look exclusively in Ascomycota, determin which class/order/families might be driving differences
Asco.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Ascomycota"))
#Analysis in "COBS_ITS_PhylaAnalysis.R" indicates Ascomycota affected by 3-way Crop*Date*SoilFrac interaction
Asco.3way<-ddply(Asco.data, .(Crop,Date,SoilFrac, class), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Asco.3way) 
str(Asco.3way)
ggplot(Asco.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Crop))+facet_wrap(~class, scales="free",ncol=3)+theme_bw()
max(Asco.3way$mean)
max(Asco.data$value)

#Most abundant Ascos
Asco.class<-ddply(Asco.data, .(class), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
str(Asco.class)
Class.rank<-Asco.class[order(-Asco.class$avg),]
Class.rank

Asco.order<-ddply(Asco.data, .(order), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
str(Asco.order)
Order.rank<-Asco.order[order(-Asco.order$avg),]
Order.rank

Asco.genus<-ddply(Asco.data, .(genus), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
str(Asco.genus)
Genus.rank<-Asco.genus[order(-Asco.genus$avg),]
Genus.rank

#Class Incertae_sedis
Incert.data<-droplevels(subset(Asco.data, Asco.data$class=="Ascomycota_Incertae sedis_Incertae sedis"))
Incert.sedis<-ddply(Incert.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Incert.sedis)
ggplot(Incert.sedis)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
#Maybe outlier for PF, SM (300 reads)
max(Incert.data$value)
Incert.data2<-subset(Incert.data, Incert.data$value<290)
Incert.sedis<-ddply(Incert.data2, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))

Incert.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Incert.data2, REML=FALSE)
Incert.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Incert.data2, REML=FALSE)
AICtab(Incert.model.full,Incert.model.main)
anova(Incert.model.main)
#Main model best fit, not interactions in full model
lsmeans(Incert.model.main)

#Class Dothideomycetes
Doth.data<-droplevels(subset(Asco.data, Asco.data$class=="Dothideomycetes"))
Doth.sedis<-ddply(Doth.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Doth.sedis)
ggplot(Doth.sedis)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
levels(Doth.sedis$order)
#5 orders, look at class-level relationship

Doth.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Doth.data, REML=FALSE)
Doth.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Doth.data, REML=FALSE)
AICtab(Doth.model.full,Doth.model.main)
#main model best fit, not significant interactions
anova(Doth.model.main)
lsmeans(Doth.model.main)

#Which order(s) driving this?
ggplot(Doth.sedis)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Crop))+facet_wrap(~order, scales="free",ncol=3)+theme_bw()

#Order Capnodiales
Cap.data<-droplevels(subset(Asco.data, Asco.data$order=="Capnodiales"))
Cap.sedis<-ddply(Cap.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Cap.sedis)
ggplot(Cap.sedis)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

Cap.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Cap.data, REML=FALSE)
Cap.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Cap.data, REML=FALSE)
AICtab(Cap.model.full,Cap.model.main)
#Main model best fit and no interactions 
anova(Cap.model.main)
difflsmeans(Cap.model.main)

#Order Dothideales
Dothi.data<-droplevels(subset(Asco.data, Asco.data$order=="Dothideales"))
Dothi.3way<-ddply(Dothi.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Dothi.3way)
ggplot(Dothi.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
#Inconsistent presence across fractions and sampling dates

Dothi.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Dothi.data, REML=FALSE)
Dothi.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Dothi.data, REML=FALSE)
AICtab(Dothi.model.full,Dothi.model.main)
anova(Dothi.model.full)
#Full model is best fit, interactions highly significant, but largely driven by uneven data representation

#Incertae sedis
Dothi.incert<-droplevels(subset(Asco.data, Asco.data$order=="Dothideomycetes_Incertae sedis_Incertae sedis"))
Incert.3way<-ddply(Dothi.incert, .(Crop,Date,SoilFrac, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
str(Incert.3way)
ggplot(Incert.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Crop))+facet_wrap(~species, scales="free",ncol=3)+theme_bw()

#Look at 2 most abundant species up close
Pseud.sp<-droplevels(subset(Asco.data, Asco.data$species=="Pseudeurotium sp BEA_2010"))
Pseud.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Pseud.sp, REML=FALSE)
Pseud.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Pseud.sp, REML=FALSE)
AICtab(Pseud.model.full,Pseud.model.main)
anova(Pseud.model.full)
anova(Pseud.model.main)
Pseud.3way<-ddply(Pseud.sp, .(Crop,Date,SoilFrac, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Pseud.3way)
ggplot(Pseud.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

Pseud.app<-droplevels(subset(Asco.data, Asco.data$species=="Pseudogymnoascus appendiculatus"))
Pseudapp.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Pseud.app, REML=FALSE)
Pseudapp.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Pseud.app, REML=FALSE)
AICtab(Pseudapp.model.full,Pseudapp.model.main)
anova(Pseudapp.model.full)
#3-way interaction significant
Pseud.3app<-ddply(Pseud.app, .(Crop,Date,SoilFrac, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Pseud.3app)
ggplot(Pseud.3app)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

#Order Pleosporales
Pleo.data<-droplevels(subset(Asco.data, Asco.data$order=="Pleosporales"))
Pleo.3way<-ddply(Pleo.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Pleo.3way)
ggplot(Pleo.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

Pleo.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Pleo.data, REML=FALSE)
Pleo.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Pleo.data, REML=FALSE)
AICtab(Pleo.model.full,Pleo.model.main)
#Main model a better fit and no interactions
anova(Pleo.model.main)
difflsmeans(Pleo.model.main)
Pleo.Crop<-ddply(Pleo.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Pleo.Crop

Pleo.Date<-ddply(Pleo.data, .(Date, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Pleo.Date

levels(Pleo.data$genus)
Pleo.genera<-ddply(Pleo.data, .(genus), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
Pleo.genera2<-Pleo.genera[order(-Pleo.genera$avg),]
Pleo.genera2[1:15,]
#Pyrenophora and Alternaria are the biggies

#Pyrenophora
Pyre.data<-droplevels(subset(Asco.data, Asco.data$genus=="Pyrenophora"))
Pyre.data
#Only 2 observations

#Alternaria
Alt.data<-droplevels(subset(Asco.data, Asco.data$genus=="Alternaria"))
str(Alt.data)
Alt.3way<-ddply(Alt.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Alt.3way)
ggplot(Alt.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
#Maybe outlier, PF_MM_Oct
max(Alt.data$value)
Alt.data2<-subset(Alt.data, Alt.data$value<3000)
max(Alt.data2$value)
Alt.3way<-ddply(Alt.data2, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Alt.3way)
ggplot(Alt.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
#Much improved with outlier removed

Alt.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Alt.data2, REML=FALSE)
Alt.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Alt.data2, REML=FALSE)
AICtab(Alt.model.full,Alt.model.main)
#Main model best fit and no significant interactions
anova(Alt.model.main)
difflsmeans(Alt.model.main)
Alt.Crop<-ddply(Alt.data2, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Alt.Crop
levels(Alt.data2$species)
#Only 1 speces with multiple synonyms/name options

#Unclassified Dothideomycetes
unDoth.data<-droplevels(subset(Asco.data, Asco.data$order=="unclassified_Dothideomycetes"))
str(unDoth.data)
unDoth.3way<-ddply(unDoth.data, .(Crop,Date,SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(unDoth.3way)
ggplot(unDoth.3way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

unDoth.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=unDoth.data, REML=FALSE)
unDoth.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=unDoth.data, REML=FALSE)
AICtab(unDoth.model.full,unDoth.model.main)
#Main model is best fit, no interactions, but no real main effects either
anova(unDoth.model.main)

#Saccharomycetes (looking at SoilFrac effect)
Sacc.data<-droplevels(subset(Asco.data, Asco.data$class=="Saccharomycetes"))
str(Sacc.data)

Sacc.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Sacc.data, REML=FALSE)
Sacc.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Sacc.data, REML=FALSE)
AICtab(Sacc.model.full,Sacc.model.main)
#main model best fit, no interactions
anova(Sacc.model.main)
difflsmeans(Sacc.model.main)

Sacc.SoilFrac<-ddply(Sacc.data, .(SoilFrac), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Sacc.SoilFrac

max(Sacc.data$value)
Sacc.data2<-subset(Sacc.data, Sacc.data$value<245)
max(Sacc.data2$value)

#Eurotiomycetes (looking at SoilFrac effect)
Euro.data<-droplevels(subset(Asco.data, Asco.data$class=="Eurotiomycetes"))
str(Euro.data)

Euro.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Euro.data, REML=FALSE)
Euro.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Euro.data, REML=FALSE)
Euro.model<-lmer(value~Date+Crop+SoilFrac+Crop*SoilFrac+(1|block), data=Euro.data, REML=FALSE)
AICtab(Euro.model.full,Euro.model.main,Euro.model)
#main model best fit, but there was a Crop*SoilFrac interaction, model including that interaction better fit
anova(Euro.model)
difflsmeans(Euro.model)

Euro.2way<-ddply(Euro.data, .(Crop,SoilFrac), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Euro.2way)
ggplot(Euro.2way)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()


