#Elizabeth Bach
#COBS ITS analysis, phyla-levle analysis
#2 January 2015

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

#Phyla-level diversity, summarizes total counts for each phylum
Phyla.data2<-ddply(data_taxa2, .(Sample, Date, Crop, block, SoilFrac, phylum), summarise, .drop=FALSE, .progress="text",total=sum(value))
head(Phyla.data2)

Phyla.totals<-cast(data_taxa2, Sample + Date + Crop + block + SoilFrac ~ phylum, sum)
head(Phyla.totals)

#Looking at data, ggplot
#Williams bootstrapping function for confidence intervals
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

#Crop differences
Crop<-ddply(Phyla.data2, .(Crop, phylum), summarise, .progress="text",
mean=mean(total),high95=boot.high(total),low95=boot.low(total), N=length(total), SE=(sd(total)/sqrt(N-1)))
head(Crop)

ggplot(Crop, aes(phylum, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=Crop, color=Crop))
#separation between cropping systems evident for Ascomycota and unclassified fungi

#proportional distribution across crops
Crop.prop<-ddply(data_taxa2, .(Crop, phylum), summarise, .progress="text", reads=sum(value), prop=reads/1070000)
Crop.prop
sum(Crop.prop$reads)
sum(Crop.prop$prop)
Crop.sum<-ddply(Crop.prop, .(Crop), summarise, .progress="text", TotReads=sum(reads))
Crop.sum
Crop.prop2<-merge(Crop.prop, Crop.sum, by="Crop", all=TRUE)
Crop.prop2
Crop.prop3<-ddply(Crop.prop2, .(Crop, phylum), summarise, .progress="text", CropProp=reads/TotReads)
Crop.prop3
sum(Crop.prop3$CropProp)

Total.prop<-ddply(data_taxa2, .(phylum), summarise, .progress="text", reads=sum(value), prop=reads/1070000)
Total.prop
sum(Total.prop$prop)

#Aggregate Fraction differences
Agg<-ddply(Phyla.data2, .(SoilFrac, phylum), summarise, .progress="text",
mean=mean(total),high95=boot.high(total),low95=boot.low(total))
head(Agg)

ggplot(Agg, aes(phylum, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=SoilFrac, color=SoilFrac))
#Lots of overlap in data range, may not be any real differences in phyla-level distribution among fractions
#Glomeromycota maybe

#Date
Date<-ddply(Phyla.data2, .(Date, phylum), summarise, .progress="text",
mean=mean(total),high95=boot.high(total),low95=boot.low(total))
head(Date)

ggplot(Date, aes(phylum, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=Date, color=Date))
#Probably no significant differences between sampling dates, Asco has the greatest shift (more in Oct)

#Crop*Date
CropDate<-ddply(Phyla.data2, .(Crop, Date, phylum), summarise, .progress="text",
mean=mean(total),high95=boot.high(total),low95=boot.low(total))
head(CropDate)

ggplot(CropDate, aes(phylum, mean))+geom_pointrange(aes(ymax=high95, ymin=low95, group=Crop, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()
#changes between July and Oct seem consistent between cropping systmes within phyla

#Statistical models
#Ascomycota analysis
Asco.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Ascomycota"))
Asco.null<-lmer(value~1+(1|block), data=Asco.data, REML=FALSE)
Asco.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Asco.data, REML=FALSE)
Asco.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Asco.data, REML=FALSE)
AICtab(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.null,Asco.model.full,Asco.model.main)
anova(Asco.model.full)
#Full model has lowest AIC, and all interactions are significant, including Date*Crop*SoilFrac
#Visualizing shifts within Ascomycota
#Crop
Asco.crop<-ddply(Asco.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Asco.crop) 
str(Asco.crop)
levels(Asco.crop$order)
max(Asco.crop$mean)
min(Asco.crop$mean)
Asco.high<-droplevels(subset(Asco.crop, Asco.crop$mean>35))
str(Asco.high)
ggplot(Asco.high)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))
levels(Asco.high$order)
ggplot(Asco.data, aes(Crop, value, group=order))+geom_boxplot()
sorted.Asco<-Asco.crop[order(Asco.crop$Crop, -Asco.crop$mean),]
str(sorted.Asco)
head(sorted.Asco)
sorted.Asco

#SoilFrac
Asco.SoilFrac<-ddply(Asco.data, .(SoilFrac, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Asco.SoilFrac) 
str(Asco.SoilFrac)
SF.Asco<-Asco.SoilFrac[order(Asco.SoilFrac$SoilFrac, -Asco.SoilFrac$mean),]
str(SF.Asco)
ggplot(Asco.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95))
SF.Asco
#+facet_wrap(~order, scales="free",ncol=2)+theme_bw()


#Basidiomycota
Basidio.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Basidiomycota"))
Basidio.null<-lmer(value~1+(1|block), data=Basidio.data, REML=FALSE)
Basidio.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Basidio.data, REML=FALSE)
Basidio.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Basidio.data, REML=FALSE)
AICtab(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.null,Basidio.model.full,Basidio.model.main)
anova(Basidio.model.full)
#Main model best fit for AIC, no interactions are significant in full model, so proceed with main only
anova(Basidio.model.main)
difflsmeans(Basidio.model.main)
Basidio.cropAvg<-ddply(Basidio.data, .(Crop, phylum), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Basidio.cropAvg
#Diving into important orders
#Crop
Basidio.crop<-ddply(Basidio.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Basidio.crop) 
str(Basidio.crop)
levels(Basidio.crop$order)
ggplot(Basidio.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))
sorted.Basidio<-Basidio.crop[order(Basidio.crop$Crop, -Basidio.crop$mean),]
str(sorted.Basidio)
head(sorted.Basidio)
sorted.Basidio
Basidio.Corn<-droplevels(subset(Basidio.data, Basidio.data$Crop=="CC"))
Basidio.high<-droplevels(subset(Basidio.crop, Basidio.crop$mean>20))
ggplot(Basidio.high)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~order, scales="free",ncol=2)+theme_bw()

#Glomeromycota
Glom.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Glomeromycota"))
Glom.null<-lmer(value~1+(1|block), data=Glom.data, REML=FALSE)
Glom.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Glom.data, REML=FALSE)
Glom.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Glom.data, REML=FALSE)
AICtab(Glom.null,Glom.model.full,Glom.model.main)
anova(Glom.null,Glom.model.full,Glom.model.main)
anova(Glom.model.full)
#Full model best fit, only Date*SoilFrac and Crop*SoilFrac significant
Glom.model<-lmer(value~Date+Crop+SoilFrac+Date*SoilFrac+Crop*SoilFrac+(1|block), data=Glom.data, REML=FALSE)
AICtab(Glom.null,Glom.model.full,Glom.model.main,Glom.model)
anova(Glom.null,Glom.model.full,Glom.model.main,Glom.model)
#Glom.model best fit overall
anova(Glom.model)

#Zygomycota
Zygo.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Zygomycota"))
Zygo.null<-lmer(value~1+(1|block), data=Zygo.data, REML=FALSE)
Zygo.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Zygo.data, REML=FALSE)
Zygo.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Zygo.data, REML=FALSE)
AICtab(Zygo.null,Zygo.model.full,Zygo.model.main)
anova(Zygo.null,Zygo.model.full,Zygo.model.main)
#null model actually best fit, followed by main
anova(Zygo.model.full)
#No significant interactions
anova(Zygo.model.main)
#No significant factors

#Unk Fungi
Unk.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="unclassified_Fungi"))
Unk.null<-lmer(value~1+(1|block), data=Unk.data, REML=FALSE)
Unk.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Unk.data, REML=FALSE)
Unk.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Unk.data, REML=FALSE)
AICtab(Unk.null,Unk.model.full,Unk.model.main)
anova(Unk.null,Unk.model.full,Unk.model.main)
#Full model best fit
anova(Unk.model.full)
#Date*Crop, Date*SoilFrac, and Crop*SoilFrac interactions sig
Unk.model<-lmer(value~Date+Crop+SoilFrac+Date*Crop+Date*SoilFrac+SoilFrac*Crop+(1|block), data=Unk.data, REML=FALSE)
AICtab(Unk.null,Unk.model.full,Unk.model.main, Unk.model)
anova(Unk.null,Unk.model.full,Unk.model.main, Unk.model)
#Unk model best fit
anova(Unk.model)

#Chytridiomycota
Chytri.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Chytridiomycota"))
str(Chytri.data)
#Only 2 OTUs, found only in P13_SM_July, P13_micro_Oct, P13_SM_July, P13_WS_July, and P24_LM_Oct