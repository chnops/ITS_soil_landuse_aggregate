#Elizabeth Bach
#COBS ITS
#Basidiomycota taxa analysis
#23 Jan 2015

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

Basidio.data<-droplevels(subset(data_taxa2, data_taxa2$phylum=="Basidiomycota"))
#Analysis in "COBS_ITS_PhylaAnalysis.R" indicates Basidiomycota only affected by Crop
Basidio.crop<-ddply(Basidio.data, .(Crop, class), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Basidio.crop) 
str(Basidio.crop)
ggplot(Basidio.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~class, scales="free",ncol=2)+theme_bw()

#Most abundant Basidiomycetes
#Excluding a couple of outliers
Basidio.data2<-subset(Basidio.data, Basidio.data$value<1000)
Basidio.abund<-ddply(Basidio.data, .(class), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
str(Basidio.abund)
Class.rank<-Basidio.abund[order(-Basidio.abund$avg),]
Class.rank

Basidio.order<-ddply(Basidio.data2, .(order), summarise, .progress="text",avg=mean(value),
high95=boot.high(value), low95=boot.low(value))
Order.rank<-Basidio.order[order(-Basidio.order$avg),]
Order.rank


#Class Agaricomycetes
Agari.data<-droplevels(subset(Basidio.data, Basidio.data$class=="Agaricomycetes"))
Agari.crop<-ddply(Agari.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Agari.crop)
ggplot(Agari.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~order, scales="free",ncol=2)+theme_bw()

#Order Agaricales
Agar.data<-droplevels(subset(Agari.data, Agari.data$order=="Agaricales"))
Agar.crop<-ddply(Agar.data, .(Crop, family), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Agar.crop)
ggplot(Agar.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~family, scales="free",ncol=2)+theme_bw()

#Family Entolomataceae
Ent.data<-droplevels(subset(Agari.data, Agari.data$family=="Entolomataceae"))
Ent.crop<-ddply(Ent.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Ent.crop)
ggplot(Ent.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
levels(Ent.data$species)
Ento.data<-droplevels(subset(Agari.data, Agari.data$genus=="Entoloma"))
Ento.crop<-ddply(Ento.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Ento.crop)
ggplot(Ento.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()
#Maybe an high outlier, removing one very high reading
Ento.data2<-droplevels(subset(Ento.data, Ento.data$value<1000))
str(Ento.data2)
Ento2.crop<-ddply(Ento.data2, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Ento2.crop)
ggplot(Ento2.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()
#Outlier was driving these relationhips

#Family Strophariaceae, all in genus Stropharia, species Stropharia rugosoannulata 
Stroph.data<-droplevels(subset(Agari.data, Agari.data$family=="Strophariaceae"))
Stroph.crop<-ddply(Stroph.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Stroph.crop)
ggplot(Stroph.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()
str(Stroph.data)
Stroph.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Stroph.data, REML=FALSE)
Stroph.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Stroph.data, REML=FALSE)
AICtab(Stroph.model.full,Stroph.model.main)
anova(Stroph.model.full)
anova(Stroph.model.main)
levels(Stroph.data$Date)
levels(Stroph.data$SoilFrac)

Stroph.SoilFrac<-ddply(Stroph.data, .(SoilFrac, Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Stroph.SoilFrac)
ggplot(Stroph.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

#Unclassified Agaricales
UnkAg.data<-droplevels(subset(Agari.data, Agari.data$order=="unclassified_Agaricomycetes"))
UnkAg.crop<-ddply(UnkAg.data, .(Crop, family), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(UnkAg.crop)
ggplot(UnkAg.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~family, scales="free",ncol=2)+theme_bw()

#Order Atheliales, all in family Atheliaceae
Athel.data<-droplevels(subset(Agari.data, Agari.data$order=="Atheliales"))
Athel.crop<-ddply(Athel.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Athel.crop)
ggplot(Athel.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
#Most reads in Cristinia
Crist.data<-droplevels(subset(Athel.data, Athel.data$genus=="Cristinia"))
Crist.data

#Order Auriculariales
Auri.data<-droplevels(subset(Agari.data, Agari.data$order=="Auriculariales"))
Auri.crop<-ddply(Auri.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Auri.crop)
ggplot(Auri.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
Auri.data2<-droplevels(subset(Auri.data, Auri.data$value<105))
str(Auri.data2)
Auri2.crop<-ddply(Auri.data2, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Auri2.crop)
ggplot(Auri2.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
Auri.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Auri.data2, REML=FALSE)
Auri.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Auri.data2, REML=FALSE)
AICtab(Auri.model.full,Auri.model.main)
anova(Auri.model.main)
lsmeans(Auri.model.main)
Auri2.SoilFrac<-ddply(Auri.data2, .(SoilFrac, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Auri2.SoilFrac)
ggplot(Auri2.SoilFrac)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()

#Order Cantharellales, all Ceratobasidiaceae
Canth.data<-droplevels(subset(Agari.data, Agari.data$order=="Cantharellales"))
Canth.crop<-ddply(Canth.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Canth.crop)
ggplot(Canth.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
#Test each genus
Cerato.data<-droplevels(subset(Agari.data, Agari.data$genus=="Ceratobasidium"))
Cerato.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Cerato.data, REML=FALSE)
Cerato.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Cerato.data, REML=FALSE)
AICtab(Cerato.model.full,Cerato.model.main)
anova(Cerato.model.main)

Cerato.crop<-ddply(Cerato.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Cerato.crop)
ggplot(Cerato.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

Than.data<-droplevels(subset(Agari.data, Agari.data$genus=="Thanatephorus"))
Than.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Than.data, REML=FALSE)
Than.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Than.data, REML=FALSE)
Than.model2<-lmer(value~Date+Crop+SoilFrac+Date*SoilFrac+Crop*SoilFrac+(1|block), data=Than.data, REML=FALSE)
AICtab(Than.model.full,Than.model.main,Than.model2)
anova(Than.model.full)

Than.crop<-ddply(Than.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Than.crop)
ggplot(Than.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

#Hymenochaetales is next order to investigate
Hymn.data<-droplevels(subset(Basidio.data, Basidio.data$order=="Hymenochaetales"))
Hymn.crop<-ddply(Hymn.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Hymn.crop)
ggplot(Hymn.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()

#Polyporales
Poly.data<-droplevels(subset(Basidio.data, Basidio.data$order=="Polyporales"))
Poly.crop<-ddply(Poly.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Poly.crop)
ggplot(Poly.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
#Ganoderma most abundant genus, one sample=544, no others >100
Gano.data<-droplevels(subset(Basidio.data, Basidio.data$genus=="Ganoderma"))
Gano.data2<-subset(Gano.data, Gano.data$value<100)
dim(Gano.data2)
#only 1 value greater than 100, excluding
Gano.crop<-ddply(Gano.data2, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Gano.crop)
ggplot(Gano.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

#Sebacinales
Seb.data<-droplevels(subset(Basidio.data, Basidio.data$order=="Sebacinales"))
Seb.crop<-ddply(Seb.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Seb.crop)
ggplot(Seb.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()

#Class Cystobasidiomycetes, only 1 order, Cystobasidiales, 1 genus Occultifur, 1 species Occultifu externus
Cysto.data<-droplevels(subset(Basidio.data, Basidio.data$class=="Cystobasidiomycetes"))
Cysto.crop<-ddply(Cysto.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Cysto.crop)
ggplot(Cysto.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

Cysto.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Cysto.data, REML=FALSE)
Cysto.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Cysto.data, REML=FALSE)
AICtab(Cysto.model.full,Cysto.model.main)
anova(Cysto.model.main)
#Model shows NS crop effect, but more abundant in Micros than other fractions
Cysto.SoilFrac<-ddply(Cysto.data, .(SoilFrac, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
Cysto.SoilFrac

#Class Microbotryomycetes, no apparent difference between crop
Micro.data<-droplevels(subset(Basidio.data, Basidio.data$class=="Microbotryomycetes"))
Micro.crop<-ddply(Micro.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Micro.crop)
ggplot(Micro.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~order, scales="free",ncol=2)+theme_bw()
#No crop differences at order level
Micro.crop<-ddply(Micro.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Micro.crop)
ggplot(Micro.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
subset(Micro.data, Micro.data$genus=="Rhodosporidium")

#Class Tremellomycetes, more abundant in CC
Treme.data<-droplevels(subset(Basidio.data, Basidio.data$class=="Tremellomycetes"))
Treme.crop<-ddply(Treme.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Treme.crop)
ggplot(Treme.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~order, scales="free",ncol=2)+theme_bw()
#consistent pattern across all orders, CC has more
#Filobasidiales is most abundant order, means 50-160
Filo.data<-droplevels(subset(Treme.data, Treme.data$order=="Filobasidiales"))
Filo.crop<-ddply(Filo.data, .(Crop, genus), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Filo.crop)
ggplot(Filo.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~genus, scales="free",ncol=2)+theme_bw()
#driven by Cryptococcus
subset(Filo.data, Filo.data$genus=="Cryptococcus g2")
Crypto.data<-droplevels(subset(Treme.data, Treme.data$genus=="Cryptococcus g2"))

Crypto.crop<-ddply(Crypto.data, .(Crop, species), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Crypto.crop)
ggplot(Crypto.crop)+geom_pointrange(aes(x=Crop,y=mean,ymin=low95,ymax=high95))+facet_wrap(~species, scales="free",ncol=2)+theme_bw()

Crypto.terr<-droplevels(subset(Treme.data, Treme.data$species=="Cryptococcus terreus"))
Crypto.model.full<-lmer(value~Date*Crop*SoilFrac+(1|block), data=Crypto.terr, REML=FALSE)
Crypto.model.main<-lmer(value~Date+Crop+SoilFrac+(1|block), data=Crypto.terr, REML=FALSE)
AICtab(Crypto.model.full,Crypto.model.main)
anova(Crypto.model.main)

#Unclassified
Unk.data<-droplevels(subset(Basidio.data, Basidio.data$class=="unclassified_Basidiomycota"))
Unk.crop<-ddply(Unk.data, .(Crop, order), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
head(Unk.crop)
#fairly low abundance
