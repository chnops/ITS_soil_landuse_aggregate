#Elizabeth Bach
#Glomeromycota panel
#15 November 2016

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

#figure
ggplot(Glom.inter)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95, color=Date))+facet_wrap(~Crop, scales="free",ncol=3)+theme_bw()

#For aggregate fraction figure
Fung.agg<-ddply(data_taxa2, .(SoilFrac, phylum), summarise, .progress="text",mean=mean(value),
high95=boot.high(value), low95=boot.low(value))
levels(Fung.agg$SoilFrac)
Fung.agg$SoilFrac<-factor(Fung.agg$SoilFrac, levels=c("Micro","SM","MM","LM","WS"), order=TRUE)

ggplot(Fung.agg)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95), size=0.75)+theme_bw()+facet_wrap(~phylum, scales="free",ncol=3)+
theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1.5), strip.text=element_text(size=10, face="bold"), axis.text=element_text(size=12))
