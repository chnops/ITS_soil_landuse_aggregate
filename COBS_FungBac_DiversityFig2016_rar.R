#Elizabeth Bach
#Diversity figure for ITS & 16S COBS data
#4 May 2016

rm(list=ls())
library(plyr)
library(vegan)
library(reshape)
library(ggplot2)
library(gridExtra)

#Fungi
#use "COBS_ITS_data_nonrar_no_outlier.csv"
data.fungi<-read.csv(file.choose())
dim(data.fungi)
data.fungi[1:10,1:10]
#add metadata
data.metadata2<-read.csv(file.choose(),na.strings=".")
head(data.metadata2)
data.fungi2<-merge(data.metadata2[,1:6], data.fungi, by.x="Sample")
head(data.fungi2[,1:10])
str(data.fungi2)

#calculating richness, shannons, and evenness

#calculating richness, shannons, and evenness
#Fisher's alpha value for logrithmic series
richness<-fisher.alpha(data.fungi[,-c(1:12)],1)
head(richness)
alpha_richness<-fisher.alpha(data.fungi[,-c(1:12)],1)

#Total OTUs
richness.2<-rowSums(data.fungi2[,-c(1:6)]>0)
head(richness.2)
hist(richness.2)
total_otu<-rowSums(data.fungi[,-c(1:12)]>0)

shannons<-diversity(data.fungi[,-c(1:12)])
hist(shannons)
evenness<-shannons/log(richness)
hist(evenness)

#removing outliers really smooths out data distribution

#ENSpie function from Ryan Williams, scaling diversity, from Chase & Knight

ENSpie<-function(x,Sample){
Sample<-Sample
enspie<-numeric(nrow(x))
for (i in 1:nrow(x)){
enspie[i]<-1/sum((x[i,]/rowSums(x[i,]))^2)
}
return(data.frame(Sample,enspie))
}

#Testing the function, there was much tweaking between these steps to get it right
ENSpie.test1<-ENSpie(data.fungi2[1,-c(1:7)])
ENSpie.test1

ENSpie.test2<-ENSpie(data.fungi2[2,-c(1:7)])
ENSpie.test2

ENSpie.test1.2<-ENSpie(data.fungi2[1:2,-c(1:7)])
ENSpie.test1.2
#seems to work!!!

ENSPie.ITS<-ENSpie(data.fungi[,-c(1:12)],data.fungi[,3])
ENSPie.ITS

#look at distribution
hist(ENSPie.ITS[,2])

effective_species<-ENSPie.ITS[,2]

div_stats.fungi<-data.frame(data.fungi[,3:7],alpha_richness,total_otu,shannons,evenness,effective_species)
head(div_stats.fungi)
str(div_stats.fungi)

#PointRange Figure
#Bootstrap functions from R. Williams
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

#SoilFrac figure for Manuscript
head(div_stats.fungi)
Div.data.fungi<-melt(div_stats.fungi, id=c("SampleName","SoilFrac","Date","Crop","block"))
head(Div.data.fungi)
levels(Div.data.fungi$SoilFrac)

Div.sum.fungi<-ddply(Div.data.fungi, .(SoilFrac,variable), summarise,.progress="text",
mean=mean(value),se=sd(value)/(sqrt(length(value)-1)),
high95=boot.high(value),
low95=boot.low(value)
)
head(Div.sum.fungi)
str(Div.sum.fungi)

sizes<-c(">2000",">2000",">2000",">2000",">2000","<250","<250","<250","<250","<250","1000-2000","1000-2000","1000-2000","1000-2000","1000-2000","250-1000","250-1000","250-1000","250-1000","250-1000","Whole Soil","Whole Soil","Whole Soil","Whole Soil","Whole Soil","WSprop","WSprop","WSprop","WSprop","WSprop")
#add back into data with WSprop
#,"Whole Soil","Whole Soil","Whole Soil","WS prop","WS prop","WS prop")
#
Div.sum.fungi2<-data.frame(sizes, Div.sum.fungi)
head(Div.sum.fungi2)
print(levels(Div.sum.fungi2$sizes))
Div.sum.fungi2$sizes=factor(Div.sum.fungi2$sizes, levels(Div.sum.fungi2$sizes)[c(1,4,3,2,5,6)])
print(levels(Div.sum.fungi2$sizes))

#Use LM, MM, etc.
print(levels(Div.sum.fungi$SoilFrac))
Div.sum.fungi$SoilFrac=factor(Div.sum.fungi$SoilFrac, levels(Div.sum.fungi$SoilFrac)[c(2,4,3,1,5,6)])
print(levels(Div.sum.fungi$SoilFrac))


agg.div.fung<-ggplot(Div.sum.fungi)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95), size=1)+facet_wrap(~variable, scales="free",ncol=5)+theme_bw()+
theme(aspect.ratio=1,text=element_text(face=2, size=20), panel.grid=element_blank(), panel.border=element_rect(size=3, colour="black"), legend.position="none", axis.ticks=element_line(size=2), axis.text.x=element_text(size=18, face="bold", colour="black", angle=30, vjust=0.8, hjust=0.75), strip.background=element_blank(), strip.text=element_text(size=20, face="bold"), axis.title=element_blank())
agg.div.fung


#16S data
#use "data_rar_WSprop.csv"
data.16S<-read.csv(file.choose(), header=TRUE)
head(data.16S[,1:10])
str(data.16S)

#calculating richness, shannons, and evenness
#Fisher's alpha value for logrithmic series
richness<-fisher.alpha(data.16S[,-c(1:8)],1)
alpha_richness<-fisher.alpha(data.16S[,-c(1:8)],1)
head(richness)
#Total OTUs
richness.2<-rowSums(data.16S[,-c(1:8)]>0)
head(richness.2)
hist(richness.2)
total_otu<-rowSums(data.16S[,-c(1:8)]>0)

shannons<-diversity(data.16S[,-c(1:8)])
evenness<-shannons/log(richness)
hist(richness)

#ENSpie function from Ryan Williams, scaling diversity, from Chase & Knight

ENSpie<-function(x,Sample){
Sample<-Sample
enspie<-numeric(nrow(x))
for (i in 1:nrow(x)){
enspie[i]<-1/sum((x[i,]/rowSums(x[i,]))^2)
}
return(data.frame(Sample,enspie))
}


ENSPie.16S<-ENSpie(data.16S[,-c(1:8)],data.16S[,1])
ENSPie.16S

#look at distribution
hist(ENSPie.16S[,2])

effective_species<-ENSPie.16S[,2]

div_stats.16S<-data.frame(data.16S[,1:5],alpha_richness,total_otu,shannons,evenness,effective_species)

head(div_stats.16S)
str(div_stats.16S)

#PointRange Figure

#SoilFrac figure for Manuscript
head(div_stats.16S)
Div.data.16S<-melt(div_stats.16S, id=c("Sample","CropBlock","SoilFrac","Date","Crop"))
head(Div.data.16S)
str(Div.data.16S)

Div.sum.16S<-ddply(Div.data.16S, .(SoilFrac,variable), summarise,.progress="text",
mean=mean(value),se=sd(value)/(sqrt(length(value)-1)),
high95=boot.high(value),
low95=boot.low(value)
)
head(Div.sum.16S)
str(Div.sum.16S)

sizes<-c(">2000",">2000",">2000",">2000",">2000","<250","<250","<250","<250","<250","1000-2000","1000-2000","1000-2000","1000-2000","1000-2000","250-1000","250-1000","250-1000","250-1000","250-1000","Whole Soil","Whole Soil","Whole Soil","Whole Soil","Whole Soil","WS prop","WS prop","WS prop","WS prop","WS prop")
#,"Whole Soil","Whole Soil","Whole Soil","Whole Soil","Whole Soil")
#,"Whole Soil","Whole Soil","Whole Soil","WS prop","WS prop","WS prop")
Div.sum.16S2<-data.frame(sizes, Div.sum.16S)
head(Div.sum.16S2)
str(Div.sum.16S2)
print(levels(Div.sum.16S2$sizes))
#Div.sum.16S2$sizes=factor(Div.sum.16S2$sizes, levels(Div.sum.16S2$sizes)[c(1,4,3,2,5,6)])
#print(levels(Div.sum.16S2$sizes))

#WS, LM, MM, etc
Div.sum.16S$SoilFrac=factor(Div.sum.16S$SoilFrac, levels(Div.sum.16S$SoilFrac)[c(2,4,3,1,5,6)])
print(levels(Div.sum.16S$SoilFrac))

agg.div.16S<-ggplot(Div.sum.16S)+geom_pointrange(aes(x=SoilFrac,y=mean,ymin=low95,ymax=high95), size=1)+facet_wrap(~variable, scales="free",ncol=5)+theme_bw()+
theme(aspect.ratio=1,text=element_text(face=2, size=20), panel.grid=element_blank(), panel.border=element_rect(size=3, colour="black"), legend.position="none", axis.ticks=element_line(size=2), axis.text.x=element_blank(), strip.background=element_blank(), strip.text=element_text(size=20, face="bold"), axis.title=element_blank())
agg.div.16S

grid.arrange(agg.div.16S, agg.div.fung, ncol=1)
