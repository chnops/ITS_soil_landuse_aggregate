#Elizabeth Bach
#COBS ITS data, Fan curated + Ryan rarificaiton for diversity
#Figure 3, WS prop comparison
#15 March 2015

rm(list=ls())
library(plyr)
library(vegan)
library(reshape)
library(ggplot2)

#use "phylo_and_fracs_unrar.csv"
data.16S<-read.csv(file.choose(), header=TRUE)
head(data.16S[,1:10])
str(data.16S)

#calculating reads and reads per sample
data.noWSprop<-droplevels(subset(data.16S, data.16S$SoilFrac!="WSprop"))
str(data.noWSprop[1:10])
reads<-rowSums(data.noWSprop[,-c(1:8)])
range(reads)

#total OTUs
richness.2<-rowSums(data.noWSprop[,-c(1:8)]>0)
head(richness.2)
hist(richness.2)
#double-checking above is a count, not sum
richness.3<-rowSums(data.noWSprop[,-c(1:8)])
head(richness.3)
#Yes, richness.2 is correct

range(richness.2)
data.OTUs<-cbind(data.noWSprop[,c(1:8)],richness.2)
head(data.OTUs)

#summarize by crop
crop.16S<-ddply(data.OTUs, .(Crop), summarize, mean=mean(richness.2))
crop.16S

#calculating richness, shannons, and evenness

richness<-fisher.alpha(data.matrix[,-c(1:5)],1)
head(richness)
shannons<-diversity(data.matrix[,-c(1:5)])
evenness<-shannons/log(richness)
hist(richness)
div_stats<-data.frame(data.matrix[,1:5],richness,shannons,evenness)

head(div_stats)
str(div_stats)

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
head(div_stats)
Div.data<-melt(div_stats, id=c("Sample","SoilFrac","Date","Crop","Block"))
head(Div.data)

Div.sum<-ddply(Div.data, .(SoilFrac,variable), summarise,.progress="text",
mean=mean(value),se=sd(value)/(sqrt(length(value)-1)),
high95=boot.high(value),
low95=boot.low(value)
)
head(Div.sum)
sizes<-c(">2000",">2000",">2000","<250","<250","<250","1000-2000","1000-2000","1000-2000","250-1000","250-1000","250-1000")
#,"Whole Soil","Whole Soil","Whole Soil","WS prop","WS prop","WS prop")
Div.sum2<-data.frame(sizes, Div.sum)
head(Div.sum2)
print(levels(Div.sum2$sizes))
Div.sum2$sizes=factor(Div.sum2$sizes, levels(Div.sum2$sizes)[c(1,4,3,2,5,6)])
print(levels(Div.sum2$sizes))

agg.div<-ggplot(Div.sum2)+geom_pointrange(aes(x=sizes,y=mean,ymin=low95,ymax=high95), size=1)+facet_wrap(~variable, scales="free",ncol=3)+theme_bw()+
theme(aspect.ratio=1,text=element_text(face=2, size=20), panel.grid=element_blank(), panel.border=element_rect(size=3, colour="black"), legend.position="none", axis.ticks=element_line(size=2), axis.text.x=element_text(size=18, face="bold", colour="black", angle=30, vjust=0.8, hjust=0.75), strip.background=element_blank(), strip.text=element_text(size=20, face="bold"), axis.title=element_blank())
agg.div
