#Elizabeth Bach
#Merging OTU table with taxonomy
#R Williams data
#16 November 2016
rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

#use "data_rar_WSprop.csv" in Williams_16Sdata
data.nosing.rar<-read.csv(file.choose(), header=TRUE, check.names=FALSE)
#Filter out WSprop
data.16S<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac!="WSprop"))
head(data.16S[,1:10])
str(data.16S[,1:10])
dim(data.16S)

data_melt<-melt(data.16S, id=c("CropBlock","SoilFrac","Date","Crop","Block","Sample","prop_agg_fraction","splitter"))
head(data_melt)
dim(data_melt)

#Use "taxonomy_COBS_16s2.csv" file
taxonomy<-read.csv(file.choose())
head(taxonomy)
head(data_melt)

data_taxa<-merge(data_melt,taxonomy,by="variable")
head(data_taxa)
dim(data_taxa)
#data_taxa2<-data.frame(data_taxa[ which(data_taxa$value>0),])
#head(data_taxa2)
#dim(data_taxa2)
write.csv(data_taxa, file="COBS_16S_rar_taxa.csv")

#Looking at data first
str(data_taxa2)
#5059 OTUs (matches rarification step!!)
#107 samples, also matches rarifaction step
#value is a character vector, need to change to numeric?
levels(data_taxa2$phylum)
#7 Phyla: Asco, Basidio, Chtryi, "Incertae sedis" (unknw phylogenic placement, but described),
"unclassified",Zygo

#OTUs per sample
sample.richness<-ddply(data_taxa2, .(Sample), summarise, .progress="text", OTUr=length(variable))
dim(sample.richness)
head(sample.richness)
min(sample.richness$OTUr)
max(sample.richness$OTUr)
