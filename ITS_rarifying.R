#Elizabeth Bach
#Fan-curated ITS sequence data,
#rarification
#16 December 2014

#removes all objects previously assigned in R session, ensures no conflicting objects from another analysis
rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

#use "COBS_ITS_rawFan.csv" for rarification
dataset<-read.csv(file.choose())
dataset[1:10,1:10]
dim(dataset)
str(dataset)

# First we remove the singletons using the dropspc() function form the labdsv package.  In the line below I bind metadata (columns 1 through 7)
# and the singleton-removed dataset.  Note that for the dropspc function I only use columns 8 through 7437 as these are all the OTUs (no metadata),
# and I remove species that occur 1 or less times

data.nosing<-cbind(dataset[,c(1:7)],dropspc(dataset[,8:7437],1))
str(data.nosing)
dim(data.nosing)
data.nosing[1:10,1:10]

#Removing singltons brings OTUs to 5057

#Now its time to figure out how many reads to rarefy by...I added a column to our dataset of the total number of reads per sample (row)

reads<-rowSums(data.nosing[,-c(1:7)])
data.nosing.reads<-cbind(reads,data.nosing)
head(data.nosing.reads[,1:10])

hist(reads)
#we can see how many samples we have by subsetting by number of reads

dim(subset(data.nosing.reads, reads > 9999))
#113 samples of 117 remain if rarify to 1000
#107 samples of 117 remian if rarify to 10000
#111 samples remain if rarify to 5000
#Rarify to 10000, don't lose much data and more robust

# lets create rarefaction curves for each sample starting around 1000
rared<-rarefy(subset(data.nosing.reads, reads > 9999)[,-c(1:8)],sample=c(1,10,25,50,75,100,250,500,700,1250,2500,5000,10000),se=FALSE)
rared_melt<-melt(rared)
names(rared_melt)<-c("sample","sample_size","OTUs")
rared_melt$sample_size<-c(rep(1,107),rep(10,107),rep(25,107),rep(50,107),rep(75,107),rep(100,107),rep(250,107),rep(500,107),rep(700,107),rep(1250,107),rep(2500,107),rep(5000,107),rep(10000,107))
head(rared_melt)

ggplot(rared_melt)+geom_line(aes(x=sample_size,y=OTUs,colour=sample,group=sample))+theme(aspect.ratio=1)+theme_bw()

# I will decide with 10000 (107 samples)
head(data.nosing[,1:5])
data.nosing.rar<-cbind(subset(data.nosing, reads > 9999)[,1:7],rrarefy(subset(data.nosing,reads > 9999)[,-c(1:7)],10000))
head(data.nosing.rar[,1:10])

#Write to .csv file so can be pulled into subsequent analyses directly
write.csv(data.nosing.rar, file="COBS_ITS_data_rar.csv")
