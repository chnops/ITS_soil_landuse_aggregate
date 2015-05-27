#Elizabeth Bach, Ryan Williams
#Merging OTU table with taxonomy
#25 May 2015

#removes all objects previously assigned in R session, ensures no conflicting objects from another analysis
rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

#use "COBS_ITS_data_rar.csv", generated in "ITS_rarifying.R"
data.nosing.rar<-read.csv(file.choose())
str(data.nosing.rar)
dim(data.nosing.rar)
head(data.nosing.rar[,1:10])

#melt data to integrat the OTU information within rows, this creats a much longer data fram as each OTU occurance in each sample now is a row, not a column
data_melt<-melt(data.nosing.rar, id=c("X.1","X","CropBlock","SoilFrac","Date","Crop","block","Sample"))
head(data_melt)
dim(data_melt)
#We will merge the taxonomic information with the melted data frame

#Use "COBS_ITS_taxonomy.txt" file to import taxonomic info, this will have OTUs as rows, taxonomic info in columns
taxonomy<-read.delim(file.choose())
head(taxonomy)
head(data_melt)
#OTU column titled "X" in the taxonomy file, title "variable" in the melted data frame
#we merge the data frames together by integrating these columns
data_taxa<-merge(data_melt,taxonomy,by.x="variable",by.y="X")
head(data_taxa)
dim(data_taxa)
#create .csv file with taxonomy that can be used in subsequent analyses
write.csv(data_taxa, file="COBS_ITS_data_taxa2.csv")

#Double checkng data matches the original COBS_ITS_data_rar.csv"
str(data_taxa2)
#5059 OTUs (matches rarification step!)
#107 samples, also matches rarifaction step

#Melted data format allows us to easily find the total number of OTUs per sample
#OTUs per sample
#remove 0s (occurances where OTUs do not occur in an indivdual sample), so only OTUs present in sample are in data frame
data_taxa2<-data.frame(data_taxa[ which(data_taxa$value>0),])
sample.richness<-ddply(data_taxa2, .(Sample), summarise, .progress="text", OTUr=length(variable))
dim(sample.richness)
head(sample.richness)
min(sample.richness$OTUr)
max(sample.richness$OTUr)
