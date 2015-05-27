#Elizabeth Bach, Ryan Williams
#COBS ITS:  Pairwise comparisons of crop*date from ADONIS
#29 Jan. 2015

rm(list=ls())
library(vegan)
library(plyr)

#use "COBS_ITS_data.rar" from COBS2012_FungalSequencing\FanData
data.nosing.rar<-read.csv(file.choose())

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:10])
data.trans.rar<-cbind(data.nosing.rar[,1:8],decostand(data.nosing.rar[,-c(1:8)],"total"))

#From "multivariate_tests.R"
#using adonis for the test; I pulled out non-significant interactions
#Abundance (total reads)

adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)
#Only Date*Crop interaction sig, so removed others from model
#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rar$Date                       1     0.859 0.85860  3.0994 0.02491 0.0001 ***
#data.trans.rar$Crop                       2     4.253 2.12641  7.6760 0.12337 0.0001 ***
#data.trans.rar$SoilFrac                   4     1.488 0.37198  1.3428 0.04316 0.0048 ** 
#data.trans.rar$Date:data.trans.rar$Crop   2     1.002 0.50095  1.8083 0.02906 0.0006 ***


#Running ADONIS for each pairwise comparison
#Looking at PF only,
data.PF<-droplevels(subset(data.trans.rar, Crop=="PF"))
adonis(data.PF[,-c(1:8)]~data.PF$Date+data.PF$SoilFrac, permutations=9999)
#Date significant effect, P=0.0001

#Looking at P only,
data.P<-droplevels(subset(data.trans.rar, Crop=="P"))
adonis(data.P[,-c(1:8)]~data.P$Date+data.P$SoilFrac, permutations=9999)
#Date significant, P=0.03

#Looking at CC only,
data.CC<-droplevels(subset(data.trans.rar, Crop=="CC"))
adonis(data.CC[,-c(1:8)]~data.CC$Date+data.CC$SoilFrac, permutations=9999)
#Date sig. P=0.0001, SoilFrac sig, P=0.01

#Looking at July only,
data.July<-droplevels(subset(data.trans.rar, Date=="July2012_trimmed"))
adonis(data.July[,-c(1:8)]~data.July$Crop+data.July$SoilFrac, permutations=9999)
#Crop sig, P=0.0001; SoilFrac sig, P=0.03

data.July<-droplevels(subset(data.trans.rar, Date=="July2012_trimmed" & (Crop=="PF"|Crop=="P")))
adonis(data.July[,-c(1:8)]~data.July$Crop+data.July$SoilFrac, permutations=9999)
#P:PF diff, P=0.0001
data.July2<-droplevels(subset(data.trans.rar, Date=="July2012_trimmed" & (Crop=="CC"|Crop=="P")))
adonis(data.July2[,-c(1:8)]~data.July2$Crop+data.July2$SoilFrac, permutations=9999)
#CC:P, different, P=0.0001
data.July3<-droplevels(subset(data.trans.rar, Date=="July2012_trimmed" & (Crop=="CC"|Crop=="PF")))
adonis(data.July3[,-c(1:8)]~data.July3$Crop+data.July3$SoilFrac, permutations=9999)
#CC:PF, different, P=0.0001

#Looking at Oct only,
data.Oct<-droplevels(subset(data.trans.rar, Date=="October2012_trimmed"))
adonis(data.Oct[,-c(1:8)]~data.Oct$Crop+data.Oct$SoilFrac, permutations=9999)
#Crop sig, P=0.0001 (SoilFrac NS)

data.Oct<-droplevels(subset(data.trans.rar, Date=="October2012_trimmed"& (Crop=="PF"|Crop=="P")))
adonis(data.Oct[,-c(1:8)]~data.Oct$Crop+data.Oct$SoilFrac, permutations=9999)
#PF:P, different, P=0.0001
data.Oct2<-droplevels(subset(data.trans.rar, Date=="October2012_trimmed" & (Crop=="CC"|Crop=="P")))
adonis(data.Oct2[,-c(1:8)]~data.Oct2$Crop+data.Oct2$SoilFrac, permutations=9999)
#CC:P, different, P=0.0001
data.Oct3<-droplevels(subset(data.trans.rar, Date=="October2012_trimmed" & (Crop=="CC"|Crop=="PF")))
adonis(data.Oct3[,-c(1:8)]~data.Oct3$Crop+data.Oct3$SoilFrac, permutations=9999)
#P:PF, different, P=0.0001
