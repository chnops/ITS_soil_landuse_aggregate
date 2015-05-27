#Elizabeth Bach, Ryan Williams
#5 Jan 2015
#multivariate stats to analyze COBS ITS communities

rm(list=ls())
library(reshape)
library(vegan)

#Use "COBS_ITS_data_rar.csv", will perform analysis using OTUs only (no taxonomic summaries)
data.nosing.rar<-read.csv(file.choose())

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:10])
data.trans.rar<-cbind(data.nosing.rar[,1:8],decostand(data.nosing.rar[,-c(1:8)],"total"))

#using adonis for the test; I pulled out non-significant interactions (Date*SoilFrac*Crop, Crop*SoilFrac, Date*SoilFrac NS)
adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)

#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rar$Date                       1     0.859 0.85860  3.0994 0.02491  1e-04 ***
#data.trans.rar$Crop                       2     4.253 2.12641  7.6760 0.12337  1e-04 ***
#data.trans.rar$SoilFrac                   4     1.488 0.37198  1.3428 0.04316  5e-03 ** 
#data.trans.rar$Date:data.trans.rar$Crop   2     1.002 0.50095  1.8083 0.02906  2e-04 ***
#Residuals                                97    26.871 0.27702         0.77950           
#Total                                   106    34.472                 1.00000           

#now transforming to presence/absence (0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rar[,1:8])
data.trans.rar<-cbind(data.nosing.rar[,1:8],decostand(data.nosing.rar[,-c(1:8)],"pa"))

#using adonis for the test; I pulled out non-significant interactions and jaccard for presence/absence
adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)
#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rar$Date                       1    0.6790 0.67900  2.8405 0.02326  1e-04 ***
#data.trans.rar$Crop                       2    3.1380 1.56899  6.5637 0.10748  1e-04 ***
#data.trans.rar$SoilFrac                   4    1.3642 0.34106  1.4268 0.04673  6e-04 ***
#data.trans.rar$Date:data.trans.rar$Crop   2    0.8286 0.41430  1.7332 0.02838  6e-04 ***
#Residuals                                97   23.1871 0.23904         0.79416           
#Total                                   106   29.1969                 1.00000          

#Exploring to see which SoilFracs are different from each other
adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999, data=subset(data.trans.rar, data.trans.rar$SoilFrac=="LM"|SoilFrac=="Micro"))

mds.dist<-metaMDSdist(decostand(data.trans.rar[,-c(1:8)],"pa" ),k=6,autotransform=FALSE)
SoilFrac.groups<-betadisper(mds.dist, data.trans.rar$SoilFrac, type="median")
TukeyHSD(SoilFrac.groups)
#micros and SM different from LM

#Looking at management systme differences within each aggregate fraction
#Looking at WS measurements only (for ecosystem effects without aggregates)
data.nosing.rarWS<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="WS"))
str(data.nosing.rarWS[,1:10])

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rarWS[,1:10])
data.trans.rarWS<-cbind(data.nosing.rarWS[,1:8],decostand(data.nosing.rarWS[,-c(1:8)],"total"))

#using adonis for the test;
adonis(data.trans.rarWS[,-c(1:8)]~data.trans.rarWS$Date*data.trans.rarWS$Crop, permutations=9999)
#no Date*Crop interaction, look at fit of main effects model only
adonis(data.trans.rarWS[,-c(1:8)]~data.trans.rarWS$Date+data.trans.rarWS$Crop, permutations=9999)
#Essentially same answer and P-vals:
#                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarWS$Date  1    0.3048 0.30480  1.0468 0.04447 0.3308    
#data.trans.rarWS$Crop  2    1.3086 0.65430  2.2472 0.19091 0.0001 ***
#Residuals             18    5.2410 0.29117         0.76462           
#Total                 21    6.8545                 1.00000         

mds.dist<-metaMDSdist(decostand(data.trans.rarWS[,-c(1:8)],"range" ),k=6,autotransform=FALSE)
Crop.groups<-betadisper(mds.dist, data.trans.rarWS$Crop, type="median")
TukeyHSD(Crop.groups)

#now transforming to presence/absence (0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rarWS[,1:8])
data.trans.rarWSpa<-cbind(data.nosing.rarWS[,1:8],decostand(data.nosing.rarWS[,-c(1:8)],"pa"))

#using adonis for the test; I pulled out non-significant interactions and jaccard for presence/absence
adonis(data.trans.rarWSpa[,-c(1:8)]~data.trans.rarWSpa$Date*data.trans.rarWSpa$Crop, permutations=9999)
#no interaction of crop and date
adonis(data.trans.rarWSpa[,-c(1:8)]~data.trans.rarWSpa$Date+data.trans.rarWSpa$Crop, permutations=9999)
#essentially same answer as full model, also same out come as scaled data
#                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarWSpa$Date  1    0.3055 0.30552  1.2864 0.05384 0.0951 .  
#data.trans.rarWSpa$Crop  2    1.0936 0.54681  2.3023 0.19273 0.0001 ***
#Residuals               18    4.2751 0.23751         0.75342           
#Total                   21    5.6743                 1.00000 

mds.distPA<-metaMDSdist(decostand(data.trans.rarWSpa[,-c(1:8)],"pa" ),k=6,autotransform=FALSE)
Crop.groupsPA<-betadisper(mds.distPA, data.trans.rarWSpa$Crop, type="median")
TukeyHSD(Crop.groupsPA)


#LM, large macroagregates, >2000 um
data.nosing.rarLM<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="LM"))
str(data.nosing.rarLM[,1:10])
#sacling to 0-1
data.trans.rarLM<-cbind(data.nosing.rarLM[,1:8],decostand(data.nosing.rarLM[,-c(1:8)],"total"))
str(data.trans.rarLM[,1:10])

#using adonis for the test;
adonis(data.trans.rarLM[,-c(1:8)]~data.trans.rarLM$Date*data.trans.rarLM$Crop, permutations=9999)

#Crop significant, Date not, interaction not
#
#                                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarLM$Date                        1    0.3618 0.36180 1.02978 0.04446 0.3730    
#data.trans.rarLM$Crop                        2    1.1626 0.58128 1.65449 0.14287 0.0001 ***
#data.trans.rarLM$Date:data.trans.rarLM$Crop  2    0.6400 0.31998 0.91076 0.07865 0.7964    
#Residuals                                   17    5.9727 0.35134         0.73402           
#Total                                       22    8.1371                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#MM, medium macroaggregates, 1000-2000 um
data.nosing.rarMM<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="MM"))
str(data.nosing.rarMM[,1:10])
#sacling to 0-1
data.trans.rarMM<-cbind(data.nosing.rarMM[,1:8],decostand(data.nosing.rarMM[,-c(1:8)],"total"))
str(data.trans.rarMM[,1:10])

#using adonis for the test;
adonis(data.trans.rarMM[,-c(1:8)]~data.trans.rarMM$Date+data.trans.rarMM$Crop, permutations=9999)

#Crop significant, date marginal, no interaction, same result with main effects model only
#                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarMM$Date  1    0.3975 0.39751  1.3692 0.05950 0.0623 .  
#data.trans.rarMM$Crop  2    1.3484 0.67421  2.3223 0.20182 0.0001 ***
#Residuals             17    4.9354 0.29032         0.73869           
#Total                 20    6.6814                 1.00000           

#SM, small macroaggregates, 250-1000 um
data.nosing.rarSM<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="SM"))
str(data.nosing.rarSM[,1:10])
#sacling to 0-1
data.trans.rarSM<-cbind(data.nosing.rarSM[,1:8],decostand(data.nosing.rarSM[,-c(1:8)],"total"))
str(data.trans.rarSM[,1:10])

#using adonis for the test;
adonis(data.trans.rarSM[,-c(1:8)]~data.trans.rarSM$Date*data.trans.rarSM$Crop, permutations=9999)

#Crop effect, no date, no interaction, main effects model shows same
#                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarSM$Date  1    0.2763 0.27632  1.0681 0.04648 0.2992    
#data.trans.rarSM$Crop  2    1.2703 0.63517  2.4552 0.21369 0.0001 ***
#Residuals             17    4.3981 0.25871         0.73983           
#Total                 20    5.9447                 1.00000   

#micro, microaggregates, <250 um
data.nosing.rarMicro<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="Micro"))
str(data.nosing.rarMicro[,1:10])
#sacling to 0-1
data.trans.rarMicro<-cbind(data.nosing.rarMicro[,1:8],decostand(data.nosing.rarMicro[,-c(1:8)],"total"))
str(data.trans.rarMicro[,1:10])

#using adonis for the test;
adonis(data.trans.rarMicro[,-c(1:8)]~data.trans.rarMicro$Date*data.trans.rarMicro$Crop, permutations=9999)

#crop effect and date effect, no interaction
#                                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarMicro$Date                           1    0.5012 0.50117  2.2875 0.09418 0.0006 ***
#data.trans.rarMicro$Crop                           2    1.2251 0.61256  2.7959 0.23023 0.0001 ***
#data.trans.rarMicro$Date:data.trans.rarMicro$Crop  2    0.5278 0.26390  1.2045 0.09918 0.1287    
#Residuals                                         14    3.0673 0.21909         0.57641           
#Total                                             19    5.3214                 1.00000 

#see "ITS_ADONIS_pairwise.R" for pairwise comparison to untangle Date*Crop interation