#Elizabeth Bach
#28 Oct 2016
#multivariate stats for 16S, by aggregate

library(reshape)
library(vegan)

#Use "COBS_16S_rar_wsprop.csv"
data.nosing.rar<-read.csv(file.choose())
data.nosing.rar[1:10,1:10]

#remove WSprop
data.aggs<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac!="WSprop"))
str(data.aggs)

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.aggs[,1:10])
data.trans.rar<-cbind(data.aggs[,1:8],decostand(data.aggs[,-c(1:8)],"total"))

#using adonis for the test; I pulled out non-significant interactions (Date*SoilFrac*Crop, Crop*SoilFrac, Date*SoilFrac NS)
adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999)

#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rar$Date                       1    0.2168 0.21675  2.1958 0.01533 0.0246 *  
#data.trans.rar$Crop                       2    2.1000 1.05001 10.6371 0.14850 0.0001 ***
#data.trans.rar$SoilFrac                   4    0.6557 0.16393  1.6607 0.04637 0.0127 *  
#data.trans.rar$Date:data.trans.rar$Crop   2    0.3105 0.15523  1.5725 0.02195 0.0562 .  
#Residuals                               110   10.8583 0.09871         0.76785           
#Total                                   119   14.1413                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#now transforming to presence/absence (0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.aggs[,1:8])
data.trans.rar2<-cbind(data.aggs[,1:8],decostand(data.aggs[,-c(1:8)],"pa"))
head(data.trans.rar2[,1:10])

#using adonis for the test; I pulled out non-significant interactions and jaccard for presence/absence
adonis(data.trans.rar2[,-c(1:8)]~data.trans.rar2$Date*data.trans.rar2$Crop+data.trans.rar2$SoilFrac, permutations=9999)

#                                           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rar2$Date                        1    0.1727 0.17266  1.8964 0.01426 0.0168 *  
#data.trans.rar2$Crop                        2    1.1797 0.58985  6.4786 0.09743 0.0001 ***
#data.trans.rar2$SoilFrac                    4    0.5127 0.12817  1.4077 0.04234 0.0127 *  
#data.trans.rar2$Date:data.trans.rar2$Crop   2    0.2287 0.11437  1.2562 0.01889 0.0952 .  
#Residuals                                 110   10.0150 0.09105         0.82709           
#Total                                     119   12.1087                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#similar result as above, except date*crop interaction not sig

#Exploring to see which SoilFracs are different from each other
adonis(data.trans.rar[,-c(1:8)]~data.trans.rar$Date*data.trans.rar$Crop+data.trans.rar$SoilFrac, permutations=9999, data=subset(data.trans.rar, data.trans.rar$SoilFrac=="LM"|SoilFrac=="Micro"))

mds.dist<-metaMDSdist(decostand(data.trans.rar[,-c(1:8)],"pa" ),k=6,autotransform=FALSE)
SoilFrac.groups<-betadisper(mds.dist, data.trans.rar$SoilFrac, type="median")
TukeyHSD(SoilFrac.groups)
#micros and LM different

#Looking at WS measurements only (for ecosystem effects without aggregates)
data.nosing.rarWS<-droplevels(subset(data.nosing.rar, data.nosing.rar$SoilFrac=="WS"))
str(data.nosing.rarWS[,1:10])

#first transforming to scaled data within each sample (between 0 and 1) for multivariate tests.  This uses decostand() from the vegan pacakge
head(data.nosing.rarWS[,1:10])
data.trans.rarWS<-cbind(data.nosing.rarWS[,1:8],decostand(data.nosing.rarWS[,-c(1:8)],"total"))

#using adonis for the test;
adonis(data.trans.rarWS[,-c(1:8)]~data.trans.rarWS$Date*data.trans.rarWS$Crop, permutations=9999)

#no interaction, crop only main effect
#main effect only model:
adonis(data.trans.rarWS[,-c(1:8)]~data.trans.rarWS$Date+data.trans.rarWS$Crop, permutations=9999)

#                      Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#data.trans.rarWS$Date  1   0.08615 0.086153 0.84002 0.03283 0.5486   
#data.trans.rarWS$Crop  2   0.48699 0.243496 2.37415 0.18557 0.0047 **
#Residuals             20   2.05123 0.102561         0.78161          
#Total                 23   2.62437                  1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Look at crop effect within each aggregate fraction
#LM
data.nosing.rarLM<-droplevels(subset(data.aggs, data.aggs$SoilFrac=="LM"))
str(data.nosing.rarLM[,1:10])
#sacling to 0-1
data.trans.rarLM<-cbind(data.nosing.rarLM[,1:8],decostand(data.nosing.rarLM[,-c(1:8)],"total"))
str(data.trans.rarLM[,1:10])

#using adonis for the test;
adonis(data.trans.rarLM[,-c(1:8)]~data.trans.rarLM$Date*data.trans.rarLM$Crop, permutations=9999)

#
#                                            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#data.trans.rarLM$Date                        1    0.1303 0.130326 0.98354 0.04058 0.3786  
#data.trans.rarLM$Crop                        2    0.5056 0.252783 1.90769 0.15743 0.0244 *
#data.trans.rarLM$Date:data.trans.rarLM$Crop  2    0.1903 0.095165 0.71819 0.05927 0.8467  
#Residuals                                   18    2.3851 0.132508         0.74272         
#Total                                       23    3.2114                  1.00000         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Crop main effect


#MM
data.nosing.rarMM<-droplevels(subset(data.aggs, data.aggs$SoilFrac=="MM"))
str(data.nosing.rarMM[,1:10])
#sacling to 0-1
data.trans.rarMM<-cbind(data.nosing.rarMM[,1:8],decostand(data.nosing.rarMM[,-c(1:8)],"total"))
str(data.trans.rarMM[,1:10])

#using adonis for the test;
adonis(data.trans.rarMM[,-c(1:8)]~data.trans.rarMM$Date+data.trans.rarMM$Crop, permutations=9999)

#                      Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#data.trans.rarMM$Date  1   0.07617 0.076167 0.76894 0.02874 0.6010   
#data.trans.rarMM$Crop  2   0.59272 0.296362 2.99191 0.22367 0.0017 **
#Residuals             20   1.98109 0.099054         0.74759          
#Total                 23   2.64998                  1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Crop main effect


#SM
data.nosing.rarSM<-droplevels(subset(data.aggs, data.aggs$SoilFrac=="SM"))
str(data.nosing.rarSM[,1:10])
#sacling to 0-1
data.trans.rarSM<-cbind(data.nosing.rarSM[,1:8],decostand(data.nosing.rarSM[,-c(1:8)],"total"))
str(data.trans.rarSM[,1:10])

#using adonis for the test;
adonis(data.trans.rarSM[,-c(1:8)]~data.trans.rarSM$Date*data.trans.rarSM$Crop, permutations=9999)

#                                            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#data.trans.rarSM$Date                        1   0.09515 0.095147 0.96435 0.03746 0.3781   
#data.trans.rarSM$Crop                        2   0.52930 0.264649 2.68231 0.20837 0.0019 **
#data.trans.rarSM$Date:data.trans.rarSM$Crop  2   0.13975 0.069876 0.70822 0.05502 0.8751   
#Residuals                                   18   1.77597 0.098665         0.69915          
#Total                                       23   2.54017                  1.00000     

#crop main effect
#micro
data.nosing.rarMicro<-droplevels(subset(data.aggs, data.aggs$SoilFrac=="Micro"))
str(data.nosing.rarMicro[,1:10])
#sacling to 0-1
data.trans.rarMicro<-cbind(data.nosing.rarMicro[,1:8],decostand(data.nosing.rarMicro[,-c(1:8)],"total"))
str(data.trans.rarMicro[,1:10])

#using adonis for the test;
adonis(data.trans.rarMicro[,-c(1:8)]~data.trans.rarMicro$Date*data.trans.rarMicro$Crop, permutations=9999)

#                                                  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#data.trans.rarMicro$Date                           1   0.10350 0.103501 1.12431 0.04208 0.2538    
#data.trans.rarMicro$Crop                           2   0.53478 0.267389 2.90460 0.21742 0.0003 ***
#data.trans.rarMicro$Date:data.trans.rarMicro$Crop  2   0.16436 0.082182 0.89273 0.06682 0.5586    
#Residuals                                         18   1.65703 0.092057         0.67368           
#Total                                             23   2.45967                  1.00000  

#Crop main effect
