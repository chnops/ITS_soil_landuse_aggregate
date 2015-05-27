#Elizabeth Bach
#COBS ITS analysis
#Rank abundance & Relative abundance code
#7 February 2015

rm(list=ls())
library(reshape)
library(lme4)
library(lmerTest)
library(bbmle)
library(ggplot2)
library(plyr)

#Use "COBS_ITS_data_taxa2.csv" generated in "COBS_ITS_taxonomy_merge.R" file (0s NOT removed)
data_taxa2<-read.csv(file.choose())
names(data_taxa2)[11]<-paste("reads")
names(data_taxa2)[2]<-paste("OTU")
head(data_taxa2)
str(data_taxa2)

#data_melt<-melt(data_taxa2, id=c("X.2","X.1","X","OTU","CropBlock","SoilFrac","Date","Crop","block","Sample","reads"))
data_melt<-melt(data_taxa2[,-c(12,17,18)],id=c("X.2","X.1","X","OTU","CropBlock","SoilFrac","Date","Crop","block","Sample","reads","phylum","class","order","family"))
data_melt$CropDate<-factor(paste(data_melt$Crop, data_melt$Date, sep=""))
str(data_melt)
head(data_melt)

data_ordered<-arrange(data_melt, Sample, -reads)
head(data_ordered)
dim(data_ordered)

#Should be 5057 taxa per sample (5057 total unique OTUs)
#Total rows=541099/107 samples=5057
#yes!!

CC12LMJuly2012<-droplevels(subset(data_ordered, data_ordered$Sample=="CC12-LM-July2012_trimmed"))
length(CC12LMJuly2012$OTU)
#Confirm yes!!

data_ordered$rank<-rep(seq(1,5057,1),107)
dim(data_ordered)
head(data_ordered)
data_ordered[5055:5060,]
range(data_ordered$rank)

data_ordered2<-data.frame(cast(data_ordered, CropBlock + SoilFrac + Crop + OTU +class ~ Date, fill=0))
head(data_ordered2)
str(data_ordered2)
dim(data_ordered2)
range(data_ordered2$July2012_trimmed)
count(data_ordered2$July2012_trimmed==1)
#data is unbalanced due to dropped samples, e.g. there are Oct samples for which we do not include a corresponding July sample
#using 0's for "missing values"

DateCrop.rank<-ddply(data_ordered2, .(OTU), mutate, delta=July2012_trimmed-October2012_trimmed)
head(DateCrop.rank)
dim(DateCrop.rank)
str(DateCrop.rank)

#who are top-ranked taxa for each CropDate, SoilFrac?
#Be sure the carry the "species" colmun with this!!

#Which taxa have highest abs(delta) between sampling events (by crop)?
range(DateCrop.rank$July2012_trimmed)
Corn<-droplevels(subset(DateCrop.rank, DateCrop.rank$Crop=="CC"))
Corn2<-subset(Corn, Corn$July2012_trimmed>0)
range(Corn2$July2012_trimmed)
dim(Corn2)
head(Corn2)
str(Corn2)
head(arrange(Corn2, July2012_trimmed))

Corn.class<-ddply(Corn2, .(class), summarize, mean.July=mean(July2012_trimmed), se.July=sd(July2012_trimmed)/sqrt(length(July2012_trimmed)-1), sum.July=sum(July2012_trimmed),
 mean.Oct=mean(October2012_trimmed), se.Oct=sd(October2012_trimmed)/sqrt(length(October2012_trimmed)-1), sum.Oct=sum(October2012_trimmed))
head(Corn.class)
arrange(Corn.class, mean.July)[1:10,]
arrange(Corn.class, sum.July)[1:10,]
arrange(Corn.class, mean.Oct)[1:10,]
arrange(Corn.class, sum.Oct)[1:10,]

arrange(Corn.July, sum)[1:10,]
#Same top taxa regardless of mean or sum

Corn.July<-ddply(Corn, .(OTU), summarize, mean=mean(July2012_trimmed))
head(Corn.July)
head(arrange(Corn.July, mean))

#PF
PF<-droplevels(subset(DateCrop.rank, DateCrop.rank$Crop=="PF"))
PF2<-subset(PF, PF$July2012_trimmed>0)
range(PF2$July2012_trimmed)
dim(PF2)
head(PF2)
str(PF2)
head(arrange(PF2, July2012_trimmed))

PF.class<-ddply(PF2, .(class), summarize, mean.July=mean(July2012_trimmed), se.July=sd(July2012_trimmed)/sqrt(length(July2012_trimmed)-1), sum.July=sum(July2012_trimmed),
 mean.Oct=mean(October2012_trimmed), se.Oct=sd(October2012_trimmed)/sqrt(length(October2012_trimmed)-1), sum.Oct=sum(October2012_trimmed))
head(PF.class)
arrange(PF.class, mean.July)[1:10,]
arrange(PF.class, sum.July)[1:10,]
arrange(PF.class, mean.Oct)[1:10,]
arrange(PF.class, sum.Oct)[1:10,]

#P
P<-droplevels(subset(DateCrop.rank, DateCrop.rank$Crop=="P"))
P2<-subset(P, P$July2012_trimmed>0)
range(P2$July2012_trimmed)
dim(P2)
head(P2)
str(P2)
head(arrange(P2, July2012_trimmed))


P.class<-ddply(P2, .(class), summarize, mean.July=mean(July2012_trimmed), se.July=sd(July2012_trimmed)/sqrt(length(July2012_trimmed)-1), sum.July=sum(July2012_trimmed),
 mean.Oct=mean(October2012_trimmed), se.Oct=sd(October2012_trimmed)/sqrt(length(October2012_trimmed)-1), sum.Oct=sum(October2012_trimmed))
head(P.class)
arrange(P.class, mean.July)[1:10,]
arrange(P.class, sum.July)[1:10,]
arrange(P.class, mean.Oct)[1:10,]
arrange(P.class, sum.Oct)[1:10,]


#Relative abundance
# #reads/total reads per sample
data.relative<-ddply(data_ordered, .(Sample,OTU, Crop, SoilFrac, Date, phylum,class,order, family), summarize, .progress="text",
relative_read=(reads/10000))

sum_relative<-ddply(data.relative, .(Sample), summarize, .progress="text", sum=sum(relative_read))
head(sum_relative)
#check, sum of all proportions in each sample is 1

head(data.relative)
max(data.relative$relative_read)
dim(data.relative)

#remove 0s
data.relative2<-subset(data.relative, data.relative$relative_read>0)
dim(data.relative2)
head(data.relative2)

data.relative_order<-arrange(data.relative2, Sample, -relative_read)
head(data.relative_order)

#phyla abundance
phylum.abund<-ddply(data.relative2, .(phylum, Crop), summarize, .progress="text",
mean=mean(relative_read*100), n=length(relative_read), SE=sd(relative_read*100)/(sqrt(n-1)))
head(phylum.abund)
phylum.abund_ordered<-arrange(phylum.abund, Crop, -mean)
phylum.abund_ordered

#class abundance
class.abund<-ddply(data.relative2, .(class, Crop, Date), summarize, .progress="text",
mean=mean(relative_read), n=length(relative_read), SE=sd(relative_read*100)/(sqrt(n-1)))
class.abund_ordered<-arrange(class.abund, Crop,Date, -mean)
class.abund_ordered

#order abundance
order.abund<-ddply(data.relative2, .(order, Crop), summarize, .progress="text",
mean=mean(relative_read))
order.abund_ordered<-arrange(order.abund, Crop, -mean)
order.abund_ordered

#family abundance
family.abund<-ddply(data.relative2, .(family, Crop), summarize, .progress="text",
mean=mean(relative_read))
family.abund_ordered<-arrange(family.abund, Crop, -mean)
family.abund_ordered

#SoilFrac
#class abundance
class.abund<-ddply(data.relative2, .(class, SoilFrac), summarize, .progress="text",
mean=mean(relative_read))
class.abund_ordered<-arrange(class.abund, SoilFrac, -mean)
class.abund_ordered

#order abundance
order.abund<-ddply(data.relative2, .(order, SoilFrac), summarize, .progress="text",
mean=mean(relative_read))
order.abund_ordered<-arrange(order.abund, SoilFrac, -mean)
order.abund_ordered

#family abundance
family.abund<-ddply(data.relative2, .(family, SoilFrac), summarize, .progress="text",
mean=mean(relative_read))
family.abund_ordered<-arrange(family.abund, SoilFrac, -mean)
family.abund_ordered

#Group OTUs by class first, then relativize
#Sum reads in sample by Class
class.sum<-ddply(data_ordered, .(class, Sample, Crop, Date, SoilFrac), summarize, .progress="text",
sum=sum(reads))
head(class.sum)

class.relative<-ddply(class.sum, .(Sample,Crop, Date, SoilFrac, class), summarize, .progress="text",
relative_read=(sum/10000))
head(class.relative)

sum_relative<-ddply(class.relative, .(Sample), summarize, .progress="text", sum=sum(relative_read))
head(sum_relative)
#check, sum of all proportions in each sample is 1

class.abund2<-ddply(class.relative, .(class, Crop, Date), summarize, .progress="text",
mean=mean(relative_read), n=length(relative_read), SE=sd(relative_read)/(sqrt(n-1)))
class.abund2_ordered<-arrange(class.abund2, Crop,Date, -mean)
class.abund2_ordered

#SoilFrac
#class abundance
class.SoilFrac<-ddply(class.relative, .(class, SoilFrac), summarize, .progress="text",
mean=mean(relative_read),n=length(relative_read), SE=sd(relative_read)/(sqrt(n-1)))
class.SoilFrac_ordered<-arrange(class.SoilFrac, SoilFrac, -mean)
class.SoilFrac_ordered

#Sum reads in sample by Plylum
phylum.sum<-ddply(data_ordered, .(phylum, Sample, Crop, Date, SoilFrac, block), summarize, .progress="text",
sum=sum(reads))
head(phylum.sum)

phylum.relative<-ddply(phylum.sum, .(Sample,Crop, Date, SoilFrac, phylum, block), summarize, .progress="text",
relative_read=(sum/10000))
head(phylum.relative)

sum_relative<-ddply(phylum.relative, .(Sample), summarize, .progress="text", sum=sum(relative_read))
head(sum_relative)
#check, sum of all proportions in each sample is 1

#Cropping system means
phylum.abund2<-ddply(phylum.relative, .(phylum, Crop), summarize, .progress="text",
mean=mean(relative_read), n=length(relative_read), SE=sd(relative_read)/(sqrt(n-1)))
phylum.abund2_ordered<-arrange(phylum.abund2, Crop, -mean)
phylum.abund2_ordered

#Model comaprisons
model.null<-lmer(relative_read~1+(1|block), data=phylum.relative, REML=FALSE)
model.full<-lmer(relative_read~Date*Crop*SoilFrac+(1|block), data=phylum.relative, REML=FALSE)
model.main<-lmer(relative_read~Date+Crop+SoilFrac+(1|block), data=phylum.relative, REML=FALSE)
AICtab(model.null,model.full,model.main)
anova(model.null,model.full,model.main)
anova(model.main)

