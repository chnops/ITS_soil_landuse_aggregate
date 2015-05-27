#Elizabeth Bach
#COBS ITS
#Figure 1 for manuscript
#Fungal community NMDS

rm(list=ls())
library(reshape)
library(grid)
library(ggplot2)
library(vegan)
library(gridExtra)

#We will use taxonomy to look at phylum-level centroids within the NMDS
#use "COBS_ITS_data_taxa.csv" generated in "COBS_ITS_taxonomy_merge.R" code
data_taxa<-read.csv(file.choose())
names(data_taxa[,1:16])
data_phyla<-data.frame(cast(data_taxa, Sample~phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_phyla)

#We will use OTU data to generate the NMDS scores, merge with Phylum-summarized data
#Use "COBS_ITS_data_rar.csv"
data.nosing.rar<-read.csv(file.choose())
head(data.nosing.rar[,1:10])
merged_taxa<-merge(data_phyla, data.nosing.rar, by="Sample")
dim(merged_taxa)
names(merged_taxa[,1:30])
rar.levels<-levels(merged_taxa$Sample)

#Environmental metadata, Use data.metadata2.csv
data.metadata2<-read.csv(file.choose())
head(data.metadata2)
str(data.metadata2)
#ensure environmental data reported at same levels (aggregate fractions) as OTU info
data.metadata3<-data.metadata2[data.metadata2$Sample %in% c(rar.levels),]
dim(data.metadata3)
head(data.metadata3)

#Presence/Absence MDS
mds.pa<-metaMDS(decostand(merged_taxa[,-c(1:14)],"pa" ),k=6,autotransform=FALSE, na.rm=TRUE)

#Looking at taxonomic correlations
IntVectors1<-envfit(mds.pa, data_phyla[,2:6], na.rm=TRUE)
IntVectors1
vectors<-data.frame(IntVectors1$vectors[1:4])
vectors
names<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Unk")
IntVectors2<-data.frame(names, vectors)
IntVectors3<-(subset(IntVectors2, pvals<0.053))
IntVectors3

#Environmental vectors
envectors1<-envfit(mds.pa, data.metadata3[,7:24], na.rm=TRUE)
head(envectors1)
vectors2<-data.frame(envectors1$vectors[1:4])
vectors2
names<-rownames(vectors2)
vectors3<-subset(data.frame(names,vectors2), pvals<0.051)
vectors3
vectors4<-vectors3[1:3,]
vectors4

#Abundance
mds.ab<-metaMDS(decostand(merged_taxa[,-c(1:14)],"total" ),k=6,autotransform=FALSE, na.rm=TRUE)

IntVectors1ab<-envfit(mds.ab, data_phyla[,2:6], na.rm=TRUE)
IntVectors1ab
vectors_ab<-data.frame(IntVectors1ab$vectors[1:4])
vectors_ab
names<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Unk")
IntVectors2ab<-subset(data.frame(names, vectors_ab), pvals<0.05)
IntVectors2ab

#Environmental vectors
envectors1ab<-envfit(mds.ab, data.metadata3[,-c(1:6)], na.rm=TRUE)
head(envectors1ab)
vectors2ab<-data.frame(envectors1ab$vectors[1:4])
vectors2ab
names<-rownames(vectors2ab)
vectors3ab<-subset(data.frame(names,vectors2ab), pvals<0.051)
vectors3ab
vectors4ab<-vectors3ab[vectors3ab$names %in% c("ph","MBN","ExtN","BD","TP","RootBiomass","MWD_um"),]
vectors4ab
vectors5ab<-vectors4ab[1:6,]
vectors5ab

#For CropDate Interaction figures
#NMDS plotting function, from R. Williams
ggplot.NMDS<-function(XX,ZZ,COLORS){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=2.5,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1.5, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

merged_taxa$CropDate<-as.factor(paste(merged_taxa$Crop, merged_taxa$Date))
levels(merged_taxa$CropDate)
#color palate from iwanthue
color.1<-rgb(97,146,131, max=255)
color.2<-rgb(132,221,200, max=255)
color.3<-rgb(223,157,57, max=255)
color.4<-rgb(178,130,73, max=255)
color.5<-rgb(218,70,65, max=255)
color.6<-rgb(194,114,99, max=255)
colors.cropdate<-c(color.5,color.6,color.1,color.2,color.3,color.4)

#Abundance figure
env.vjust<-c("1","1","-0.25","1","0.5","-0.25","0")
env.hjust<-c("0.5","0.5","0","0.5","1","0.75","-0.1")
Abund.CropDate<-ggplot.NMDS(mds.ab, (merged_taxa$CropDate), colors.cropdate)+
geom_point(data=IntVectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=1.25,size=4,fontface="bold")+
geom_text(data=vectors4ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust,vjust=env.vjust,size=4,fontface="bold")+geom_segment(data=vectors4ab, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position=c(0.95, 0.15), legend.background=element_blank(), legend.text=element_text(size=12, face="bold"),legend.key=element_blank(),legend.title=element_blank(),panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
Abund.CropDate

#Presence/Absence figure
env.vjust2<-c("-0.5","-0.5","1")
env.hjust2<-c("0.5","0.5","0.5")
PA.CropDate<-ggplot.NMDS(mds.pa, (merged_taxa$CropDate), colors.cropdate)+
geom_point(data=IntVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=1.25,size=4,fontface="bold")+
geom_text(data=vectors4,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust2,vjust=env.vjust2,size=4,fontface="bold")+geom_segment(data=vectors4, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none",panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
PA.CropDate

#SoilFrac figures
#Aggregate colors (from iWantHue)
color.1<-rgb(22,95,103, max=255)
color.2<-rgb(160,11,115, max=255)
color.3<-rgb(155,238,179, max=255)
color.4<-rgb(56,20,15, max=255)
color.5<-rgb(18,101,47, max=255)
colors.agg<-c(color.1,color.2,color.3,color.4,color.5)

#NMDS plotting function altered to only include LM and Micro centroids, for clarity
ggplot.NMDS2<-function(XX,ZZ,COLORS){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ella <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ella <- rbind(df_ella, cbind(as.data.frame(with(NMDS[NMDS$Treatment=="LM",],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group="LM"))
  }

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ellb <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ellb <- rbind(df_ellb, cbind(as.data.frame(with(NMDS[NMDS$Treatment=="Micro",],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group="Micro"))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=2.5,alpha=0.75) +
    geom_path(data=df_ella, aes(x=MDS1, y=MDS2,colour=group), size=1.5, linetype=5)+geom_path(data=df_ellb, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

#Presence Absence, hjust, vjust positions above work fine
PA.Agg<-ggplot.NMDS2(mds.pa, (merged_taxa$SoilFrac), colors.agg)+geom_point(data=IntVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",inherit_aes=FALSE, size=3)+
geom_text(data=IntVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4, fontface="bold", vjust=1.25)+
geom_text(data=vectors4,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust2,vjust=env.vjust2,size=4,fontface="bold")+geom_segment(data=vectors4, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none", panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
PA.Agg

#code for legend when needed
#legend.position=c(0.9, 0.1), legend.background=element_blank(), legend.text=element_text(size=12, face="bold"),legend.key=element_blank(),legend.title=element_blank(),
#legend.position="none"

#Abundance
env.vjust2<-c("1","1","-0.4","1.25","-0.2","0.1")
env.hjust2<-c("0.5","0.5","1","0.75","0","1")

Abund.Agg<-ggplot.NMDS2(mds.ab, (merged_taxa$SoilFrac), colors.agg)+
geom_point(data=IntVectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=1.25,size=4,fontface="bold")+
geom_text(data=vectors5ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust2,vjust=env.vjust2,size=4,fontface="bold")+geom_segment(data=vectors5ab, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none", legend.background=element_blank(), legend.text=element_text(size=12, face="bold"),legend.key=element_blank(),legend.title=element_blank(),panel.background=element_blank(), axis.text=element_text(size=1, face="bold", colour="black"), axis.title=element_text(size=1.5, face="bold", colour="black"))
Abund.Agg

#Manuscript Figure, taxa+vectors
#Presence/Absence, CropDate effect
env.vjust2<-c("-0.75","-0.5","1")
env.hjust2<-c("0.5","0.75","0.75")
taxa.vjust<-c("0.1","-0.5","0")
taxa.hjust<-c("-0.1","0.5","-0.25")
PA.CropDate<-ggplot.NMDS(mds.pa, (merged_taxa$CropDate), colors.cropdate)+
geom_point(data=IntVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=taxa.vjust, hjust=taxa.hjust,size=4,fontface="bold")+
geom_text(data=vectors4,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust2,vjust=env.vjust2,size=4,fontface="bold")+geom_segment(data=vectors4, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none",panel.background=element_blank(), axis.text=element_text(size=6, face="bold", colour="black"), axis.title=element_text(size=8, face="bold", colour="black"))

#Abundance, CropDate effect
env.vjust<-c("1","1","-0.25","1","0.5","-0.2")
env.hjust<-c("0.5","0.5","0","0.5","1","0.75")
taxa.vjust2<-c("-0.5","1.25","0")
taxa.hjust2<-c("0.5","0.5","-0.25")
Abund.CropDate<-ggplot.NMDS(mds.ab, (merged_taxa$CropDate), colors.cropdate)+
geom_point(data=IntVectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=taxa.vjust2,hjust=taxa.hjust2,size=4,fontface="bold")+
geom_text(data=vectors5ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust,vjust=env.vjust,size=4,fontface="bold")+geom_segment(data=vectors5ab, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none",panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))

#Presence/Absence, aggregate effect
env.vjust2<-c("-0.75","-0.5","1")
env.hjust2<-c("0.5","0.75","0.75")
taxa.vjust<-c("0.1","-0.5","0")
taxa.hjust<-c("-0.1","0.5","-0.25")
PA.Agg<-ggplot.NMDS2(mds.pa, (merged_taxa$SoilFrac), colors.agg)+geom_point(data=IntVectors3, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",inherit_aes=FALSE, size=3)+
geom_text(data=IntVectors3,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),size=4, fontface="bold", vjust=taxa.vjust, hjust=taxa.hjust)+
geom_text(data=vectors4,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust2,vjust=env.vjust2,size=4,fontface="bold")+geom_segment(data=vectors4, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none", panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))

#Abundance, aggregate effect
env.vjust<-c("1","1","-0.25","1","0.5","-0.2")
env.hjust<-c("0.5","0.5","0","0.5","1","0.75")
taxa.vjust2<-c("-0.75","1.25","0")
taxa.hjust2<-c("0.5","0.5","-0.25")
Abund.Agg<-ggplot.NMDS2(mds.ab, (merged_taxa$SoilFrac), colors.agg)+
geom_point(data=IntVectors2ab, aes(x=arrows.NMDS1,y=arrows.NMDS2),colour="black",size=3,inherit_aes=FALSE)+geom_text(data=IntVectors2ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),vjust=taxa.vjust2, hjust=taxa.hjust2,size=4,fontface="bold")+
geom_text(data=vectors5ab,aes(x=arrows.NMDS1,y=arrows.NMDS2,label=names),hjust=env.hjust,vjust=env.vjust,size=4,fontface="bold")+geom_segment(data=vectors5ab, aes(x=0,xend=arrows.NMDS1,y=0,yend=arrows.NMDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none", legend.background=element_blank(), legend.text=element_text(size=12, face="bold"),legend.key=element_blank(),legend.title=element_blank(),panel.background=element_blank(), axis.text=element_text(size=1, face="bold", colour="black"), axis.title=element_text(size=1.5, face="bold", colour="black"))

grid.arrange(PA.CropDate, Abund.CropDate,PA.Agg, Abund.Agg, ncol=2)

