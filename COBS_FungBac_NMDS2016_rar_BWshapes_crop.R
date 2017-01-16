#Elizabeth Bach
#COBS: 16S+ITS joint figures
#28 April 2016

rm(list=ls())
library(reshape)
library(grid)
library(ggplot2)
library(vegan)
library(gridExtra)

#use "COBS_ITS_data_taxa.csv" generated in "COBS_ITS_taxonomy_merge.R" code
data_taxa<-read.csv(file.choose())
names(data_taxa[,1:16])
data_phyla<-data.frame(cast(data_taxa, Sample~phylum, value="value", fun.aggregate=sum, add.missing=TRUE))
head(data_phyla)
#Use "COBS_ITS_data_rar.csv"
data.nosing<-read.csv(file.choose())
head(data.nosing[,1:10])
merged_taxa<-merge(data_phyla, data.nosing, by="Sample")
dim(merged_taxa)
names(merged_taxa[,1:30])
#remove whole soil for clarity
rar.levels<-levels(merged_taxa$Sample)
str(merged_taxa)
merged_taxa2<-droplevels(subset(merged_taxa,merged_taxa$SoilFrac!="WS"))
dim(merged_taxa2)
str(merged_taxa2)
samples.ITS<-(merged_taxa2$Sample)

levels(merged_taxa2$SoilFrac)
merged_taxa2$SoilFrac<-factor(merged_taxa2$SoilFrac, levels=c("Micro","SM","MM","LM"))
levels(merged_taxa2$SoilFrac)

#Environmental metadata, Use data.metadata2.csv
data.metadata2<-read.csv(file.choose(),na.strings=".")
head(data.metadata2)
str(data.metadata2)
all.levels<-levels(data.metadata2$Sample)
diff(all.levels,rar.levels)
drop.rar<-setdiff(all.levels,rar.levels)
data.metadata3<-data.metadata2[data.metadata2$Sample %in% c(rar.levels),]
dim(data.metadata3)
head(data.metadata3)

#remove whole soil for clarity

#Presence/Absence MDS
mds.pa2<-metaMDS(decostand(merged_taxa2[,-c(1:12)],"pa" ),k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#stress=0.25, on the cusp, 3-D probably better, but not sure if we want to go down that hole

#Looking at taxonomic correlations
IntVectors1<-envfit(mds.pa2, data_phyla[,2:7], na.rm=TRUE)
IntVectors1
vectors<-data.frame(IntVectors1$vectors[1:4])
vectors
names<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Unk","Zygomycota")
IntVectors2<-data.frame(names, vectors)
IntVectors3<-(subset(IntVectors2, pvals<0.053))
IntVectors3
#Environmental vectors
envectors1<-envfit(mds.pa2, data.metadata3[,7:24], na.rm=TRUE)
head(envectors1)
vectors2<-data.frame(envectors1$vectors[1:4])
vectors2
names<-rownames(vectors2)
vectors3<-subset(data.frame(names,vectors2), pvals<0.051)
vectors3
vectors4<-vectors3[1:3,]
vectors4

#SoilFrac
#Aggregate colors (from iWantHue)
color.1<-rgb(144,100,34, max=255)
color.2<-rgb(146,118,187, max=255)
color.3<-rgb(144,193,63, max=255)
color.4<-rgb(172,74,97, max=255)
colors.agg<-c(color.1,color.2,color.3,color.4)

#Aggregate shapes
shapes<-c(15,16,17,18)

#Agg colors B/W
colors.bw<-c("black","gray60","grey20","gray40")

#Crop shapes
shapes.crop<-c(15,16,17)
#Crop colors B/W
colors.crop.bw<-c("black","gray50","grey20")

ggplot.NMDS2<-function(XX,ZZ,COLORS,SHAPES){
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

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_path(data=df_ella, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE, size=1.5, linetype=5)+geom_path(data=df_ellb, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE,size=1.5, linetype=5)+geom_point(aes(shape=Treatment, colour=Treatment), size=3) +
    theme_bw()+theme(aspect.ratio=1)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+scale_color_manual(values=COLORS)+scale_shape_manual(values=SHAPES)+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

#Presence Absence

#Fungi, base figure, no metadata
PA.Agg1<-ggplot.NMDS2(mds.pa2, (merged_taxa2$SoilFrac), colors.bw, shapes)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), legend.position="none",panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))+
annotate("text", x=0.6,y=1.15, label="Fungi", size=12)
PA.Agg1




#16S data
#Use "data_rar_WSprop.csv"

data.16S<-read.csv(file.choose(), header=TRUE)
head(data.16S)
dim(data.16S)
str(data.16S)

#remove WS, WSprop
levels(data.16S$SoilFrac)
data.16S2<-droplevels(subset(data.16S,data.16S$SoilFrac!="WS"))
data.16S3<-droplevels(subset(data.16S2,data.16S2$SoilFrac!="WSprop"))
dim(data.16S3)
str(data.16S3)
levels(data.16S3$SoilFrac)

data.16S3$SoilFrac<-factor(data.16S3$SoilFrac, levels=c("LM","MM","SM","Micro"), order=TRUE)
levels(data.16S3$SoilFrac)
str(data.16S3[,1:10])

Treatment<-data.16S3$SoilFrac
levels(Treatment)
MDS1<-data.frame(scores(mds.16S.pa2))$NMDS1
MDS2<-data.frame(scores(mds.16S.pa2))$NMDS2

NMDS<-data.frame(MDS1,MDS2,Treatment)
str(NMDS)

#what samples are in 16S that are not in ITS?)
setdiff(samples.ITS, samples.16S)
#OK maybe not worth worry about?  On the other hand, it will be obvious, so I guess explicitly state in the methods?

#Presence/Absence MDS
mds.16S.pa2<-metaMDS(decostand(data.16S3[,-c(1:8)],"pa" ),k=2,autotransform=FALSE, na.rm=TRUE)
mds.16S.pa2
#Stress=0.12, very good, 2-D is the way to go!

#Presence Absence

#base figure, no metadata
PA.Agg16S<-ggplot.NMDS2(mds.16S.pa2, (data.16S3$SoilFrac), colors.bw, shapes)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), legend.position=c(0.2,0.15),legend.key=element_blank(),axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))+
annotate("text", x=0.2,y=0.22, label="Bacteria", size=12)+labs(color="Aggregate Size", shape="Aggregate Size")+scale_fill_discrete(labels=c("Aggregate Size",">2 mm","<0.25 mm","1-2 mm","0.25-1 mm"))

PA.Agg16S



#Combined Figure:
grid.arrange(PA.Agg16S,PA.Agg1,ncol=2)




