#Elizabeth Bach
#COBS 2012 Fungi:  AMF root colonization Figures
#data collected by G. Narvaez and K. Murray
#July 8, 2014

rm(list-ls())
library(ggplot2)
library(plyr)
library(reshape)
library(gridExtra)

#Use "COBS_AMF_JulyOct.csv"
COBS.data<-read.csv(file.choose())

#summerize data to generate mean values for each cropping system in July and Oct.
COBS.data$Month<-factor(COBS.data$Month, levels=unique(as.character(COBS.data$Month)))
COBS.data$Crop<-as.factor(COBS.data$Crop)
COBS.data$inter<-paste(COBS.data$Crop, COBS.data$Month, sep=":")
COBS.data$inter<-factor(COBS.data$inter, levels=unique(as.character(COBS.data$inter)))
head(COBS.data)

#Total percent colonization
Pcol.sum<-ddply(COBS.data, .(inter), summarise, mean=(mean(P.col)), SE=(sd(P.col)/sqrt((length(P.col)-1))))
Pcol.names<-colsplit(Pcol.sum$inter, split=":", names=c("Crop", "Date"))
Pcol.matrix<-cbind(Pcol.sum, Pcol.names)
head(Pcol.matrix)
shapes<-c(21,22,24)
ylab<-expression(paste("Total AMF Root Colonization (%)"))

#Point graph, black & white, Corn>>Prairies P<0.0001, Prairies show Month*Date interaction, P=0.02
Pcol.2012<-ggplot()+geom_pointrange(data=Pcol.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position=c(0.9,0.1), legend.background=element_blank(), legend.text=element_text(size=16, face="bold"),legend.key=element_blank(),legend.title=element_blank(),axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Pcol.2012)

#for multipanel fig (no legend)
Pcol.2012a<-ggplot()+geom_pointrange(data=Pcol.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position="none",axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)+coord_cartesian(ylim=c(0,100))
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Pcol.2012a)


#Percent Vesicle infection
Pves.sum<-ddply(COBS.data, .(inter), summarise, mean=(mean(P.ves)), SE=(sd(P.ves)/sqrt((length(P.ves)-1))))
Pves.names<-colsplit(Pves.sum$inter, split=":", names=c("Crop", "Date"))
Pves.matrix<-cbind(Pves.sum, Pves.names)
head(Pves.matrix)
shapes<-c(21,22,24)
ylab<-expression(paste("Vesicle Colonization (%)"))

#Point graph, black & white, Prairies>>Corn P<0.0001, Prairies show Month*Date interaction, P<0.0001
Pves.2012<-ggplot()+geom_pointrange(data=Pves.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position=c(0.9,0.1), legend.background=element_blank(), legend.text=element_text(size=16, face="bold"),legend.key=element_blank(),legend.title=element_blank(),axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Pves.2012)

#for multipanel fig (no legend)
Pves.2012b<-ggplot()+geom_pointrange(data=Pves.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position="none",axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)+coord_cartesian(ylim=c(0,100))
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Pves.2012b)

#Percent Hyphae infection
Phyphae.sum<-ddply(COBS.data, .(inter), summarise, mean=(mean(P.hyphae)), SE=(sd(P.hyphae)/sqrt((length(P.hyphae)-1))))
Phyphae.names<-colsplit(Phyphae.sum$inter, split=":", names=c("Crop", "Date"))
Phyphae.matrix<-cbind(Phyphae.sum, Phyphae.names)
head(Phyphae.matrix)
shapes<-c(21,22,24)
ylab<-expression(paste("Hyphae Colonization (%)"))

#Point graph, black & white, Corn>>Prairies P<0.0001, Prairies show main effect of month, P=0.05
Phyphae.2012<-ggplot()+geom_pointrange(data=Phyphae.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position=c(0.9,0.1), legend.background=element_blank(), legend.text=element_text(size=16, face="bold"),legend.key=element_blank(),legend.title=element_blank(),axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Phyphae.2012)

#For multipanel figure (legend position adjusted)
Phyphae.2012c<-ggplot()+geom_pointrange(data=Phyphae.matrix, aes(x=Date, y=mean, ymin=mean-SE, ymax=mean+SE, shape=Crop), fill="black", color="black", size=1.25)+theme_bw()+
theme(aspect.ratio=1, axis.text=element_text(size=16, face="bold", colour="black"), axis.line=element_line(size=2), panel.border=element_blank(), panel.grid=element_blank(), legend.position=c(0.2,1.2), legend.background=element_blank(), legend.text=element_text(size=16, face="bold"),legend.key=element_blank(),legend.title=element_blank(),axis.ticks=element_line(size=2), strip.background=element_blank(),axis.title.x=element_blank(), axis.title.y=element_text(size=18, face="bold", colour="black"))+
ylab(ylab)+scale_shape_manual(values=shapes)+coord_cartesian(ylim=c(0,100))
#+annotate("text", label="C)", x=0.7, y=95, cex=8)
print(Phyphae.2012c)


Multipanel figure
grid.arrange(Pcol.2012a, Pves.2012b, Phyphae.2012c, ncol=3)

