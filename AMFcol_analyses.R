#Elizabeth Bach
#COBS 2012 Fungi:  AMF root colonization data subset
#data collected by G. Narveaz and K. Murray
#Analysis for Fungal manuscript
#July 8, 2014

rm(list=ls())
library(lme4)
library(lmerTest)
library(bbmle)


#Use "COBS_AMF_JulyOct.csv"
COBS.data<-read.csv(file.choose())
#note, data is bimodal, with prairies and corn producing sepearate peaks
#Ran anova to discuss cropping system difference, then sub-set data to look at month differences in prairies and corn seperately

#For total AMF root colonization
hist(COBS.data$P.col)
#first is a null model with no effects other than the random effect for block
test.null<-lmer(P.col~1+(1|Block),data=COBS.data,REML=FALSE)

#next is the full model with all factorial combinations of month, LP, and CS
test.model.full<-lmer(P.col~Month*Crop+(1|Block), data=COBS.data, REML=FALSE)

#This model just includes the main effects
test.model.main<-lmer(P.col~Month+Crop+(1|Block), data=COBS.data, REML=FALSE)

#compare with AIC and anova for fitted models 
AICtab(test.null,test.model.full,test.model.main)
anova(test.null,test.model.full,test.model.main)
#full model is best fit, but AIC is only 0.4 lower than the main effects only model
anova(test.model.full)
#Month*Crop interaction is not significant, P=0.1, so may be justified in using main effects only
anova(test.model.main)
difflsmeans(test.model.main, ddf="Satterthwaite",type=3,method.grad="simple")

#For %vesicle infection
hist(COBS.data$P.ves)
test.null<-lmer(P.ves~1+(1|Block),data=COBS.data,REML=FALSE)
test.model.full<-lmer(P.ves~Month*Crop+(1|Block), data=COBS.data, REML=FALSE)
test.model.main<-lmer(P.ves~Month+Crop+(1|Block), data=COBS.data, REML=FALSE)

AICtab(test.null,test.model.full,test.model.main)
anova(test.null,test.model.full,test.model.main)
#Full model is best fit again, AIC score 5 lower than main effect only model
anova(test.model.full)
#There is a significant Month*Crop interaction: P=0.01
difflsmeans(test.model.full, ddf="Satterthwaite",type=3,method.grad="simple")

#For %Hyphae infection
hist(COBS.data$P.hyphae)
test.null<-lmer(P.hyphae~1+(1|Block),data=COBS.data,REML=FALSE)
test.model.full<-lmer(P.hyphae~Month*Crop+(1|Block), data=COBS.data, REML=FALSE)
test.model.main<-lmer(P.hyphae~Month+Crop+(1|Block), data=COBS.data, REML=FALSE)

AICtab(test.null,test.model.full,test.model.main)
anova(test.null,test.model.full,test.model.main)
#Main effects model is best fit, by 4, no Month*Crop interaction (P=0.9)
anova(test.model.main)
difflsmeans(test.model.main, ddf="Satterthwaite",type=3,method.grad="simple")

#Subset data to examine prairies and corn seperately
COBS.corn<-subset(COBS.data, Crop=="Corn")
head(COBS.corn)
COBS.prairies<-subset(COBS.data, Crop==c("prairie","prairieFert"))
head(COBS.prairies)

#Corn
hist(COBS.corn$P.col)
#Total colonization
test.model.main<-lmer(P.col~Month+(1|Block), data=COBS.corn, REML=FALSE)
anova(test.model.main)
#Month is NS
#vesicle infection
test.model.main<-lmer(P.ves~Month+(1|Block), data=COBS.corn, REML=FALSE)
anova(test.model.main)
difflsmeans(test.model.main, ddf="Satterthwaite",type=3,method.grad="simple")
#Oct>July, P=0.03
#Hyphae infection
test.model.main<-lmer(P.hyphae~Month+(1|Block), data=COBS.corn, REML=FALSE)
anova(test.model.main)
#Month is NS

#Prairies
hist(COBS.prairies$P.col)
#Total colonization
test.model.full<-lmer(P.col~Month*Crop+(1|Block), data=COBS.prairies, REML=FALSE)
anova(test.model.full)
difflsmeans(test.model.full, ddf="Satterthwaite",type=3,method.grad="simple")
#Month*Crop interaction P=0.02
#vesicle infection
test.model.full<-lmer(P.ves~Month*Crop+(1|Block), data=COBS.prairies, REML=FALSE)
anova(test.model.full)
difflsmeans(test.model.full, ddf="Satterthwaite",type=3,method.grad="simple")
#Month*Crop highly significant, P<0.0001
#hyphae infection
test.model.full<-lmer(P.hyphae~Month*Crop+(1|Block), data=COBS.prairies, REML=FALSE)
anova(test.model.full)
#Month*Crop NS
