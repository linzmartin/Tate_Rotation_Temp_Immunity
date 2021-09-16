############################
#Code written by Sadiq, with my edits
#############################
# clear workspace
rm(list = ls())


#load necessary packages
library(survival)
library(multcomp)
library(lme4)
library(MASS)


##########################################################
#Survival of Bt infection at different temperatures
##########################################################


library(openxlsx)
library(readxl)
#############survival
#read in file using import in Rstudio

#choose 1 of the following data sources:
#mydata<- read.xlsx("Survival_data.xlsx") #this is data for Sadiq's experiment 1
        #mydata<- X2020_09_29temperature_bt_infection #file name has an X in the beginning

mydata<- read.xlsx("2021-2-22_btu_temperature_survival_exp-2.xlsx") #data for Sadiq's second experiment
##############################
attach(mydata) #attach data to be able to call columns directly
summary(mydata) #overview to check if data was read in correctly
hours <- as.numeric(hours)
status <-as.numeric(status)
temp <-as.factor(temp)
sex<-as.factor(sex)


#plot survival curves

fit<-survfit(Surv(hours,status)~treatment+temp)
# curves for treatment (IS and Btu) and temperature (24, 30, 34)
plot(survfit(Surv(hours,status)~treatment+temp),  #decide which variables to show
     lty=c(1,1,1,2,2,2), lwd=2,                   #type of line used for each group
     col=c("cornflowerblue","darkgoldenrod","red"), #color used for each group
     xlim=c(7,22), ylim=c(0,1),                   #limits of both axes
     ylab="survival", xlab="hours post infection")#labels of the axes 

#add a legend
legend("bottomleft",                              #placement of legend
       c("24C","30C","34C","Btu","control"),      #labels in legend
       lty=c(1,1,1,1,2),cex=0.9,lwd=2,bty="n",    #line types, text size, line width
       col=c("cornflowerblue","darkgoldenrod","red","black","black")) #colors

fitsumm<-summary(fit)
fitsumm

#extract median survival times:
library(survminer)
surv_median(fit)
med_surv_times_extract_df<-as.data.frame(surv_median(fit))

#################################################################
####survival curves without IS control and seperate for both sexes

#create new data subset
btu <- mydata[ which(treatment=="BtU"),]
summary(btu) #check correctness of new subset
detach(mydata) #change attached dataset
attach(btu)
#plot survival curves by temperature and sex
        plot(survfit(Surv(hours,status)~sex+temp),  #decide which variables to show
        lty=c(1,1,1,2,2,2), lwd=2,                   #type of line used for each group
        col=c("cornflowerblue","darkgoldenrod","red"), #color used for each group
        xlim=c(7,22), ylim=c(0,1),                   #limits of both axes
        ylab="survival", xlab="hours post infection")#labels of the axes 

#add a legend
legend("bottomleft",                              #placement of legend
       c("24C","30C","34C","female","male"),      #labels in legend
       lty=c(1,1,1,1,2),cex=0.9,lwd=2,bty="n",    #line types, text size, line width
       col=c("cornflowerblue","darkgoldenrod","red","black","black")) #colors

btufit<-survfit(Surv(hours,status)~temp,data=btu)
btufit
survdiff(btufit)
survdiff(Surv(hours, status) ~ temp, data = btu)
####statistical analysis
#Cox proportional hazard


###testing for effect of temperature and sex on survival (including potential interaction)
#full model
model<-coxph(Surv(hours, status,)~(temp*sex))
summary(model)

#simplify model by excluding non significant terms (testing only the effect of temperature)
model<-coxph(Surv(hours, status,)~(temp))
summary(model)

#post hoc test (which temperatures are different from each other?)
summary(glht(model,mcp(temp="Tukey")))

#test whether hazards for each group are proportional and do not cross over time 
#(is it ok to perform Cox ph analysis?)
modelzph<-cox.zph(model)
modelzph
plot(modelzph[1])
abline(h=0,lty=3)

#reattach original data set if wanted
detach(btu)
attach(mydata)