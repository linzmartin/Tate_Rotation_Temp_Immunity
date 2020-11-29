<<<<<<< HEAD
#SI with compartment P model

=======
#SIP model

#transmission depends on bacterial growth at time t
#writing loops to use w/in host model outputs
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
###########################
# clear workspace
rm(list = ls())
#set wd to your project folder
setwd("C:/Users/linzm/Documents/Tate_Lab_Rotation/Tate_Rotation_Temp_Immunity")
###########################
#load packages
library (deSolve)
library(ggplot2)
library(tidyr)
library(nlstools)
library(dplyr)
library(tidyr)
<<<<<<< HEAD
library(openxlsx)
library(readxl)
library(ggplot2)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(purrr)
###########################
#import the data
withinhostmodeldata <- read.xlsx("HIB_within_host_model_output_longversion.xlsx") #import within host model HIB data
Health <- subset(withinhostmodeldata, variable=="H") #subset to only get health variable

tempparamdata<-read.xlsx("Temp_dependent_params_withinhost.xlsx") #import Temp param data
###########################################
#extract time values for time to 50% health from data frame:
time_values<-data.frame()
for (i in unique(Health$temp)){
  Healthsub<-subset(Health,temp==i)
  health50percent<-which(abs(Healthsub$value-0.5)==min(abs(Healthsub$value-0.5)))
  timetoH50<-Healthsub$time[health50percent]
  times<-cbind(timetoH50,temp=i)
  time_values<-rbind(time_values,times)
}
time_values$daystoH50<-(time_values$timetoH50)/0.0027397 #convert time from fraction of years to days (24 hrs = 0.0027397 yrs)

#combine temp, TB, and time to 50% health into one data frame for use:
paramtable<-as.data.frame(cbind(Time50H=time_values$daystoH50,Temp=time_values$temp,TB=tempparamdata$TB))

savingparams<-as.data.frame(cbind(Temp=time_values$temp, t50_value=time_values$daystoH50,TB=tempparamdata$TB ))
write.xlsx(savingparams, file = "Temp_dependent_params_betweenhost.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)
#################################################
#3 compartment between host model:
SIS_between_host_model_with_compartmentP<- function (t, x, params) {
=======
###########################
SIS_between_host_model<- function (t, x, params) {
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
  #state variables
  S <- x[1] #Susceptible
  I <- x[2] #Infected
  P <- x[3] #Dead and spore producing
  state <-c(S,I,P)
  
  #parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  d <- params["d"]
  delta <- params["delta"]
  r<-params["r"]
<<<<<<< HEAD
  K<-params["K"]
  g<-params["g"]
  #Kp<-params["Kp"]
  mu<-params["mu"]
  TB<-params["TB"] #bacterial growth rate fraction (from within host / Ratkowsky models)
  #f<-params["f"] #fraction of topt Median Survival Time / bacterial growth rate fraction?

  #model equations
  #### developing model: #####
  #dSdt <- r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S)
  #dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I)
  #dPdt <- g*((d+delta)*I)*(1-P/Kp)-(beta*P)
  
  #dSdt <- r*(1-(S+I)/K)-(beta*S*g*P*(1-P/Kp))+(gamma*I)-(d*S) #susceptible beetles
  #dIdt <- (beta*S*g*P*(1-P/Kp))-(gamma*I)-(f*(d+delta)*I) #infected beetles
  #dPdt <- (f*(d+delta)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  ##############
  
  #better ODEs for model:
  dSdt <- r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S) #susceptible beetles
  dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I) #infected beetles
  dPdt <- (g*TB*(d+delta)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  
  dndt <- c(dSdt,dIdt,dPdt)
  list(dndt) #must be list format
}
#######################################################
#create initial empty data frame to store model output
SIP_output.df <- data.frame()#must do this every time prior to running new models!
#run model at diff temps and store output:
gammaval_output <-data.frame()
for (i in 1:length(paramtable$Temp)){
  t50 <- paramtable[i,1]
  TB <- paramtable[i,3]
  
  xstart<-c(S=100,I=100,P=0) #start with 100 individuals in each compartment
  times<-seq(0,300,by=1) #300 days, by 1 day
  
  parms <- c (beta = 0.008, #Beta = contact rate * transmission probability
              #Beta = % of cases from overall pop that result in infection
              gamma = 1/t50, #gamma = 1/(infectious period)
              d=0.0005,#natural death rate
              delta=0.0005,#additional death rate due to infection
              r=0.8,#0.8#reproductive rate of susceptible beetles
              K=20000,#carrying capacity of live beetles
              g=100, #sporulation rate of dead infected beetle #this has huge impact on shape
              #Kp=1000, #carrying capacity of spores per infected beetle
              mu=0.05, #decay of dead beetles producing spores
              TB = TB) #fraction of Topt bacterial growth

  #run model at each temp w/ parameters
  ode(
    func=SIS_between_host_model_with_compartmentP,
    y=xstart,
    times=times,
    parms=parms
  ) %>% 
    as.data.frame() -> SIP_out 
  
  SIP_out<-mutate(SIP_out,temp=paramtable[i,2]) #add temp to output
  SIP_output.df<-rbind(SIP_output.df,SIP_out) #store in data frame
  
  #save gamma values as well for reference
  gammavalue<-(1/t50)
  gammaout<-data.frame(gammavalue)
  gammaout<-cbind(gammaout, temp=paramtable[i,2]) #add temp
  gammaval_output <- rbind(gammaval_output,gammaout)
}

#next, organize data
SIP_output.df<-SIP_output.df %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(SIP_output.df)
as.data.frame(SIP_output.df)

##################
#save data frame to excel spreadsheet in working directory:
write.xlsx(SIP_output.df, file = "SIP_3compartment_between_host_model_output.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

write.xlsx(gammaval_output, file="Temp_dependent_gamma_values_between_host_model.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

SIP_subset<-subset(SIP_output.df,variable=="S"|variable=="I")
########################
#Graph output at diff temps & save graphs to working directory:
modelgraph<-ggplot(SIP_output.df, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~(variable),scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")
modelgraph

#graph SI only:
modelgraphSIonly<-ggplot(SIP_subset, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~(variable),scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")
modelgraphSIonly

jpeg(file="SIP3C_24to34.jpeg",width=600,height=500)
modelgraph
dev.off()

model_with_colors<-SIP_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(aes(color=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(size="legend",shape="legend")+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_colour_discrete(name="Temp(°C)")
model_with_colors

jpeg(file="SIP3C_24to34_withcolors.jpeg",width=600,height=500)
model_with_colors
dev.off()

SIP_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line(aes(linetype=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_linetype_discrete(name="Temperature(°C)")


############################################################
#############################################################
#individ model w/o effects of temp:
ode(
  func=SIS_between_host_model_with_compartmentP,
=======
  #N<-params["N"]
  K<-params["K"]
  #b<-params["b"]
  g<-params["g"]
  Kp<-params["Kp"]
  #f<-params["f"]
  
  #model equations
  dSdt <- r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S)
  dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I)
  #dPdt <- g*((d+delta)*I)*(1-P/Kp)-(f*beta*S*P)
  dPdt <- g*((d+delta)*I)*(1-P/Kp)-(beta*P)
  dndt <- c(dSdt,dIdt,dPdt)
  list(dndt) #must be list format
}


#temp affects infectious period 
#at lower temps, takes longer to decline to 50% health
#at higher temps, reach 50% health faster, have longer IP
#all are dead by 2 days if infected & do not recover
#infectious_period =  2days - (time to 50% health) # infectious period
t50 <- 0.8 #days out of 365 #hypothetical
infectious_period <- t50 #time infected to time of death

#Sinit <- 100 #1000      # susceptible hosts initial
#Iinit <- 100 #1000           # infectious hosts initial
#N <- Sinit + Iinit 
#Pinit
#Pinit<-0
#xstart <- c (S = Sinit/N, I = Iinit/N,P = Pinit)
parms <- c (beta = 0.008, #Beta = contact rate * transmission probability
            #Beta = % of cases from overall pop that result in infection
            gamma = 1/1.1, #gamma = 1/(infectious period) #has impact
            d=0.0005,#death rate
            delta=0.0003,#additional death rate due to infection
            r=0.8,#reproductive rate of susceptible beetles #little impact
            K=20000,#carrying capacity of live beetles
            g=200, #sporulation rate of dead infected beetles #this has huge impact on shape
            Kp=1500) #carrying capacity of spores #little impact
  #decay of spores as they infect new hosts

xstart<-c(S=100,I=100,P=0)
#times <- seq(0, 50/365, by=1/365)
times<-seq(0,25,by=1)
#times<-seq(0,25,by=0.5)
output = lsoda (xstart, times, SIS_between_host_model, parms)
plot (S ~ time, data = output, type='b', col = 'blue') 
plot (I ~ time, data = output, type='b', col = 'red')
plot(P~time,data=output,type="b",col="orange")


ode(
  func=SIS_between_host_model,
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
  y=xstart,
  times=times,
  parms=parms
) %>% 
<<<<<<< HEAD
  as.data.frame() -> SISdatawithP
=======
  as.data.frame() -> SISdata
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

<<<<<<< HEAD
plot(S~time,data=SISdatawithP,type="l",xlab="time",ylab="S")
plot(I~time,data=SISdatawithP,type="l",xlab="time",ylab="I")
plot(P~time,data=SISdatawithP,type="l",xlab="time",ylab="P")
=======
plot(S~time,data=SISdata,type="l",xlab="time",ylab="S")
plot(I~time,data=SISdata,type="l",xlab="time",ylab="I")
plot(P~time,data=SISdata,type="l",xlab="time",ylab="P")
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4

mtext(outer=TRUE,line=-1.5,"Between host model")
par(op)

<<<<<<< HEAD
SISdatawithP %>%
=======
SISdata %>%
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable", x="Time (days)")
