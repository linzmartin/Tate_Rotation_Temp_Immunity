#SIS Model - between hosts - non-obligate killer
#transmission depends on bacterial growth at time t
#writing loops to use w/in host model outputs
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
##############################
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
#################################################
#this is the 2-compartment between-host SIS model:
=======
###########################
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
SIS_between_host_model<- function (t, x, params) {
  #state variables
  S <- x[1] #Susceptible
  I <- x[2] #Infected
  state <-c(S,I)
  
  #parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  d <- params["d"]
  delta <- params["delta"]
  r<-params["r"]
<<<<<<< HEAD
=======
  #N<-params["N"]
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
  K<-params["K"]
  b<-params["b"]
  
  #model equations
  dSdt <- r*(1-(S+I)/K)-(b*beta*S*I)+(gamma*I)-(d*S)
  dIdt <- (b*beta*S*I)-(gamma*I)-(d+delta)*I
<<<<<<< HEAD
  
=======
  #N=S+I
>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4
  dndt <- c(dSdt,dIdt)
  list(dndt) #must be list format
}

<<<<<<< HEAD
###############
#create initial empty data frame to store model output
SIS_output.df <- data.frame() #must do this every time prior to running new models!
#run model at diff temps and store output:
for (i in 1:length(paramtable$Temp)){
  t50 <- paramtable[i,1]
  TB <- paramtable[i,3]
  
  xstart<-c(S=100,I=100) #start with 100 individuals in each compartment
  times<-seq(0,300,by=1) #300 days, by 1 day

  parms <- c (beta = 0.008, #Beta = contact rate * transmission probability
              #Beta = % of cases from overall pop that result in infection
              gamma = 1/t50, #gamma = 1/(infectious period)
              d=0.0005,#natural death rate
              delta=0.0005,#additional death rate due to infection
              r=0.8,#reproductive rate of susceptible beetles
              K=20000,#carrying capacity of live beetles
              b=TB) #coeff of opt bacterial growth
  #run model at each temp w/ parameters
  ode(
    func=SIS_between_host_model,
    y=xstart,
    times=times,
    parms=parms
  ) %>% 
    as.data.frame() -> SIS_out 
  
  SIS_out<-mutate(SIS_out,temp=paramtable[i,2]) #add temp to output
  SIS_output.df<-rbind(SIS_output.df,SIS_out) #store in data frame
}

#next, organize data
SIS_output.df<-SIS_output.df %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(SIS_output.df)
as.data.frame(SIS_output.df)
##################
#save data frame to excel spreadsheet in working directory:
write.xlsx(SIS_output.df, file = "SIS_2compartment_between_host_model_output.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)
########################
#Graph output at diff temps & save graphs to working directory:
modelgraph<-ggplot(SIS_output.df, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~variable,scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")
modelgraph

jpeg(file="SIS2C_24to34C.jpeg",width=600,height=500)
modelgraph
dev.off()

model_with_colors<-SIS_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(aes(color=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(size="legend",shape="legend")+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_colour_discrete(name="Temp(°C)")
model_with_colors

jpeg(file="SIS2C_24to34C_withcolors.jpeg",width=600,height=500)
model_with_colors
dev.off()

SIS_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line(aes(linetype=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_linetype_discrete(name="Temperature(°C)")
=======

#temp affects infectious period 
#at lower temps, takes longer to decline to 50% health
#at higher temps, reach 50% health faster, have longer IP
#all are dead by 2 days if infected & do not recover
#infectious_period =  2days - (time to 50% health) # infectious period
t50 <- 1.0 #hypothetical
infectious_period <- t50 #time infected to time of death
#Beta = contact rate * transmission probability
#Beta = % of cases from overall pop that result in infection
beta_value <- 0.8 #0.005
gamma_value <- 1 / infectious_period
d_value <- 0.5
delta_value <- 0.8
r<-1.3
K<-2000
b<-0.8
#gamma_value<-1/5 #1
#Ro <- beta_value / gamma_value
parms <- c (beta = beta_value, gamma = gamma_value, d=d_value,delta=delta_value,
            r=r,K=K,b=b)
Sinit <- 50 #1000      # susceptible hosts initial
Iinit <- 100 #1000           # infectious hosts initial
N <- Sinit + Iinit 
xstart <- c (S = Sinit/N, I = Iinit/N)
times <- seq (0, 20, by=1)
output = lsoda (xstart, times, SIS_between_host_model, parms)
plot (S ~ time, data = output, type='b', col = 'blue') 
plot (I ~ time, data = output, type='b', col = 'red')


ode(
  func=SIS_between_host_model,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> SISdata

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(S~time,data=SISdata,type="l",xlab="time",ylab="S")

plot(I~time,data=SISdata,type="l",xlab="time",ylab="I")

mtext(outer=TRUE,line=-1.5,"Between host model")
par(op)

SISdata %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable", x="Time (days)")


>>>>>>> bf161c5b7449982909bf8cffa040566eafd9f7a4

