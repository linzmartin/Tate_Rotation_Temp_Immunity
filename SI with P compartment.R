#SI with compartment P model

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
#################################################
#3 compartment between host model:
SIS_between_host_model_with_compartmentP<- function (t, x, params) {
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
  K<-params["K"]
  g<-params["g"]
  Kp<-params["Kp"]
  mu<-params["mu"]
  f<-params["f"] #fraction of topt Median Survival Time / bacterial growth rate fraction?

  #model equations
  #dSdt <- r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S)
  #dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I)
  #dPdt <- g*((d+delta)*I)*(1-P/Kp)-(beta*P)
  dSdt <- r*(1-(S+I)/K)-(beta*S*g*P*(1-P/Kp))+(gamma*I)-(d*S) #susceptible beetles
  dIdt <- (beta*S*g*P*(1-P/Kp))-(gamma*I)-(f*(d+delta)*I) #infected beetles
  dPdt <- (f*(d+delta)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  
  dndt <- c(dSdt,dIdt,dPdt)
  list(dndt) #must be list format
}
#######################################################
#create initial empty data frame to store model output
SIP_output.df <- data.frame() #must do this every time prior to running new models!
#run model at diff temps and store output:
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
              Kp=1000, #carrying capacity of spores per infected beetle
              mu=0.05, #decay of dead beetles producing spores
              f = TB) #fraction of Topt bacterial growth

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
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> SISdatawithP

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(S~time,data=SISdatawithP,type="l",xlab="time",ylab="S")
plot(I~time,data=SISdatawithP,type="l",xlab="time",ylab="I")
plot(P~time,data=SISdatawithP,type="l",xlab="time",ylab="P")

mtext(outer=TRUE,line=-1.5,"Between host model")
par(op)

SISdatawithP %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable", x="Time (days)")
