#SIP model

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
###########################
SIS_between_host_model<- function (t, x, params) {
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
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> SISdata

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(S~time,data=SISdata,type="l",xlab="time",ylab="S")
plot(I~time,data=SISdata,type="l",xlab="time",ylab="I")
plot(P~time,data=SISdata,type="l",xlab="time",ylab="P")

mtext(outer=TRUE,line=-1.5,"Between host model")
par(op)

SISdata %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable", x="Time (days)")
