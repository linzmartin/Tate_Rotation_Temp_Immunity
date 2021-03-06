#SIS Model - between hosts
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
  state <-c(S,I)
  
  #parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  d <- params["d"]
  delta <- params["delta"]
  r<-params["r"]
  #N<-params["N"]
  K<-params["K"]
  b<-params["b"]

  #model equations
  dSdt <- r*(1-(S+I)/K)-(b*beta*S*I)+(gamma*I)-(d*S)
  dIdt <- (b*beta*S*I)-(gamma*I)-(d+delta)*I
  #N=S+I
  dndt <- c(dSdt,dIdt)
  list(dndt) #must be list format
}


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



