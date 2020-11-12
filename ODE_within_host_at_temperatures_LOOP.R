#Within host model for T. castaneum infection w/ Bt x Temp
#this is using the correct model from ODE_within_host_at_temperatures.R to show how T affects HIB
####################
# clear workspace
rm(list = ls())

###########################
#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(deSolve)
############################
Bt_Tcast_within_host_infection <- function (t, x, params) {
  #state variables
  H <- x[1] #Health Status
  I <- x[2] #Immune Effectors
  B <- x[3] #Bacteria
  state <-c(H, I, B)
  
  #parameters
  p <- params["p"]
  m <- params["m"]
  psi <- params["psi"]
  d <- params["d"]
  c <- params["c"]
  beta <- params["beta"]
  gamma <- params["gamma"]
  alpha <- params["alpha"]
  w <- params["w"]
  z <- params["z"]
  r <- params["r"]
  KH <- params["KH"]
  KB <- params["KB"]
  KI <- params["KI"]
  sigma<-params["sigma"]
  TCI <- params["TCI"]
  TII <- params["TII"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  dIdt <- (TCI*gamma*I*(KI-I)) + (TII*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z))
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}
#####################################

#times in ode is default unit of years - must make a fraction to get days
times <- seq(from=0,to=3/365,by=1/365/4)

parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
r=6000,d=0.003,c=0.005,sigma=500000)
#TCI=,TII=,TB=) 

#Bt -- Fraction of Topt Growth Rate:
#24C = 0.4821786
#30C = 0.7883752
#34C = 0.9553218

#Immune Genes - Fraction of Topt Rate of Rel expression
#HK 24C = 0.9718167
#HK 30C = 0.6631556
#HK 34C = 0.3677866

#BT 24C = 0.9156949
#BT 30C = 0.3954252
#BT 34C = 0.1603162

xstart<-c(H=1,I=1000,B=10000)

#for B. cererus, r = 0.010112 min^-1.
#therefore, 1/0.010112 = 98.89 min --> x60 (per hour) = 5933.54 --> 6000
#0.010112*60*6


##########################
ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(H~time,data=out,type="l",xlab="time",ylab="H")
plot(I~time,data=out,type="l",xlab="time",ylab="I")
plot(B~time,data=out,type="l",xlab="time",ylab="B")
mtext(outer=TRUE,line=-1.5,"Within host model")
par(op)

#######################
#graph the three together
#https://kinglab.eeb.lsa.umich.edu/480/nls/de.html#solving_odes_in_r
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  scale_y_log10()+
  labs(y="State Variable")+
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))

#####################################################################
###################################################################
##################################################################
#individ temps:

Bt_Tcast_within_host_infection <- function (t, x, params) {
  #state variables
  H <- x[1] #Health Status
  I <- x[2] #Immune Effectors
  B <- x[3] #Bacteria
  state <-c(H, I, B)
  
  #parameters
  p <- params["p"]
  m <- params["m"]
  psi <- params["psi"]
  d <- params["d"]
  c <- params["c"]
  beta <- params["beta"]
  gamma <- params["gamma"]
  alpha <- params["alpha"]
  w <- params["w"]
  z <- params["z"]
  r <- params["r"]
  KH <- params["KH"]
  KB <- params["KB"]
  KI <- params["KI"]
  sigma<-params["sigma"]
  TCI <- params["TCI"]
  TII <- params["TII"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  dIdt <- (TCI*gamma*I*(KI-I)) + (TII*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z))
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}
###############################
times <- seq(from=0,to=3/365,by=1/365/4)
#24C
parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
          gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
          r=6000,d=0.003,c=0.005,sigma=500000,
          TCI=0.9718167,TII=0.9156949,TB=0.4821786) 
xstart<-c(H=1,I=1000,B=10000)

ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out

vec<-c(rep(24,13))
out["temp"]<-vec
T24C<-out

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(H~time,data=out,type="l",xlab="time",ylab="H")
plot(I~time,data=out,type="l",xlab="time",ylab="I")
plot(B~time,data=out,type="l",xlab="time",ylab="B")
mtext(outer=TRUE,line=-1.5,"Within host model - T=24C")
par(op)

###################
#30C
times <- seq(from=0,to=3/365,by=1/365/4)

parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
          gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
          r=6000,d=0.003,c=0.005,sigma=500000,
          TCI=0.6631556,TII=0.3954252,TB=0.7883752) 
xstart<-c(H=1,I=1000,B=10000)

ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out


vec<-c(rep(30,13))
out["temp"]<-vec
T30C<-out

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(H~time,data=out,type="l",xlab="time",ylab="H")
plot(I~time,data=out,type="l",xlab="time",ylab="I")
plot(B~time,data=out,type="l",xlab="time",ylab="B")
mtext(outer=TRUE,line=-1,"Within host model - T=30C")
par(op)


##########################
#34C

times <- seq(from=0,to=3/365,by=1/365/4)
parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
          gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
          r=6000,d=0.003,c=0.005,sigma=500000,
          TCI=0.3677866,TII=0.1603162,TB=0.9553218) 
xstart<-c(H=1,I=1000,B=10000)

ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out


vec<-c(rep(34,13))
out["temp"]<-vec
T34C<-out

op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
          mar=c(3,3,2,2),mgp=c(2,1,0))

plot(H~time,data=out,type="l",xlab="time",ylab="H")
plot(I~time,data=out,type="l",xlab="time",ylab="I")
plot(B~time,data=out,type="l",xlab="time",ylab="B")
mtext(outer=TRUE,line=-1,"Within host model - T=34C")
par(op)
##################################

total <- rbind(T24C,T30C,T34C)

total %>%
  gather(variable,value,-time,-temp) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line(aes(linetype=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable")+
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))+
  scale_linetype_discrete(name="Temperature(°C)")



?facet_wrap()

total %>%
  gather(variable,value,-time,-temp)
