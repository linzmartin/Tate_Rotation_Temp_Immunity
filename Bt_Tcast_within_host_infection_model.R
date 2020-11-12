#Within host model for T. castaneum infection w/ Bt x Temp
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
  H <- x[1]
  I <- x[2]
  B <- x[3]
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
  

  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  #dIdt <- (gamma*I*(1-I/KI)) + (alpha*I*B*(1-I/KI))-(I*(B*w + z))
  dIdt <- (gamma*I*(KI-I)) + (alpha*I*B*(KI-I))-(I*(B*w + z))
  dBdt <- r*B*(1-B/KB)-B*d-c*I*B
    #dBdt <- r*B*(1-B/KB)*(KB/(1+exp(-0.01*B*0.5)))-B*d-c*I*B
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt)
}

#parms take 1 - works for model but I is steep
parms <-c(psi=0.8, KH=1, beta=50,p=0.2,m=0.2,
          gamma=100,KI=2000,alpha=100,KB=3000000,w=1,z=10,
          r=200,d=20,c=10) #initial conditions
times <- seq(from=0,to=3/365,by=1/365/4)
xstart<-c(H=1,I=100,B=10000)


##
#parms take 1 - works for model but I is steep
parms <-c(psi=0.8, KH=100, beta=5,p=0.2,m=0.2,
          gamma=50,KI=20,alpha=50,KB=30000,w=1,z=10,
          r=20,d=2,c=0.5) #initial conditions
times <- seq(from=0,to=3/365,by=1/365/4)
xstart<-c(H=1,I=100,B=100)
#times <- seq(from=0,to=3,by=1/365/4)
#times<-seq(from=0,to=24,by=2)
times
#times <- seq(from=0,to=60/365,by=1/365/4)
#times <- seq(from=0, to=3, by=1/3/12)
#times<-seq(0,4/24,by=1/24/60)
#times
xstart<-c(H=1,I=100,B=100)

xstart<-c(H=1,I=100,B=10000)
##########################



parms <-c(psi=0.8, KH=1, beta=50,p=0.2,m=0.2,
          gamma=50,KI=200,alpha=50,KB=300000,w=1,z=10,
          r=200,d=2,c=0.5) #initial conditions
times <- seq(from=0,to=3/365,by=1/365/4)
#times <- seq(from=0, to=48, by=1/48)

#times
xstart<-c(H=1,I=100,B=1000)
#########################################
########
#take 2
parms <-c(psi=0.8, KH=1, beta=5,p=0.2,m=0.2,
          gamma=50,KI=200,alpha=50,KB=30,w=1,z=10,
          r=20,d=2,c=0.5) #initial conditions
#times <- seq(from=0,to=3,by=1/365/4)
#times <- seq(from=0,to=3/365,by=1/365/4)
#times <- seq(from=0,to=60/365,by=1/365/4)
#times <- seq(from=0, to=3, by=1/3/12)
#times <- seq(from=0, to=3, by=1/3/12)
times<-seq(from=0,to=36,by=12)
times<-seq(from=0,to=36,by=4)
times<-seq(from=0,to=36,by=1)
#times
xstart<-c(H=1,I=10,B=1e5)

#################


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














out

#plot results:
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  #geom_line(size=2)+
  geom_point() +
  geom_line(size=1)+
  theme_classic()+
  labs(x='time (yr)',y='number of individuals')
##
#or graph by:
newout<- out %>% 
  gather(variable,value,-time)
newout<-as.data.frame(newout)
#is.data.frame(newout)

#individual graphs:
H<-filter(newout,variable=="H")
B<-filter(newout,variable=="B")
I<-filter(newout,variable=="I")

ggplot(H, aes(x=time,y=value)) +
  geom_point()
ggplot(B, aes(x=time,y=value)) +
  geom_point()
ggplot(I, aes(x=time,y=value)) +
  geom_point()
  
#################################
op <- par(fig=c(0,1,0,1),mfrow=c(2,2),
            mar=c(3,3,2,2),mgp=c(2,1,0))
plot(H~time,data=out,type="l",main="H",xlab="time",ylab="dH/dt")
plot(I~time,data=out,type="l",main="I",xlab="time",ylab="dI/dt")
plot(B~time,data=out,type="l",main="B",xlab="time",ylab="dB/dt")
mtext(outer=TRUE,line=-1,"Within host model")
par(op)

#######################
#graph the three together
#https://kinglab.eeb.lsa.umich.edu/480/nls/de.html#solving_odes_in_r
out %>%
    gather(variable,value,-time) %>%
    ggplot(aes(x=time,y=value,color=variable))+geom_line()+
    facet_wrap(~variable,scales="free_y",ncol=1)+
    scale_y_log10()+
    guides(color=FALSE)+theme_bw()




###########################################

Bt_Tcast_within_host_infection <- function (time, x, params) {
  params<-as.list(c(x,params))
  #state variables
  H <- x[1]
  I <- x[2]
  B <- x[3]
  state <-c(H, I, B)
  
  #parameters
  p <- params["p"]
  d <- params["d"]
  c <- params["c"]
  beta <- params["beta"]
  gamma <- params["gamma"]
  alpha <- params["alpha"]
  w <- params["w"]
  z <- params["z"]
  r <- params["r"]
  psi <- params["psi"]
  KH <- params["KH"]
  KB <- params["KB"]
  KI <- params["KI"]
  m <- params["m"]
  mu <- params["mu"]
  sigma1 <- params["sigma1"]
  
  with(params,{
  #model equations
  #dHdt <- (beta*B*(1-exp(H)))-(p*H)-(I*(1-H))+exp(-theta*H) original
  #dHdt <- (beta*B*(1-exp(H/B)))-(p*H)-(I*(1-H))+exp(-theta*H) 
  dHdt <- (psi*H*(KH-H))-(beta*B*(1-exp(H)))-(p*H)-(I*m*(1-H))
  #dHdt <- (psi*H*(1-KH))-(beta*B*(1-exp(H)))-(p*H)-(I*m*(1-H))
  
  #dIdt <- gamma*(1-I/K) + alpha*B*(1-I/K)-I*(B*w + z) #original
  #dIdt <- (gamma) + (alpha*B*(1-I/K))-(I*(B*w + z))
  dIdt <- (gamma*I*(KI-1)) + (alpha*I*(KI-1))*(KB/(1+exp(-mu*(B-sigma1))))-(I*(B*w + z)) #this might work??
  
  dBdt <- r*B*(1-B/KB)-B*(d+c*I)
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt)
  })
}


print(as.list(names(parms)))


library(deSolve)
#init <- c(S = N-Infected[1], I = Infected[1], R = 0)
parms <-c(psi=20, KH=10, beta=51,p=0.2,
          gamma=5,KI=20,alpha=5,KB=30,mu=1,sigma1=10, w=1,z=10,
          r=20,d=2,c=0.5) #initial conditions

RSS <- function(parms) {
  names(parms) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}

Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
##      beta     gamma 
## 0.6746089 0.3253912

t <- 1:70 # time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
col <- 1:3 # colour

matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col)
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")
## Warning in xy.coords(x, y, xlabel, ylabel, log = log): 1 y value <= 0
## omitted from logarithmic plot

points(Day, Infected)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.05)
title("SIR model 2019-nCoV China", outer = TRUE, line = -2)






parms <-c(psi=20, KH=10, beta=51,p=0.2,
          gamma=5,KI=20,alpha=5,KB=30,mu=1,sigma1=10, w=1,z=10,
          r=20,d=2,c=0.5) #initial conditions
#times <- seq(from=0,to=3,by=1/365/4)
times <- seq(from=0,to=3/365,by=1/365/4)
#times <- seq(from=0,to=60/365,by=1/365/4)
#times <- seq(from=0, to=3, by=1/3/12)

times
xstart<-c(H=1,I=100,B=100)


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
mtext(outer=TRUE,line=-2,"Within host model")
par(op)


#################################

#BTmid
##this is the function:
Ratkowsky <- function(temp, Tmin, Tmax, b, c){
  rate <- (b*(temp-Tmin)*(1-exp(c*(temp-Tmax))))^2
  return(rate)
}

Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
Btoptrate
Bt24rate<-Ratkowsky(temp=24,Tmin=6,Tmax=49,b=0.004,c=0.14)
Bt24rate
Bt30rate<-Ratkowsky(temp=30,Tmin=6,Tmax=49,b=0.004,c=0.14)
Bt30rate
Bt34rate<-Ratkowsky(temp=34,Tmin=6,Tmax=49,b=0.004,c=0.14)
Bt34rate

