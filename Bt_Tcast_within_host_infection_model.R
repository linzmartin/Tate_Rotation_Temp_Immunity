library(dplyr)
library(tidyr)
library(ggplot2)
library(deSolve)

Bt_Tcast_within_host_infection <- function (t, x, params) {
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
}
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
