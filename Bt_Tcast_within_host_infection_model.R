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
  theta <- params["theta"]
  gamma <- params["gamma"]
  alpha <- params["alpha"]
  w <- params["w"]
  z <- params["z"]
  r <- params["r"]
  K <- params["K"]

  #model equations
  #dHdt <- (beta*B*exp(-H))-(p*H)-(I*(1-H))+exp(-theta*H) #H not limited to 1
  #dHdt <- (beta*B*(1-exp(H)))-(p*H)-(I*(1-H))+exp(-theta*H) #this works
  #dHdt <- (beta*B*(-exp(H)))-(p*H)-(I*(1-H))+exp(-theta*H) #H becomes negative
  #dHdt <- (beta*B*(exp(H)))-(p*H)-(I*(1-H))+exp(-theta*H) #doesn't work at all
  
  
  #dHdt <- (beta*B*(1-exp(H)))-(p*H)-(I*(1-H))+exp(-theta*H)
  dHdt <- (beta*B*(1-exp(H/B)))-(p*H)-(I*(1-H))+exp(-theta*H) 
  
  #dIdt <- gamma*(1-I/(1-exp(-K))) + alpha*B*(1-I/(1-exp(-K)))-I*(B*w + z)
  #dIdt <- gamma*(1-I/(K*1/B)) + alpha*B*(1-I/(K*1/B))-I*(B*w + z) #works but is harsh
  #dIdt <- gamma*(1-I/(K)) + alpha*B*(1-I/(K-B))-I*(B*w + z)
  #dIdt <- gamma*(1-I/K) + alpha*B*(1-I/K)-I*(B*w + z) #original 
  
  dIdt <- gamma + (alpha*B*(1-I/K))-(I*(B*w + z)) #this might work??
  
  dBdt <- r*B*(1-B/K)-B*(d+c*I)
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt)
}
parms <-c(gamma=500,K=2e5,alpha=50,
          beta=510,p=0.2,theta=0.8,
          w=1,z=10,r=200,d=2,c=0.5) #initial conditions
#times <- seq(from=0,to=3,by=1/365/4)
#times <- seq(from=0,to=3/365,by=1/365/4)
times <- seq(from=0,to=60/365,by=1/365/4)
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
