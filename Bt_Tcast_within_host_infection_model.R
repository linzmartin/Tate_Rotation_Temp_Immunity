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
  dHdt <- beta*B*exp(-H)-p*H-I*(1-H)+exp(-theta*H)
  dIdt <- gamma*(1-I/K) + alpha*B*(1-I/K)-I*(B*w + z)
  dBdt <- r*B*(1-B/K)-B*(d+c*I)
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt)
}


parms <-c(gamma=500,K=1,alpha=500,beta=510,p=0.1,theta=0.4,w=1,z=1,r=200,d=0.2,c=0.8) #initial conditions
times <- seq(from=0,to=60/365,by=1/365/4)
xstart<-c(H=1,I=50,B=1)



ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out

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

newout<- out %>% 
  gather(variable,value,-time)
newout<-as.data.frame(newout)
is.data.frame(newout)

  ggplot(aes(x=time,y=B,color=B))+
  #geom_line(size=2)+
  geom_point() +
  geom_line(size=1)+
  theme_classic()+
  labs(x='time (yr)',y='number of individuals')

ggplot(aes(data=newout,x=time,y=variable))
  ggplot(aes(data=newout,x=time,y=B,color=B))
