Bt_Tcast_within_host_infection <- function (t, x, params) {
  #state variables
  H <- x[1]
  #L <- x[2]
  B <- x[3]
  X <- x[4]
  #state <-c(H, L, B, X)
  state <-c(H, B, X)
  
  #parameters
  M <- params["M"]
  d <- params["d"]
  Toll <- params["Toll"]
  IMD <- params["IMD"]
  fB <- params["fB"]
  c <- params["c"]
  beta <- params["beta"]
  y <- params["y"]
  s <- params["s"]

  #model equations
  dHdt <- M*s-M*d*H-beta*H*X
  #dLdt <- beta*H*X
  dBdt <- M*fB-M*Toll*IMD*c*B
  dXdt <- y*M*fB-M*Toll*IMD*c*X
  
  #single vector as list
  #dzdt <- c(dHdt,dLdt,dBdt,dXdt)
  dzdt <- c(dHdt,dBdt,dXdt)
  list(dzdt)
}


parms <-c(s=10,M=10,d=10,Toll=15,IMD=19,fB=74,c=1,beta=110,y=50) #initial conditions
times <- seq(from=0,to=30/365,by=1/365/4)
xstart<-c(H=1000,B=5e10,X=5e10)

ode(
  func=Bt_Tcast_within_host_infection,
  y=xstart,
  times=times,
  parms=parms
) %>% 
  as.data.frame() -> out

out


