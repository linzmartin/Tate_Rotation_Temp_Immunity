#SIR Practice code
###############################################################
#code from: https://kinglab.eeb.lsa.umich.edu/480/nls/de.html
#this will be practice manipulation of ODE / SIR equations in R
install.packages("deSolve")
library(deSolve)

#?ode
#Ode is the function for general solver for ODEs
'''
ode(y, times, func, parms, 
    method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
               "euler", "rk4", "ode23", "ode45", "radau", 
               "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration"),
    ...)
    '''

#needs initial values of state variables, times, and ODE function
#default method is Livermore solver of ordinary differential equations with automatic method switching (LSODA)
#################################################################

#SIR for closed Pop - ignore births and deaths
# B = mu = 0


closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  N <- S+I+R
  ## now code the model equations
  dSdt <- -beta*S*I/N
  dIdt <- beta*S*I/N-gamma*I
  dRdt <- gamma*I
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}
#order in this func must match expected ODE: time, func, params
#if all eqs w/in model do not depend on time, is autonomous
#ode solver expects values of ODE RHS to be elements of list

#Example: measles
parms <- c(beta=400,gamma=365/13) #initial conditions
times <- seq(from=0,to=60/365,by=1/365/4)
xstart <- c(S=999,I=1,R=0)

library(tidyverse)
ode(
  func=closed.sir.model,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out #store result in data frame

#plot results:
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (yr)',y='number of individuals')


##
#how does it depend on Beta transmission rate and infectious period?
betavals <- c(20,50,500) #vary Beta from 20 to 500
ipvals <- c(5,10,30)/365 #vary IP from 5 to 30 days
gammavals <- 1/ipvals

expand.grid(beta=betavals,gamma=gammavals)%>%
  group_by(beta,gamma) %>%
  do(
    {
      ode(func=closed.sir.model,y=xstart,times=times,
          parms=c(beta=.$beta,gamma=.$gamma)) %>%
        as.data.frame()
    }
  ) %>%
  ggplot(aes(x=time,y=I))+
  geom_line()+
  facet_grid(beta~gamma,scales='free_y',labeller=label_both)+
  theme_bw()


#######################################################

#for open pop
#have raw births B and deaths mu
open.sir.model <- function (t, x, params) {
  B <- params["B"]
  beta <- params["beta"]
  mu <- params["mu"]
  gamma <- params["gamma"]
  N <- x[1]+x[2]+x[3]
  dSdt <- B - beta*x[1]*x[2]/N - mu*x[1]
  dIdt <- beta*x[1]*x[2]/N - (mu+gamma)*x[2]
  dRdt <- gamma*x[2] - mu*x[3]
  list(c(dSdt,dIdt,dRdt))
}

#Example
parms <- c(B=20,mu=1/50,beta=400,gamma=365/13)

ode(
  func=open.sir.model,
  y=xstart,
  times=seq(from=0,to=25,by=1/365),
  parms=parms
) %>% 
  as.data.frame() -> out

#plot S, I, R against time
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line()+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  scale_y_log10()+
  guides(color=FALSE)+theme_bw()

#plot I against S
out %>%
  ggplot(aes(x=S,y=I))+geom_path()+
  scale_x_log10()+scale_y_log10()+theme_bw()

##################################################
#code from: https://cran.r-project.org/web/packages/shinySIR/vignettes/Vignette.html
install.packages("shinySIR")
library(shinySIR)
#can apply simple ODEs to this, doesn't need to be SIR

mySIRS <- function(t, y, parms) {
  with(as.list(c(y, parms)),{
    
    # Change in Susceptibles
    dS <- - beta * S * I + delta * R
    
    # Change in Infecteds
    dI <- beta * S * I - gamma * I
    
    # Change in Recovereds
    dR <- gamma * I - delta * R
    
    return(list(c(dS, dI, dR)))
  })
}

run_shiny(model = "SIRS (w/out demography)", 
          neweqns = mySIRS,
          ics = c(S = 9999, I = 1, R = 0),
          parm0 = c(beta = 5e-5, gamma = 1/7, delta = 0.1),
          parm_names = c("Transmission rate", "Recovery rate", "Loss of immunity"),
          parm_min = c(beta = 1e-5, gamma = 1/21, delta = 1/365),
          parm_max = c(beta = 9e-5, gamma = 1 , delta = 1))

default_models()

#####################################################
#system of ODE Equations
#code from http://rstudio-pubs-static.s3.amazonaws.com/32888_197d1a1896534397b67fb04e0d4899ae.html
library(deSolve)
library(ggplot2)
install.packages("reshape2")
library(reshape2)

# time sequence 
time <- seq(0, 50, by = 0.01)

# parameters: a named vector
parameters <- c(r = 2, k = 0.5, e = 0.1, d = 1)

# initial condition: a named vector
state <- c(V = 1, P = 3)

# R function to calculate the value of the derivatives at each time value
# Use the names of the variables as defined in the vectors above
lotkaVolterra <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dV = r * V - k * V * P
    dP = e * k * V * P - d * P
    return(list(c(dV, dP)))
  })
}
## Integration with 'ode'
out <- ode(y = state, times = time, func = lotkaVolterra, parms = parameters)

## Ploting
out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
library(reshape2)
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column

p <- ggplot(out.m, aes(time, value, color = variable)) + geom_point()
print(p)
