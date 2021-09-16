# Tate Lab - Temperature x Immunity Models
This project aims to evaluate the effects of temperature on flour beetle immunity and spread of disease in a flour beetle population over time.
This project creates (1) a within-host model for Bt infection in T. castaneum and (2) a between-host model for infection in a T. castaneum population. Then, these models incorporate changes in temperature and thermal optimum differences between host and the Bt bacterial infection to explore the nature of these relationships over time.
------------------
This was developed by the Tate Lab at Vanderbilt University, by Lindsay E. Martin under the supervision of Dr. Ann Tate.  Data is provided by current and former members of the Tate Lab, as well as pulled from prior publications (where noted). 
------------------

## Data files:
Gene expression data: "Gene_temp_data.xlsx", sheet="Rate_Data""
  This includes microbe-independent rates and microbe-dependent rates of gene  expression for immune genes at different temperatures at 20, 24, 30, and 24C.
  
Repro_data <- read_xlsx("Tcast_reproduction_Park1948data.xlsx") #reproduction v. temp data (Park 1948)


Germination_data <- read.csv("Knaysi_1964_Fig1_Graph_grabber_data.csv", stringsAsFactors = FALSE) 

#####  

## Required packages:
library(readxl)
library(ggplot2)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(deSolve)
library(xlsx)
library(DescTools)
library(reshape2)
library(metR)



## R code and models:


## Functions:

Ratkowsky bacterial growth model:
Ratkowsky <- function(temp, Tmin, Tmax, b, c){
  rate <- (b*(temp-Tmin)*(1-exp(c*(temp-Tmax))))^2
  return(rate)
}

within host model:
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
  
  TMD <- params["TMD"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  dIdt <- (TMD*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z)) #new changes: get rid of microbe-independent rate & TMI parameter
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}

Bt reproduction / spore germination:
Repro_gaussian <- function(temp, rmax, topt, a){
  reprorate <- rmax*exp(-0.5*(abs(temp-topt)/a)^2)
  return(reprorate)
}



## Output data files used for further analysis:



## Versions:


## Troubleshooting Tips:

If the function select() produces the error "unused argument", try re-running under the explicit package of dplyr via dplyr::select() - this may be an error where select is being run under the MASS package instead of dplyr.  Alternatively, try refreshing and only loading dplyr and not the MASS package.


