#loop practice / writing:
####################
# clear workspace
rm(list = ls())
###########################
#load packages
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
library(dplyr)
library(tidyr)
library(ggplot2)
library(deSolve)
#################################
#creating temp data/ model for calculating fraction of Topt
#################################
#import the data and remove NAs
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Gene_Data_Modified") #Emily's data
# lookinga at ddCT2 (fold change) vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
Gene_data<-na.omit(Gene_data) #remove NAs first
###################################
#Flinn model for Immune Gene Data
fits <- Gene_data %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~nls_multstart(ddCT2~flinn_1991(temp = Temp, a, b, c),
                                               data = .x,
                                               iter = c(5,5,5),
                                               start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') - 10,
                                               start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') + 10,
                                               lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                               upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)))

#head(fits)
#calculate ddCT2 based on model at temps to find max ddCT and %s
neededtemps <- Gene_data %>%
  do(., data.frame(Temp = c(seq(from=24,to=34,by=0.1)), stringsAsFactors = FALSE))

predsforcalc <- fits %>%
  mutate(., p = map(fit, augment, newdata = neededtemps)) %>%
  unnest(p) %>%
  group_by(., Treatment) %>%
  rename(., ddCT2 = .fitted) %>%
  ungroup()

ddCT2atTemp<-select(predsforcalc,-data,-fit)
ddCT2atTemp<-as.data.frame(ddCT2atTemp) 

rate_data <- data.frame() #create empty data frame
#then, calculate constitutive and induced rates 
for (i in 1:length(ddCT2atTemp$Temp)){
  HK <- filter(ddCT2atTemp, Treatment=="heat killed")
  NI <- filter(ddCT2atTemp, Treatment=="non-injected")
  BtU <- filter(ddCT2atTemp, Treatment=="BtU")
  
  constitutive_rate<-(HK$ddCT2-NI$ddCT2)/4 #time frame is 4 hours
  induced_rate<-(BtU$ddCT2-NI$ddCT2)/4 #time frame is 4 hours
  Tempofi<-ddCT2atTemp$Temp
  
  data_frame <- data.frame("Temp"=Tempofi, constitutive_rate,induced_rate)
  rate_data <- rbind(rate_data, data_frame) #add rates empty data frame after each calculation
}
#get rid of repeated data:
length(neededtemps[,1])
rate_data<-rate_data[1:length(neededtemps[,1]),] #now, have T 24-24 only once


#find ddCT2rate at Topt 25.2 = max rate
constit_max_rate_Topt <- subset(rate_data, Temp==25.2,select = c("constitutive_rate"))
constit_max_rate_Topt<-constit_max_rate_Topt[,1]

induced_max_rate_Topt<-subset(rate_data, Temp==25.2,select = c("induced_rate"))
induced_max_rate_Topt<-induced_max_rate_Topt[,1]

#then find fraction: ddCT2 rate at Ti / rate at Topt
fractions<-data.frame()
for (i in 1:length(rate_data$Temp)){
  fraction_induced<-(rate_data$induced_rate/induced_max_rate_Topt) %>% round(.,digits=3)
  fraction_constitutive<-(rate_data$constitutive_rate/constit_max_rate_Topt) %>% round(.,digits=3)
  
  df <- data.frame("Temp"=Tempofi, fraction_constitutive,fraction_induced)
  fraction_data <- rbind(fractions,df) #add rates empty data frame after each calculation
}
#get rid of repeated data:
fraction_data<-fraction_data[1:length(neededtemps[,1]),]
###########################################################

#now, use a loop to:
#(1) select fractions as parameters for each temp 
#(2) run model at each temp w/ paramters
#(3) store in matrix

#then plot models on graph

##############
#this is the function
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
#need to find these params from data #TCI=,TII=,TB=) 


