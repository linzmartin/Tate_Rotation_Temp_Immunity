#Setting up data from within host model for between host model use
#saving extra values for HIB model at smaller time intervals
#saving Temp parameters
############################################
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
library(openxlsx)
library(readxl)
library(ggplot2)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(purrr)
library(tidyr)
library(nlstools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(deSolve)
##############################
#########################################################
#import Data & remove NA values
Gene_data <- read.xlsx("Gene_temp_data.xlsx", sheet="Rate_Data") #Emily's data
Gene_data<-na.omit(Gene_data) #remove NAs first

###############################################
#Flinn model for Immune Gene Expression Rates
fits <- Gene_data %>%
  group_by(., Rate_Type) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                               data = .x,
                                               iter = c(5,5,5),
                                               start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                               start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                               lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                               upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)))

head(fits)

#########################
#Find max gene expression Rates and %s of optimums

##first, we need to calculate the max values:
neededtemps <- Gene_data %>%
  do(., data.frame(Temp = c(seq(from=24,to=34,by=0.1)), stringsAsFactors = FALSE))

predsforcalc <- fits %>%
  mutate(., p = map(fit, augment, newdata = neededtemps)) %>%
  unnest(p) %>%
  group_by(., Rate_Type) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

RateatTemp<-select(predsforcalc,-data,-fit)
RateatTemp<-as.data.frame(RateatTemp) 

#max temps:
select(fits, data, fit) 

indep_topt <- calc_params(fits$fit[[1]]) %>% select(topt) #microbe independent
indep_Topt <- as.numeric(indep_topt) %>% round(.,digits=1)
dep_topt <- calc_params(fits$fit[[2]]) %>% select(topt) #microbe dependent
dep_Topt <- as.numeric(dep_topt) %>% round(.,digits=1)

#find rate at Topt = max rate
indep_max_rate_Topt <- subset(RateatTemp, Rate_Type=="Microbe_Independent_Rate" & Temp==indep_Topt,select = c("Rate"))
indep_max_rate_Topt<-indep_max_rate_Topt[,1]
dep_max_rate_Topt<-subset(RateatTemp, Rate_Type=="Microbe_Dependent_Rate" & Temp==dep_Topt,select = c("Rate"))
dep_max_rate_Topt<-dep_max_rate_Topt[,1]

##########
#next, we can calculate predictions for every 2C (or degrees of interest):
tempsforgraphing<-Gene_data %>%
  do(., data.frame(Temp = c(seq(from=24,to=34,by=2)), stringsAsFactors = FALSE))

predsforgraph <- fits %>%
  mutate(., p = map(fit, augment, newdata = tempsforgraphing)) %>%
  unnest(p) %>%
  #group_by(., Rate_Type) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

RateatTempgraph<-select(predsforgraph,-data,-fit)
RateatTempgraph<-as.data.frame(RateatTempgraph) 
Tempofi_graph<-RateatTempgraph$Temp

Dep_sub<-subset(RateatTempgraph, Rate_Type=="Microbe_Dependent_Rate")
Indep_sub<-subset(RateatTempgraph, Rate_Type=="Microbe_Independent_Rate")

#then find fraction: rate at Ti / rate at Topt
dependent_fractions<-data.frame()
for (i in 1:length(Dep_sub$Temp)){
  fraction_dependent<-(Dep_sub$Rate/dep_max_rate_Topt) %>% round(.,digits=3)
  ddf <- data.frame("Temp"=Tempofi_graph, fraction_dependent)
  dependent_fraction_data <- rbind(dependent_fractions,ddf) #add rates empty data frame after each calculation
}
dependent_fraction_data<-dependent_fraction_data[1:length(tempsforgraphing[,1]),]

indep_fractions<-data.frame()
for (i in 1:length(Indep_sub$Temp)){
  fraction_indep<-(Indep_sub$Rate/indep_max_rate_Topt) %>% round(.,digits=3)
  indep_df <- data.frame("Temp"=Tempofi_graph, fraction_indep)
  indep_fraction_data <- rbind(indep_fractions,indep_df) #add rates empty data frame after each calculation
}
indep_fraction_data<-indep_fraction_data[1:length(tempsforgraphing[,1]),]

fraction_data<-cbind(indep_fraction_data, dependent_fraction_data[,2]) 
names(fraction_data)[2]<-"Microbe Independent Fraction"
names(fraction_data)[3]<-"Microbe Dependent Fraction"
fraction_data #view to verify - fractions near optimum temp should be closer to 1
###########################################################
#functions needed:
#Ratkowsky for B. cereus growth rate (substitute for BtU growth rate)
Ratkowsky <- function(temp, Tmin, Tmax, b, c){
  rate <- (b*(temp-Tmin)*(1-exp(c*(temp-Tmax))))^2
  return(rate)
}

#within host model
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
  TMI <- params["TMI"]
  TMD <- params["TMD"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  dIdt <- (TMI*gamma*I*(KI-I)) + (TMD*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z))
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}
##################################################################
#now, use a loop to run the model at the desired temps using T parameters calculated above:
output.df <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!

for (i in 1:length(fraction_data$Temp)){
  times <- seq(from=0,to=3/365,by=1/365/12) #original times
  xstart<-c(H=1,I=1000,B=10000)
  #select fractions as parameters for each temp
  TMI <- fraction_data[i,2]
  TMD <- fraction_data[i,3]
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
            gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
            r=6000,d=0.003,c=0.005,sigma=500000,TMI=TMI,TMD=TMD,TB=TB) #original
  #run model at each temp w/ parameters
  ode(
    func=Bt_Tcast_within_host_infection,
    y=xstart,
    times=times,
    parms=parms
  ) %>% 
    as.data.frame() -> out 
  
  out<-mutate(out,temp=fraction_data[i,1]) #add temp to HIB vs. time output
  output.df<-rbind(output.df,out) #store in data frame
}

#next, organize data
output.df<-output.df %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(output.df)
as.data.frame(output.df)
##################
#save data frame
write.xlsx(output.df, file = "HIB_within_host_model_output_longversion.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)



#####################################

#T params
paramdf<-data.frame()
for (i in 1:length(fraction_data$Temp)){
  #select fractions as parameters for each temp
  TMI <- fraction_data[i,2]
  TMD <- fraction_data[i,3]
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  paramcols<-cbind(TMI,TMD,TB)
  
  paramcols<-cbind(paramcols,temp=fraction_data[i,1]) #add temp to HIB vs. time output
  paramdf<-rbind(paramdf,paramcols) #store in data frame
}

write.xlsx(paramdf, file = "Temp_dependent_params_withinhost.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

