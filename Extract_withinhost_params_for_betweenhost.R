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

constit_topt <- calc_params(fits$fit[[1]]) %>% select(topt)
constit_Topt <- as.numeric(constit_topt) %>% round(.,digits=1)
induced_topt <- calc_params(fits$fit[[2]]) %>% select(topt)
induced_Topt <- as.numeric(induced_topt) %>% round(.,digits=1)

#find rate at Topt = max rate
constit_max_rate_Topt <- subset(RateatTemp, Rate_Type=="Constitutive_Rate" & Temp==constit_Topt,select = c("Rate"))
constit_max_rate_Topt<-constit_max_rate_Topt[,1]
induced_max_rate_Topt<-subset(RateatTemp, Rate_Type=="Induced_Rate" & Temp==induced_Topt,select = c("Rate"))
induced_max_rate_Topt<-induced_max_rate_Topt[,1]
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

Induced_sub<-subset(RateatTempgraph, Rate_Type=="Induced_Rate")
Constitutive_sub<-subset(RateatTempgraph, Rate_Type=="Constitutive_Rate")

#then find fraction: rate at Ti / rate at Topt
induced_fractions<-data.frame()
for (i in 1:length(Induced_sub$Temp)){
  fraction_induced<-(Induced_sub$Rate/induced_max_rate_Topt) %>% round(.,digits=3)
  df <- data.frame("Temp"=Tempofi_graph, fraction_induced)
  induced_fraction_data <- rbind(induced_fractions,df) #add rates empty data frame after each calculation
}
induced_fraction_data<-induced_fraction_data[1:length(tempsforgraphing[,1]),]

constit_fractions<-data.frame()
for (i in 1:length(Constitutive_sub$Temp)){
  fraction_constitutive<-(Constitutive_sub$Rate/constit_max_rate_Topt) %>% round(.,digits=3)
  constit_df <- data.frame("Temp"=Tempofi_graph, fraction_constitutive)
  constit_fraction_data <- rbind(constit_fractions,constit_df) #add rates empty data frame after each calculation
}
constit_fraction_data<-constit_fraction_data[1:length(tempsforgraphing[,1]),]

fraction_data<-cbind(constit_fraction_data, induced_fraction_data[,2]) 
names(fraction_data)[2]<-"Constitutive fraction"
names(fraction_data)[3]<-"Induced Fraction"
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
##################################################################
#now, use a loop to run the model at the desired temps using T parameters calculated above:
output.df <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!

for (i in 1:length(fraction_data$Temp)){
  times <- seq(from=0,to=3/365,by=1/365/12) #original times
  xstart<-c(H=1,I=1000,B=10000)
  #select fractions as parameters for each temp
  TCI <- fraction_data[i,2]
  TII <- fraction_data[i,3]
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
            gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
            r=6000,d=0.003,c=0.005,sigma=500000,TCI=TCI,TII=TII,TB=TB) #original
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
  TCI <- fraction_data[i,2]
  TII <- fraction_data[i,3]
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  paramcols<-cbind(TCI,TII,TB)
  
  paramcols<-cbind(paramcols,temp=fraction_data[i,1]) #add temp to HIB vs. time output
  paramdf<-rbind(paramdf,paramcols) #store in data frame
}

write.xlsx(paramdf, file = "Temp_dependent_params_withinhost.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

