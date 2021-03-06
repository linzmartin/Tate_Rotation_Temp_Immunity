#BtU Tcast within host model - LOOP through different temperatures
#graphs at end
#complete
####################
# clear workspace
rm(list = ls())
#set wd to your project folder
setwd("C:/Users/linzm/Documents/Tate_Lab_Rotation/Tate_Rotation_Temp_Immunity")
###########################
#load packages
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
#################################
#creating temp data/ model for calculating fraction of Topt
###############################################
#import the data and remove NAs
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Gene_Data_Modified") #Emily's data
# lookinga at ddCT2 (fold change) vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
Gene_data<-na.omit(Gene_data) #remove NAs first
###############################################
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


#########################
#calculate ddCT2 based on model at temps to find max ddCT and %s

##first, we need to calculate the max values:
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
#length(neededtemps[,1])
rate_data<-rate_data[1:length(neededtemps[,1]),] #now, have T 24-24 only once

#find ddCT2rate at Topt 25.2 = max rate
constit_max_rate_Topt <- subset(rate_data, Temp==25.2,select = c("constitutive_rate"))
constit_max_rate_Topt<-constit_max_rate_Topt[,1]
induced_max_rate_Topt<-subset(rate_data, Temp==25.2,select = c("induced_rate"))
induced_max_rate_Topt<-induced_max_rate_Topt[,1]
##########
#next, we can calculate predictions for every 2C (or degrees of interest):
tempsforgraphing<-Gene_data %>%
  do(., data.frame(Temp = c(seq(from=24,to=34,by=2)), stringsAsFactors = FALSE))

predsforgraph <- fits %>%
  mutate(., p = map(fit, augment, newdata = tempsforgraphing)) %>%
  unnest(p) %>%
  group_by(., Treatment) %>%
  rename(., ddCT2 = .fitted) %>%
  ungroup()

ddCT2atTempgraph<-select(predsforgraph,-data,-fit)
ddCT2atTempgraph<-as.data.frame(ddCT2atTempgraph) 

rate_data_forgraph <- data.frame() #create empty data frame
#then, calculate constitutive and induced rates 
for (i in 1:length(ddCT2atTempgraph$Temp)){
  HK <- filter(ddCT2atTempgraph, Treatment=="heat killed")
  NI <- filter(ddCT2atTempgraph, Treatment=="non-injected")
  BtU <- filter(ddCT2atTempgraph, Treatment=="BtU")
  
  constitutive_rate_graph<-(HK$ddCT2-NI$ddCT2)/4 #time frame is 4 hours
  induced_rate_graph<-(BtU$ddCT2-NI$ddCT2)/4 #time frame is 4 hours
  Tempofi_graph<-ddCT2atTempgraph$Temp
  
  data_frame_graph <- data.frame("Temp"=Tempofi_graph, constitutive_rate_graph,induced_rate_graph)
  rate_data_forgraph <- rbind(rate_data_forgraph, data_frame_graph) #add rates empty data frame after each calculation
}
#get rid of repeated data:
#length(neededtemps[,1])
rate_data_forgraph<-rate_data_forgraph[1:length(tempsforgraphing[,1]),] #now, have T 24-24 only once


#then find fraction: ddCT2 rate at Ti / rate at Topt
fractions<-data.frame()
for (i in 1:length(rate_data_forgraph$Temp)){
  fraction_induced<-(rate_data_forgraph$induced_rate/induced_max_rate_Topt) %>% round(.,digits=3)
  fraction_constitutive<-(rate_data_forgraph$constitutive_rate/constit_max_rate_Topt) %>% round(.,digits=3)
  
  df <- data.frame("Temp"=Tempofi_graph, fraction_constitutive,fraction_induced)
  fraction_data <- rbind(fractions,df) #add rates empty data frame after each calculation
}
#get rid of repeated data:
fraction_data<-fraction_data[1:length(tempsforgraphing[,1]),]
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
output.df <- data.frame() #Temp = Temp.vector, model = NA

for (i in 1:length(fraction_data$Temp)){
  times <- seq(from=0,to=3/365,by=1/365/4) #original times
  xstart<-c(H=1,I=1000,B=10000)
  
  #select fractions as parameters for each temp
  TCI <- fraction_data[i,2]
  TII <- fraction_data[i,3]
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  
  parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
            gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
            r=6000,d=0.003,c=0.005,sigma=500000,TCI=TCI,TII=TII,TB=TB) #original params
  ###
  #parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.000001,
   #        gamma=1,KI=20000,alpha=0.0085,KB=2000000,w=0.0005,z=0.001,
    #      r=6000,d=0.003,c=0.5,sigma=500000,TCI=TCI,TII=TII,TB=TB) #parms for high I mods - conservative changes
  #parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.000001,
  #         gamma=100,KI=20000,alpha=0.85,KB=2000000,w=0.0005,z=0.001,
   #        r=6000,d=0.003,c=1,sigma=500000,TCI=TCI,TII=TII,TB=TB) #params w/ high I mods - extreme changes
  #times <- seq(from=0,to=12,by=1) #for params w/ extreme high I to test model
  ###
  
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
###########################################
#then plot models on graph
modelgraph<-ggplot(output.df, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~variable,scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)", 
                    breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                            0.006849315,0.008219178),
                  labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                 limits=c(0,0.009))+ 
  labs(y=NULL)

jpeg(file="within_host_model_24to34C.jpeg",width=1000,height=500)
modelgraph
dev.off()

modelgraph #view
#####
#another type of graph
jpeg(file="within_host_model_temps_colored.jpeg",width=1200,height=800)
model_with_colors<-output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(aes(color=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(size="legend",shape="legend")+theme_bw()+
  labs(y="State Variable")+
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))+
  scale_colour_discrete(name="Temp(°C)")
model_with_colors
dev.off()

model_with_colors #view
#############################################
output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line(aes(linetype=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="State Variable")+
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))+
  scale_linetype_discrete(name="Temperature(°C)")
