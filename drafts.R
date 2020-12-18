#Thermal Optima Differences
#for range of differences
############################################
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
#############################
#import the data and remove NAs
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Rate_Data") #Emily's data
# lookinga at Rate vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
Gene_data<-na.omit(Gene_data) #remove NAs first

################################
temp_shifts <- data.frame(temp=c(seq(from=1, to=8,by=1)),stringsAsFactors = FALSE)

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
  
  #TMI <- params["TMI"]
  TMD <- params["TMD"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  #dIdt <- (TMI*gamma*I*(KI-I)) + (TMD*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z))
  dIdt <- (TMD*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z)) #new changes: get rid of microbe-independent rate
  #dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}


appendtimes<-data.frame() #must run this first
HIB_at_times_alltemps<-data.frame()
outputdataframetotals<-data.frame()
for (i in 1:length(temp_shifts$temp)){
  shifted_gene_temp_data <- Gene_data[1]+i 
  shifted_gene_data<-cbind(shifted_gene_temp_data,Gene_data[2:5])

  fits <- shifted_gene_data %>% #Use "shifted_gene_data" if wanting to model shifted data
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
  
  dep_Topt #check temps
  indep_Topt 
  
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
  for (j in 1:length(Dep_sub$Temp)){
    fraction_dependent<-(Dep_sub$Rate/dep_max_rate_Topt) %>% round(.,digits=3)
    ddf <- data.frame("Temp"=Tempofi_graph, fraction_dependent)
    dependent_fraction_data <- rbind(dependent_fractions,ddf) #add rates empty data frame after each calculation
  }
  dependent_fraction_data<-dependent_fraction_data[1:length(tempsforgraphing[,1]),]
  
  
  indep_fractions<-data.frame()
  for (k in 1:length(Indep_sub$Temp)){
    fraction_indep<-(Indep_sub$Rate/indep_max_rate_Topt) %>% round(.,digits=3)
    indep_df <- data.frame("Temp"=Tempofi_graph, fraction_indep)
    indep_fraction_data <- rbind(indep_fractions,indep_df) #add rates empty data frame after each calculation
  }
  indep_fraction_data<-indep_fraction_data[1:length(tempsforgraphing[,1]),]
  
  fraction_data<-cbind(indep_fraction_data, dependent_fraction_data[,2]) 
  
  names(fraction_data)[2]<-"Microbe Independent Fraction"
  names(fraction_data)[3]<-"Microbe Dependent Fraction"
  
  ####
  

  output.df <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!
  for (l in 1:length(fraction_data$Temp)){
    times <- seq(from=0,to=3/365,by=1/365/4) #original times
    xstart<-c(H=1,I=1000,B=10000)
    
    #select fractions as parameters for each temp
    TMD <- fraction_data[l,3]
    
    Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
    Btrate <- Ratkowsky(temp=fraction_data[l,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
    TB <- (Btrate/Btoptrate)
    
    #---> new parms w/o TMI & increase alpha:
    parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.00001,
              gamma=1,KI=20000,alpha=0.0085,KB=2000000,w=0.0005,z=0.001,
              r=6000,d=0.003,c=0.05,sigma=500000,TMD=TMD,TB=TB)
    
    #run model at each temp w/ parameters
    ode(
      func=Bt_Tcast_within_host_infection,
      y=xstart,
      times=times,
      parms=parms
    ) %>% 
      as.data.frame() -> out 
    
    out<-mutate(out,temp=fraction_data[l,1]) #add temp to HIB vs. time output
    output.df<-rbind(output.df,out) #store in data frame
  }
    #next, organize data
  output.df<-output.df %>%
    gather(variable,value,-temp,-time) %>% group_by(temp)
  head(output.df)
  output.df<-as.data.frame(output.df)
  output.df<-cbind(output.df,toptshift=i)
  #store data in bigger frame:
  outputdataframetotals<-rbind(outputdataframetotals,output.df)
  
  #####
  Health <- subset(output.df, variable=="H")
  time_values<-data.frame()
  for (m in unique(Health$temp)){
    Healthsub<-subset(Health,temp==m)
    health50percent<-which(abs(Healthsub$value-0.5)==min(abs(Healthsub$value-0.5)))
    timetoH50<-Healthsub$time[health50percent]
    times<-cbind(timetoH50,temp=m)
    time_values<-rbind(time_values,times)
  }
  time_values$daystoH50<-(time_values$timetoH50)/0.0027397 #convert time from fraction of years to days (24 hrs = 0.0027397 yrs)
  time_values<-cbind(toptshift=i,time_values)
  appendtimes<-rbind(time_values,appendtimes)
  
  HIB_at_times<-data.frame()
  #fullHIBrange<-data.frame()
  for (n in unique(output.df$temp)){
    Health <- subset(output.df, variable=="H")
    Immunity <- subset(output.df, variable=="I")
    Bacteria <- subset(output.df, variable=="B")
    #HealthImmunityBacteria<-cbind(Health,Immunity,Bacteria,temp=n,toptshift=i
    
    
    H <-which(abs(Health$value-0.001369863)==min(abs(Health$value-0.001369863)))
    I <-which(abs(Immunity$value-0.001369863)==min(abs(Immunity$value-0.001369863)))
    B <-which(abs(Bacteria$value-0.001369863)==min(abs(Bacteria$value-0.001369863)))
    times12<-cbind(H,I,B,time=12,temp=n,toptshift=i)
    #time_12hrs<-rbind(time_12hrs,times12)
    
    H <-which(abs(Health$value-0.002739726)==min(abs(Health$value-0.002739726)))
    I <-which(abs(Immunity$value-0.002739726)==min(abs(Immunity$value-0.002739726)))
    B <-which(abs(Bacteria$value-0.002739726)==min(abs(Bacteria$value-0.002739726)))
    times24<-cbind(H,I,B,time=24,temp=n,toptshift=i)
    #time_24hrs<-rbind(time_24hrs,times24)
    
    H <-which(abs(Health$value-0.005479452)==min(abs(Health$value-0.005479452)))
    I <-which(abs(Immunity$value-0.005479452)==min(abs(Immunity$value-0.005479452)))
    B <-which(abs(Bacteria$value-0.005479452)==min(abs(Bacteria$value-0.005479452)))
    times48<-cbind(H,I,B,time=48,temp=n,toptshift=i)
    #time_48hrs<-rbind(time_48hrs,times48)
    
    #merge all
    merged<-rbind(times12,times24, times48)
    HIB_at_times<-rbind(HIB_at_times,merged)
  }
  HIB_at_times_alltemps<-rbind(HIB_at_times_alltemps,HIB_at_times)
}

#head(outputdataframetotals)

#Modify data table to work with:
appendtimes<-select(appendtimes,-2)
appendtimes$toptshift <- (appendtimes$toptshift + 25.2)
appendtimes<-rename(appendtimes, Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)
head(appendtimes) #check

appendtimes %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=daystoH50),group_by(temp))+
  geom_point()+geom_smooth(se=FALSE)+
  facet_wrap(~temp)+ 
  theme_bw()+ guides(legend)+
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Days to 50% Health",x="Optimal Host Temperature (°C)",
       subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C")


appendtimes %>%
  ggplot(aes(x=temp,y=daystoH50),group_by(Optimal_Temp_of_Host_Immune_Gene_Expression))+
  geom_point()+geom_smooth(se=FALSE)+
  facet_wrap(~Optimal_Temp_of_Host_Immune_Gene_Expression)+ 
  theme_bw()+ guides(legend)+
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Days to 50% Health",x="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C")
       #subtitle="Environmental Temperature Ranges from 24-34°C")

############
Toptdiff<-cbind(appendtimes,Pathogen_Host_Topt_Diff=36-appendtimes$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(Toptdiff)

Toptdiff %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=daystoH50),group_by(temp))+
  geom_point()+geom_smooth(se=FALSE)+
  facet_wrap(~temp)+ 
  theme_bw()+
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Days to 50% Health",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C")
Toptdiff %>%
  ggplot(aes(x=temp,y=daystoH50),group_by(Pathogen_Host_Topt_Diff))+
  geom_point()+geom_smooth(se=FALSE)+
  facet_wrap(~Pathogen_Host_Topt_Diff)+ 
  theme_bw()+ 
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Days to 50% Health",x="Environmental Temperature (°C)",
       subtitle="Thermal Optimum Difference (Pathogen-Host, °C)",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C")
#subtitle="Environmental Temperature Ranges from 24-34°C")
###########################################
head(outputdataframetotals)

outputdataframetotals$toptshift <- (outputdataframetotals$toptshift + 25.2)
outputdataframetotals<-rename(outputdataframetotals, Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)

outputdataframetotals<-cbind(outputdataframetotals,Pathogen_Host_Topt_Diff=36-outputdataframetotals$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(outputdataframetotals) 

all_H <- subset(outputdataframetotals, variable=="H", select = -variable)
all_I <- subset(outputdataframetotals, variable=="I", select = -variable)
all_B <- subset(outputdataframetotals, variable=="B", select = -variable)

all_H %>% 
  ggplot(aes(x=time,y=value,color=temp,group=temp))+
  geom_point()+
  facet_wrap(~Pathogen_Host_Topt_Diff)+ geom_line()+geom_smooth(se=FALSE)+
  theme_bw() + 
  scale_x_continuous(name="Time (days)", 
                   breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                            0.006849315,0.008219178),
                   labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                   limits=c(0,0.009))+
  labs(color="Environmental Temp (°C)",y="Health",
       title="Host-Pathogen Thermal Optima Difference vs. Health",
       subtitle="Thermal Optima Difference (Pathogen minus Host, °C)",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C")+
  theme(legend.position = "bottom",legend.key.size = unit(0.4,"cm"))
  
all_I %>% 
  ggplot(aes(x=time,y=value,color=temp,group=temp))+
  geom_point()+
  facet_wrap(~Pathogen_Host_Topt_Diff)+ geom_line()+geom_smooth(se=FALSE)+
  theme_bw() + 
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))+
  labs(color="Environmental Temp (°C)",y="Immunity",
       title="Host-Pathogen Thermal Optima Difference vs. Immunity",
       subtitle="Thermal Optima Difference (Pathogen minus Host, °C)",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C")+
  theme(legend.position = "bottom",legend.key.size = unit(0.4,"cm"))

all_B %>% 
  ggplot(aes(x=time,y=value,color=temp,group=temp))+
  geom_point()+
  facet_wrap(~Pathogen_Host_Topt_Diff)+ geom_line()+geom_smooth(se=FALSE)+
  theme_bw() + 
  scale_x_continuous(name="Time (days)", 
                     breaks=c(0.000,0.001369863,0.002739726,0.004109589,0.005479452,
                              0.006849315,0.008219178),
                     labels = c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),
                     limits=c(0,0.009))+
  labs(color="Environmental Temp (°C)",y="Bacteria",
       title="Host-Pathogen Thermal Optima Difference vs. Bacteria",
       subtitle="Thermal Optima Difference (Pathogen minus Host, °C)",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C")+
  theme(legend.position = "bottom",legend.key.size = unit(0.4,"cm"))

#######################################################################
#######################################################################

#using w/in host range of Topt differences for Btwn-host model:

appendtimes<-data.frame() #must run this first
outputdataframetotals.long<-data.frame()
paramdfall<-data.frame()
for (i in 1:length(temp_shifts$temp)){
  shifted_gene_temp_data <- Gene_data[1]+i 
  shifted_gene_data<-cbind(shifted_gene_temp_data,Gene_data[2:5])
  
  fits <- shifted_gene_data %>% #Use "shifted_gene_data" if wanting to model shifted data
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
  
  dep_Topt #check temps
  indep_Topt 
  
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
  for (j in 1:length(Dep_sub$Temp)){
    fraction_dependent<-(Dep_sub$Rate/dep_max_rate_Topt) %>% round(.,digits=3)
    ddf <- data.frame("Temp"=Tempofi_graph, fraction_dependent)
    dependent_fraction_data <- rbind(dependent_fractions,ddf) #add rates empty data frame after each calculation
  }
  dependent_fraction_data<-dependent_fraction_data[1:length(tempsforgraphing[,1]),]
  
  
  indep_fractions<-data.frame()
  for (k in 1:length(Indep_sub$Temp)){
    fraction_indep<-(Indep_sub$Rate/indep_max_rate_Topt) %>% round(.,digits=3)
    indep_df <- data.frame("Temp"=Tempofi_graph, fraction_indep)
    indep_fraction_data <- rbind(indep_fractions,indep_df) #add rates empty data frame after each calculation
  }
  indep_fraction_data<-indep_fraction_data[1:length(tempsforgraphing[,1]),]
  
  fraction_data<-cbind(indep_fraction_data, dependent_fraction_data[,2]) 
  
  names(fraction_data)[2]<-"Microbe Independent Fraction"
  names(fraction_data)[3]<-"Microbe Dependent Fraction"
  
  
  
  #####
  #####
  paramdf<-data.frame()
  for (q in 1:length(fraction_data$Temp)){
    #select fractions as parameters for each temp
    TMI <- fraction_data[q,2]
    TMD <- fraction_data[q,3]
    Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14) #refer to model from Auger et al. 2008
    Btrate <- Ratkowsky(temp=fraction_data[q,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
    TB <- (Btrate/Btoptrate)
    paramcols<-cbind(TMI,TMD,TB)
    
    paramcols<-cbind(paramcols,temp=fraction_data[q,1]) #add temp to HIB vs. time output
    paramdf<-rbind(paramdf,paramcols) #store in data frame
  }
  paramdf<-cbind(paramdf,toptshift=i,Hosttopt=25.2+i,Toptdiff=36-(25+i))

  paramdfall<-rbind(paramdfall,paramdf)
  
  #####
  #####
  
  
  
  ####
  
  
  output.df.long <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!
  for (l in 1:length(fraction_data$Temp)){
    times <- seq(from=0,to=3/365,by=1/365/12) #divide by 12 for smaller time intervals
    xstart<-c(H=1,I=1000,B=10000)
    
    #select fractions as parameters for each temp
    TMD <- fraction_data[l,3]
    
    Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
    Btrate <- Ratkowsky(temp=fraction_data[l,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
    TB <- (Btrate/Btoptrate)
    
    #---> new parms w/o TMI & increase alpha:
    parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.00001,
              gamma=1,KI=20000,alpha=0.0085,KB=2000000,w=0.0005,z=0.001,
              r=6000,d=0.003,c=0.05,sigma=500000,TMD=TMD,TB=TB)
    
    #run model at each temp w/ parameters
    ode(
      func=Bt_Tcast_within_host_infection,
      y=xstart,
      times=times,
      parms=parms
    ) %>% 
      as.data.frame() -> out 
    
    out<-mutate(out,temp=fraction_data[l,1]) #add temp to HIB vs. time output
    output.df.long<-rbind(output.df.long,out) #store in data frame
  }
  #next, organize data
  output.df.long<-output.df.long %>%
    gather(variable,value,-temp,-time) %>% group_by(temp)
  #head(output.df.long)
  output.df.long<-as.data.frame(output.df.long)
  output.df.long<-cbind(output.df.long,toptshift=i)
  #store data in bigger frame:
  outputdataframetotals.long<-rbind(outputdataframetotals.long,output.df.long)

}

outputdataframetotals.long

outputdataframetotals.long$toptshift <- (outputdataframetotals.long$toptshift + 25.2)
outputdataframetotals.long<-rename(outputdataframetotals.long, Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)

outputdataframetotals<-cbind(outputdataframetotals.long,Pathogen_Host_Topt_Diff=36-outputdataframetotals.long$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(outputdataframetotals.long) 


#save data frames
write.xlsx(outputdataframetotals.long, file = "HIB_within_host_model_output_with_shifts_long.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

head(paramdfall)
write.xlsx(paramdfall, file = "Temp_dependent_params_withinhost_Decchanges_with_Topt_shifts.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

#########################################################################
###################################################################
###########################################################################
###########################################################################
###########################################################################

##################### SECTION 2: BETWEEN HOST MODEL #######################

# clear workspace
rm(list = ls())
#####################
#import the data
withinhostmodeldata <- read.xlsx("HIB_within_host_model_output_with_shifts_long.xlsx") #import within host model HIB data
Health <- subset(withinhostmodeldata, variable=="H") #subset to only get health variable

tempparamdata<-read.xlsx("Temp_dependent_params_withinhost_Decchanges_with_Topt_shifts.xlsx") #import Temp param data

Repro_data <- read_xlsx("Tcast_reproduction_Park1948data.xlsx") #reproduction v. temp data (Park 1948)
###########################################


##################### SECTION 2: BETWEEN HOST MODEL ########################

#import the data
withinhostmodeldata <- read.xlsx("HIB_within_host_model_output_longversion.xlsx") #import within host model HIB data
Health <- subset(withinhostmodeldata, variable=="H") #subset to only get health variable

tempparamdata<-read.xlsx("Temp_dependent_params_withinhost_Decchanges.xlsx") #import Temp param data

Repro_data <- read_xlsx("Tcast_reproduction_Park1948data.xlsx") #reproduction v. temp data (Park 1948)
###########################################
#extract time values for time to 50% health from data frame:
time_values<-data.frame()
for (i in unique(Health$temp)){
  Healthsub<-subset(Health,temp==i)
  health50percent<-which(abs(Healthsub$value-0.5)==min(abs(Healthsub$value-0.5)))
  timetoH50<-Healthsub$time[health50percent]
  times<-cbind(timetoH50,temp=i)
  time_values<-rbind(time_values,times)
}
time_values$daystoH50<-(time_values$timetoH50)/0.0027397 #convert time from fraction of years to days (24 hrs = 0.0027397 yrs)

#combine temp, TB, and time to 50% health into one data frame for use:
paramtable<-as.data.frame(cbind(Time50H=time_values$daystoH50,Temp=time_values$temp,TB=tempparamdata$TB))

savingparams<-as.data.frame(cbind(Temp=time_values$temp, t50_value=time_values$daystoH50,TB=tempparamdata$TB ))
write.xlsx(savingparams, file = "Temp_dependent_params_betweenhost.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

####################################################3
#incorporate Temp-dependent parameter for T. cast reproduction rate (Data from Park 1948)
Repro_gaussian <- function(temp, rmax, topt, a){
  reprorate <- rmax*exp(-0.5*(abs(temp-topt)/a)^2)
  return(reprorate)
}

gausfit<-nls_multstart(Mean_eggs~gaussian_1987(temp = Temp, rmax, topt, a),
                       data = Repro_data,
                       iter = c(4,4,4),
                       start_lower = get_start_vals(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987') - 10,
                       start_upper = get_start_vals(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987') + 10,
                       lower = get_lower_lims(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987'),
                       upper = get_upper_lims(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987'),
                       supp_errors = 'Y',
                       convergence_count = FALSE)
summary(gausfit)
#parms<-calc_params(gausfit)
coeff_list<-coeffs(gausfit)

Repro_rmax<-coeff_list[[1]]
Repro_topt<-coeff_list[[2]]
Repro_a<-coeff_list[[3]]

#################################################
#3 compartment between host model:
SIS_between_host_model_with_compartmentP<- function (t, x, params) {
  #state variables
  S <- x[1] #Susceptible
  I <- x[2] #Infected
  P <- x[3] #Dead and spore producing
  state <-c(S,I,P)
  
  #parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  d <- params["d"]
  delta <- params["delta"]
  r<-params["r"]
  K<-params["K"]
  g<-params["g"]
  mu<-params["mu"]
  TB<-params["TB"] #bacterial growth rate fraction (from within host / Ratkowsky models)
  TR<-params["TR"]
  
  #model equations
  #better ODEs for model:
  dSdt <- TR*r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S) #susceptible beetles
  dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I) #infected beetles
  dPdt <- (g*TB*(d+delta)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  
  dndt <- c(dSdt,dIdt,dPdt)
  list(dndt) #must be list format
}
#######################################################
#create initial empty data frame to store model output
SIP_output.df <- data.frame()#must do this every time prior to running new models!
#run model at diff temps and store output:
gammaval_output <-data.frame()
for (i in 1:length(paramtable$Temp)){
  t50 <- paramtable[i,1]
  TB <- paramtable[i,3]
  
  Repro_rate_at_T<-Repro_gaussian(temp=paramtable[i,2],rmax=Repro_rmax,topt=Repro_topt,a=Repro_a)
  TR <- (Repro_rate_at_T/Repro_rmax)
  
  xstart<-c(S=100,I=100,P=0) #start with 100 individuals in each compartment
  times<-seq(0,300,by=1) #300 days, by 1 day
  
  parms <- c (beta = 0.008, #Beta = contact rate * transmission probability
              #Beta = % of cases from overall pop that result in infection
              gamma = 1/t50, #gamma = 1/(infectious period)
              d=0.0005,#natural death rate
              delta=0.0005,#additional death rate due to infection
              r=0.8,#0.8#reproductive rate of susceptible beetles
              K=20000,#carrying capacity of live beetles
              g=100, #sporulation rate of dead infected beetle #this has huge impact on shape
              #Kp=1000, #carrying capacity of spores per infected beetle
              mu=0.05, #decay of dead beetles producing spores
              TB = TB,#fraction of Topt bacterial growth
              TR = TR) #fraction of Topt T. castaneum reproduction rate
  
  #run model at each temp w/ parameters
  ode(
    func=SIS_between_host_model_with_compartmentP,
    y=xstart,
    times=times,
    parms=parms
  ) %>% 
    as.data.frame() -> SIP_out 
  
  SIP_out<-mutate(SIP_out,temp=paramtable[i,2]) #add temp to output
  SIP_output.df<-rbind(SIP_output.df,SIP_out) #store in data frame
  
  #save gamma values as well for reference
  gammavalue<-(1/t50)
  gammaout<-data.frame(gammavalue)
  gammaout<-cbind(gammaout, temp=paramtable[i,2]) #add temp
  gammaval_output <- rbind(gammaval_output,gammaout)
}

#next, organize data
SIP_output.df<-SIP_output.df %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(SIP_output.df)
as.data.frame(SIP_output.df)

##################
#save data frame to excel spreadsheet in working directory:
write.xlsx(SIP_output.df, file = "SIP_3compartment_between_host_model_output.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

write.xlsx(gammaval_output, file="Temp_dependent_gamma_values_between_host_model.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

SIP_subset<-subset(SIP_output.df,variable=="S"|variable=="I")
########################
#Graph output at diff temps & save graphs to working directory:
modelgraph<-ggplot(SIP_output.df, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~(variable),scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")
modelgraph

#graph SI only:
modelgraphSIonly<-ggplot(SIP_subset, aes(x=time,y=(value),colour = temp, group=temp))+
  geom_line()+
  facet_wrap(~(variable),scales="free_y",labeller = labeller(.multi_line = FALSE))+
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")
modelgraphSIonly

jpeg(file="SIP3C_24to34.jpeg",width=600,height=500)
modelgraph
dev.off()

model_with_colors<-SIP_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(aes(color=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(size="legend",shape="legend")+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_colour_discrete(name="Temp(°C)")
model_with_colors

jpeg(file="SIP3C_24to34_withcolors.jpeg",width=600,height=500)
model_with_colors
dev.off()

SIP_output.df %>%
  ggplot(aes(x=time,y=value,color=variable))+geom_line(aes(linetype=(as.factor(temp))))+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_linetype_discrete(name="Temperature(°C)")

