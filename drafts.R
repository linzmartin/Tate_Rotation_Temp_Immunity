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
  as.data.frame(output.df)
  
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
  for (n in unique(output.df$temp)){
    Health <- subset(output.df, variable=="H")
    Immunity <- subset(output.df, variable=="I")
    Bacteria <- subset(output.df, variable=="B")
    
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



#appendtimes<-appendtimes[-timetoH50]





