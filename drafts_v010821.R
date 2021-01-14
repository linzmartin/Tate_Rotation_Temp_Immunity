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
library(deSolve)
library(xlsx)
library(DescTools)
#library(reshape2)
#library(metR)
#############################
#import the data and remove NAs
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Rate_Data") #Emily's data
# looking at Rate vs. temp 
# this is condensed data from Emily's work - find in my box folder
Gene_data<-na.omit(Gene_data) #remove NAs first

################################

#create list of numbers to add to thermal optimum of host immunity to create Topt shift
temp_shifts <- data.frame(temp=c(seq(from=1, to=8,by=1)),stringsAsFactors = FALSE)

#####
#store functions that will be needed:

#Ratkowsky bacterial growth model
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
  
  TMD <- params["TMD"]
  TB <- params["TB"]
  
  #model equations
  dHdt <- (psi*H*(KH-H))-(beta*B*H)-(p*H)-(I*m*H)
  dIdt <- (TMD*alpha*I*(KB/(1+exp(-0.001*(B-sigma)))))-(I*(B*w + z)) #new changes: get rid of microbe-independent rate & TMI parameter
  dBdt <- TB*r*B*(1-B/KB)-B*d-c*I*B
  
  dndt <- c(dHdt,dIdt,dBdt)
  list(dndt) #must be list format for ode
}

###################
#run the Within Host Model with different shifts in Topt 

appendtimes<-data.frame() #must run these empty frames first to initialize
outputdataframetotals<-data.frame()
HIB_at_times_alltemps<-data.frame()
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
  #####################################################
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
  #####################################################
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
    
    #---> new parms w/o TMI & increase alpha to 0.0085:
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
  
  ###
  HIB_at_times<-data.frame()
  for (n in unique(output.df$temp)){
    Health <- subset(output.df, variable=="H")
    Immunity <- subset(output.df, variable=="I")
    Bacteria <- subset(output.df, variable=="B")
    
    H <-Health[which((abs(Health$time-0.001369863)==min(abs(Health$time-0.001369863))) & Health$temp==n),]
    I <-Immunity[which((abs(Immunity$time-0.001369863)==min(abs(Immunity$time-0.001369863))) & Health$temp==n),]
    B <-Bacteria[which((abs(Bacteria$time-0.001369863)==min(abs(Bacteria$time-0.001369863))) & Health$temp==n),]
    H_val<-H$value
    I_val<-I$value
    B_val<-B$value
    times12<-cbind(H_val,I_val,B_val,time=12,temp=n,toptshift=i)
    
    H <-Health[which((abs(Health$time-0.002739726)==min(abs(Health$time-0.002739726))) & Health$temp==n),]
    I <-Immunity[which((abs(Immunity$time-0.002739726)==min(abs(Immunity$time-0.002739726))) & Health$temp==n),]
    B <-Bacteria[which((abs(Bacteria$time-0.002739726)==min(abs(Bacteria$time-0.002739726))) & Health$temp==n),]
    H_val<-H$value
    I_val<-I$value
    B_val<-B$value
    times24<-cbind(H_val,I_val,B_val,time=24,temp=n,toptshift=i)
    
    H <-Health[which((abs(Health$time-0.005479452)==min(abs(Health$time-0.005479452))) & Health$temp==n),]
    I <-Immunity[which((abs(Immunity$time-0.005479452)==min(abs(Immunity$time-0.005479452)))& Health$temp==n),]
    B <-Bacteria[which((abs(Bacteria$time-0.005479452)==min(abs(Bacteria$time-0.005479452)))& Health$temp==n),]
    H_val<-H$value
    I_val<-I$value
    B_val<-B$value
    times48<-cbind(H_val,I_val,B_val,time=48,temp=n,toptshift=i)
    
    #merge all to one data frame
    merged<-rbind(times12,times24, times48)
    HIB_at_times<-rbind(HIB_at_times,merged)
  }
  HIB_at_times_alltemps<-rbind(HIB_at_times_alltemps,HIB_at_times)
  
}

#Modify data table to work with it more easily:
appendtimes<-select(appendtimes,-2)
appendtimes$toptshift <- (appendtimes$toptshift + 25.2)
appendtimes<-rename(appendtimes, Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)
head(appendtimes) #check

head(outputdataframetotals)

###
###############
#Plot the outputs:
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

#### contour plot:
#https://jkzorz.github.io/2020/02/29/contour-plots.html

appendtimes %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=temp))+
  geom_raster(aes(fill=daystoH50))+
  #geom_contour(aes(z=daystoH50),colour="white")+
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Environmental Temperature (°C)",x="Optimal Host Temperature (°C)",
       subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C",
       fill="Time to 50% Health (days)")+
  theme(legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "top", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11), legend.key = element_blank()) + 
  scale_fill_continuous(low = "#BFE1B0", high = "#137177") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) 
####


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

#contour plot:
Toptdiff %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=daystoH50))+
  #geom_contour(aes(z=daystoH50),colour="white")+
  labs(title="Effect of Topt Increase on Time to 50% Health",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Time to 50% Health (days)")+
  theme(legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "top", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11), legend.key = element_blank()) + 
  scale_fill_continuous(low = "#BFE1B0", high = "#137177") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

###########################################
#separate data by variable
head(outputdataframetotals)

outputdataframetotals$toptshift <- (outputdataframetotals$toptshift + 25.2)
outputdataframetotals<-rename(outputdataframetotals, Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)

outputdataframetotals<-cbind(outputdataframetotals,Pathogen_Host_Topt_Diff=36-outputdataframetotals$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(outputdataframetotals) 

all_H <- subset(outputdataframetotals, variable=="H", select = -variable)
all_I <- subset(outputdataframetotals, variable=="I", select = -variable)
all_B <- subset(outputdataframetotals, variable=="B", select = -variable)
########
#Host health:
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

##################################
#Health:

#contour plot:
#health at all times (not super informative)
all_H %>% 
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=value)) +#uses all Health values at all times
  labs(title="Effect of Topt Increase on Health",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C, Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Health [at all times]")+
  theme(legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "top", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11), legend.key = element_blank()) + 
  scale_fill_continuous(low = "#BFE1B0", high = "#137177") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
##
#select times to display health values - this is more informative

#format data
head(HIB_at_times_alltemps)

HIB_at_times_alltemps$toptshift <- HIB_at_times_alltemps$toptshift+25.2
HIB_at_times_alltemps <- rename(HIB_at_times_alltemps,Optimal_Temp_of_Host_Immune_Gene_Expression=toptshift)
HIB_at_times_alltemps <- cbind(HIB_at_times_alltemps,Pathogen_Host_Topt_Diff=36-HIB_at_times_alltemps$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(HIB_at_times_alltemps)
##
Health_for_btwnhost <- select(HIB_at_times_alltemps,-2,-3)
Health_for_btwnhost<-subset(Health_for_btwnhost,time==24)
head(Health_for_btwnhost)

###
Health_for_btwnhost %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=temp))+
  geom_raster(aes(fill=H_val)) +#uses all Health values at all times
  labs(title="Effect of Increasing Host Thermal Optimum on Host Health",
       y="Environmental Temperature (°C)",x="Optimal Temperature for Host Immune Gene Expression (°C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \n Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Health at 24hrs")


Health_subset <- select(HIB_at_times_alltemps,-2,-3)
Health_48hrs<-subset(Health_subset,time==48)
head(Health_48hrs)

Health_48hrs %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=temp))+
  geom_raster(aes(fill=H_val)) +#uses all Health values at all times
  labs(title="Effect of Increasing Host Thermal Optimum on Host Health",
       y="Environmental Temperature (°C)",x="Optimal Temperature for Host Immune Gene Expression (°C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \n Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Health at 48hrs")

#Thermal optima differences - contour plot for health at 24 and 48 hours post infection
Health24hrs_toptdiff<-Health_for_btwnhost

Health24hrs_toptdiff %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=H_val)) +#uses all Health values at all times
  labs(title="Host-Pathogen Thermal Optima Differences vs. Host Health",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \n Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Health at 24hrs")



Health_48hrs %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=H_val)) +#uses all Health values at all times
  labs(title="Host-Pathogen Thermal Optima Differences vs. Host Health",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \n Host Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Health at 48hrs")

##################################
#change in health over time = health at time 0 - health at time 48

head(HIB_at_times_alltemps)
deltavalues<-data.frame()
deltavalues_alltoptshifts<-data.frame()
for (w in unique(HIB_at_times_alltemps$Optimal_Temp_of_Host_Immune_Gene_Expression)){
  #health at time zero = 1 for all cases in model
  deltavalues<-data.frame()
  for (q in unique(HIB_at_times_alltemps$temp)){
   Hat24<-subset(HIB_at_times_alltemps, time==24 & Optimal_Temp_of_Host_Immune_Gene_Expression==w & temp==q)
   H24val<-Hat24$H_val
   deltavalues<-cbind(deltavalue=1-H24val,temp=q,Optimal_Temp_of_Host_Immune_Gene_Expression=w)
   deltavalues_alltoptshifts<-rbind(deltavalues_alltoptshifts,deltavalues)
  }
}
#view(deltavalues_alltoptshifts) #check the output

deltavalues_alltoptshifts<-rename(deltavalues_alltoptshifts,Health_Decline_to_24hrs=deltavalue)
head(deltavalues_alltoptshifts)

deltavalues_alltoptshifts %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=temp))+
  geom_raster(aes(fill=Health_Decline_to_24hrs)) +#Health at time zero (1) - Health at 24hrs post infection = Health decline to 24hrs
  labs(title="Host-Pathogen Thermal Optima Differences vs. \n Change in Host Health over Time",
       y="Environmental Temperature (°C)",x="Optimal Temperature for Host Immune Gene Expression (°C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Decline in Health from \n0 to 24hrs post infection")

### topt diff:
deltavalues_alltoptshifts<-cbind(deltavalues_alltoptshifts,Pathogen_Host_Topt_Diff=36-deltavalues_alltoptshifts$Optimal_Temp_of_Host_Immune_Gene_Expression)
head(deltavalues_alltoptshifts)

deltavalues_alltoptshifts %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=Health_Decline_to_24hrs)) +#Health at time zero (1) - Health at 24hrs post infection = Health decline to 24hrs
  labs(title="Host-Pathogen Thermal Optima Differences vs. \n Change in Host Health over Time",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Decline in Health from \n0 to 24hrs post infection")


####################################
#Immunity:
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


Immunity_subset <- select(HIB_at_times_alltemps,-1,-3)
head(Immunity_subset)

Immunity_time12 <- subset(Immunity_subset, time==12)
Immunity_time12 %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=I_val)) +
  labs(title="Host-Pathogen Thermal Optima Differences vs. Host Immunity",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Host Immunity \nat 12hrs")

Immunity_time24 <- subset(Immunity_subset, time==24)
Immunity_time24 %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=I_val)) +
  labs(title="Host-Pathogen Thermal Optima Differences vs. Host Immunity",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Host Immunity \nat 24hrs")

Immunity_time48 <- subset(Immunity_subset, time==48)
Immunity_time48 %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=I_val)) +
  labs(title="Host-Pathogen Thermal Optima Differences vs. Host Immunity",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Host Immunity \nat 48hrs")

####################################
#Bacteria:
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

Bacteria_subset <- select(HIB_at_times_alltemps,-1,-2) 
head(Bacteria_subset)
#view(Bacteria_subset)

#want to make a contour plot of the max bacterial load at each temp
all_max_B_values <- data.frame()
for (s in unique(all_B$Optimal_Temp_of_Host_Immune_Gene_Expression)){
  max_b <- data.frame()
  B_sub_2<-data.frame()
  
  for (t in unique(all_B$temp)){
    B_sub_2 <- subset(all_B, temp==t & Optimal_Temp_of_Host_Immune_Gene_Expression==s)
    max_b <- B_sub_2[which(B_sub_2$value==max(B_sub_2$value)),]

    all_max_B_values <- rbind(all_max_B_values,max_b)
  }
}
head(all_max_B_values) #check results

all_max_B_values %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=value)) +
  labs(title="Effects of Host-Pathogen Thermal Optima \nDifferences on Max Bacterial Load",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Max Bacterial Load")

all_max_B_values %>%
  ggplot(aes(x=Optimal_Temp_of_Host_Immune_Gene_Expression,y=temp))+
  geom_raster(aes(fill=value)) +
  labs(title="Effects of Host-Pathogen Thermal Optima \nDifferences on Max Bacterial Load",
       y="Environmental Temperature (°C)",x="Optimal Temperature for Host Immune Gene Expression (°C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Maximum \nBacterial \nLoad in Host")

all_max_B_values %>%
  ggplot(aes(x=Pathogen_Host_Topt_Diff,y=temp))+
  geom_raster(aes(fill=value)) +
  labs(title="Effects of Host-Pathogen Thermal Optima \nDifferences on Max Bacterial Load",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       #subtitle="Environmental Temperature Ranges from 24-34°C",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Max Bacterial Load")

####
#the area under the curve for bacterial growth over time:
all_AUC_values <- data.frame()
for (i in unique(all_B$Pathogen_Host_Topt_Diff)){
  subset_topt <- subset(all_B, Pathogen_Host_Topt_Diff==i)
  AUC_values <- data.frame()
  for (j in unique(all_B$temp)){
    subset_temp <- subset(subset_topt, subset_topt$temp==j)
  
    aucval<-AUC(x=subset_temp$time,y=subset_temp$value)
  
    AUC_values <- cbind(aucval,Toptdiff=i,temp=j)
    AUC_values<-as.data.frame(AUC_values)
    all_AUC_values<-rbind(all_AUC_values,AUC_values)
  }
}

all_AUC_values %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=aucval))+
  labs(title="Effects of Host-Pathogen Thermal Optima\nDifferences on Bacterial Population Growth\nOver Time (Area Under Curve)",
       y="Environmental Temperature (°C)",x="Thermal Optimum Difference (Pathogen-Host, °C)",
       caption="Pathogen Thermal Optimum = 36°C \nHost Thermal Optimum Ranges from 26.2-33.2°C",
       fill="Area Under Curve")

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
head(Health_for_btwnhost) #using Health at 24hrs post infection for btwn host model
head(paramdfall)
write.xlsx(paramdfall, file = "Temp_dependent_params_withinhost_Decchanges_with_Topt_shifts.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)

#############################################################################
############################################################################
##################### SECTION 2: BETWEEN HOST MODEL ########################

#import the data
#also using data tables above from section 1 (don't clear ls!!!)
Health <- Health_for_btwnhost
Repro_data <- read_xlsx("Tcast_reproduction_Park1948data.xlsx") #reproduction v. temp data (Park 1948)
Germination_data <- read.csv("Knaysi_1964_Fig1_Graph_grabber_data.csv", stringsAsFactors = FALSE) 
###########################################

#combine temp, TB, and time to 50% health into one data frame for use:
Health_for_btwnhost <-cbind(Health_for_btwnhost,toptshift=Health_for_btwnhost$Optimal_Temp_of_Host_Immune_Gene_Expression-25.2)
paramtable<-as.data.frame(cbind(Hat24hrs=Health_for_btwnhost$H_val,Temp=Health_for_btwnhost$temp,TB=paramdfall$TB,toptshift=Health_for_btwnhost$toptshift))

savingparams<-as.data.frame(cbind(Temp=Health_for_btwnhost$temp, Hat24hrs=Health_for_btwnhost$H_val,TB=paramdfall$TB, toptshift=Health_for_btwnhost$toptshift))
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
coeff_list<-coefficients(gausfit)

Repro_rmax<-coeff_list[[1]]
Repro_topt<-coeff_list[[2]]
Repro_a<-coeff_list[[3]]

###########################################
#incorporate temp dependent parameter for spore germination TG
Ratkowsky_germination <- function(temp, Tmin, Tmax, a,b){
  rate <- (a*(temp-Tmin)*(1-exp(b*(temp-Tmax))))^2
  return(rate)
}

ratkowsky_germ_fit <-nls_multstart(Percent_Dark_Spores~ratkowsky_1983(temp = Temperature, tmin, tmax, a, b),
                                   data = Germination_data,
                                   iter = c(4,4,4,4),
                                   start_lower = get_start_vals(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'ratkowsky_1983') - 10,
                                   start_upper = get_start_vals(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'ratkowsky_1983') + 10,
                                   lower = get_lower_lims(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'ratkowsky_1983'),
                                   upper = get_upper_lims(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'ratkowsky_1983'),
                                   supp_errors = 'Y',
                                   convergence_count = FALSE)
summary(ratkowsky_germ_fit)

coeff_list<-coeffs(ratkowsky_germ_fit)
coeff_list
Germ_tmin<-coeff_list[[1]]
Germ_tmax<-coeff_list[[2]]
Germ_a<-coeff_list[[3]]
Germ_b<-coeff_list[[4]]

all_germ_params<-calc_params(ratkowsky_germ_fit)
germination_rmax<-all_germ_params[[1]]
germination_rmax

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
  Hat24<-params["Hat24"]
  TG<-params["TG"]
  
  #model equations
  #better ODEs for model:
  dSdt <- TR*r*(1-(S+I)/K)-(beta*S*P*TG)+(gamma*I)-(d*S) #susceptible beetles
  dIdt <- (beta*S*P*TG)-(gamma*I)-((d+delta*Hat24)*I) #infected beetles
  dPdt <- (g*TB*(d+delta*Hat24)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  
  dndt <- c(dSdt,dIdt,dPdt)
  list(dndt) #must be list format
}
#######################################################

SIP_outputdataframetotals<-data.frame()
SIP_out<-data.frame()
SIP_output.df <- data.frame()

for (i in 1:length(paramtable$toptshift)){
    TB <- paramtable[i,3]
    Hat24<- paramtable[i,1]
    
    Repro_rate_at_T<-Repro_gaussian(temp=paramtable[i,2],rmax=Repro_rmax,topt=Repro_topt,a=Repro_a)
    TR <- (Repro_rate_at_T/Repro_rmax)

    Germ_rate_at_T <- Ratkowsky_germination(temp=paramtable[i,2],Tmin=Germ_tmin,Tmax=Germ_tmax,a=Germ_a,b=Germ_b)
    TG <- (Germ_rate_at_T/germination_rmax)
    
    xstart<-c(S=100,I=100,P=0) #start with 100 individuals in each compartment
    times<-seq(0,300,by=1) #300 days, by 1 day
    
    parms <- c (beta = 0.008, #Beta = contact rate * transmission probability
                #Beta = % of cases from overall pop that result in infection
                gamma = 0.5, #gamma = 1/(infectious period)
                d=0.0005,#natural death rate
                delta=0.0005,#additional death rate due to infection
                r=0.8,#0.8#reproductive rate of susceptible beetles
                K=300,#carrying capacity of live beetles
                g=100, #sporulation rate of dead infected beetle #this has huge impact on shape
                #Kp=1000, #carrying capacity of spores per infected beetle
                mu=0.05, #decay of dead beetles producing spores
                TB = TB,#fraction of Topt bacterial growth
                TR = TR, #fraction of Topt T. castaneum reproduction rate
                TG = TG, #fraction of Topt Bacillus spore germination rate
                Hat24 = Hat24) #Health at 24hrs from within-host model output
    
    #run model at each temp w/ parameters
    ode(
      func=SIS_between_host_model_with_compartmentP,
      y=xstart,
      times=times,
      parms=parms
    ) %>% 
      as.data.frame() -> SIP_out 
    
    SIP_out<-cbind(SIP_out,temp=paramtable[i,2],
                    toptshift=paramtable[i,4], 
                    N=(SIP_out$S + SIP_out$I),
                    prop_I=SIP_out$I/(SIP_out$I + SIP_out$S),
                    prop_S=SIP_out$S/(SIP_out$I + SIP_out$S))
    SIP_out <- as.data.frame(SIP_out)
  
  #store data in bigger frame:
  SIP_outputdataframetotals<-rbind(SIP_outputdataframetotals,SIP_out)  
} 
SIP_outputdataframetotals<-SIP_outputdataframetotals %>%
  gather(variable,value,-temp,-time,-toptshift) %>% group_by(toptshift) %>% group_by(temp)

SIP_outputdataframetotals<-cbind(SIP_outputdataframetotals,Optimal_Temp_of_Host_Immune_Gene_Expression=(SIP_outputdataframetotals$toptshift+25.2),Toptdiff=(36-(SIP_outputdataframetotals$toptshift+25.2)))
head(SIP_outputdataframetotals)

#####
#next, subset data

SIP_outputdataframetotals_subsetSIP <- subset(SIP_outputdataframetotals, variable=="S"|variable=="I"|variable=="P")
SIP_outputdataframetotals_subset_toptexample <- subset(SIP_outputdataframetotals_subset, toptshift==1 | toptshift==5)
SIP_outputdataframetotals_subset_proportions <- subset(SIP_outputdataframetotals, variable=="prop_I"|variable=="prop_S"|variable=="N")

#plot the data

SIP_outputdataframetotals_subsetSIP %>% 
  group_by(toptshift) %>%
  ggplot(., aes(x=time,y=value,colour=temp,group=temp))+
  geom_line()+
  facet_grid(variable~toptshift)+ #scales are not ="free_y"
  scale_x_continuous(name="Time (days)")+
  labs(y="Number of Individual Beetles")

SIP_outputdataframetotals_subsetSIP %>% 
  ggplot(aes(x=time,y=value,group=temp),group_by(toptshift))+
  geom_line(aes(color=(as.factor(temp))))+
  facet_grid(variable~toptshift)+
  guides(size="legend",shape="legend")+theme_bw()+
  labs(y="Number of Individual Beetles")+
  scale_x_continuous(name="Time (days)")+
  scale_colour_discrete(name="Temp(°C)")

SIP_outputdataframetotals_subset_proportions%>% 
  group_by(toptshift) %>%
  ggplot(aes(x=time,y=value,colour=temp,group=temp))+
  geom_line()+
  facet_grid(variable~toptshift,scales="free_y")+
  labs(y="Number of Individual Beetles",x="Time (days)")

SIP_outputdataframetotals_subset_toptexample %>% 
  ggplot(aes(x=time,y=value,colour=temp,group=temp))+
  geom_line()+
  facet_grid(variable~toptshift)+
  labs(y="Number of Individual Beetles",x="Time (days)")

########################

#contour plots

#metrics: pop size (=N), 
  #disease prevalence (total number of cases at time point / total pop = prop I),
  #spore abundance (=P)
  #vs. temp or gamma
  #using metrics from equilibrium - ~200 days

  #subset data: 50 days, 100 days, 200 days

#############################
#Disease prevalence: 

PropI_50days <- subset(SIP_outputdataframetotals_subset_proportions, time==50 & variable=="prop_I")
PropI_100days <- subset(SIP_outputdataframetotals_subset_proportions, time==100 & variable=="prop_I")
PropI_200days <- subset(SIP_outputdataframetotals_subset_proportions, time==200 & variable=="prop_I")

PropI_50days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Disease Prevalence at 50 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nDisease Prevalence = proportion infected / total live population at 50 days",
       fill="Disease \nPrevalance \nat 50 days")

PropI_100days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Disease Prevalence at 100 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nDisease Prevalence = proportion infected / total live population at 50 days",
       fill="Disease \nPrevalance \nat 100 days")

PropI_200days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Disease Prevalence at 200 days (at Equilibrium)",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nDisease Prevalence = proportion infected / total live population at 50 days",
       fill="Disease \nPrevalance \nat 200 days")

##############
#Population size:
Popsize_50days <- subset(SIP_outputdataframetotals_subset_proportions, time==50 & variable=="N")
Popsize_100days <- subset(SIP_outputdataframetotals_subset_proportions, time==100 & variable=="N")
Popsize_200days <- subset(SIP_outputdataframetotals_subset_proportions, time==200 & variable=="N")

Popsize_50days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Population Size at 50 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nPopulation Size = S+I, open population",
       fill="Population Size\nat 50 days")

Popsize_100days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Population Size at 100 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nPopulation Size = S+I, open population",
       fill="Population Size\nat 100 days")

Popsize_200days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Population Size at 200 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nPopulation Size = S+I, open population",
       fill="Population Size\nat 200 days")
######################
#Spore Abundance:

Spores_50days <- subset(SIP_outputdataframetotals_subsetSIP, time==50 & variable=="P")
Spores_100days <- subset(SIP_outputdataframetotals_subsetSIP, time==100 & variable=="P")
Spores_200days <- subset(SIP_outputdataframetotals_subsetSIP, time==200 & variable=="P")

Spores_50days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Spore Abundance at 50 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nSpore Abundance = P",
       fill="Spore Abundance\nat 50 days")

Spores_100days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Spore Abundance at 100 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nSpore Abundance = P",
       fill="Spore Abundance\nat 100 days")

Spores_200days %>%
  ggplot(aes(x=Toptdiff,y=temp))+
  geom_raster(aes(fill=value))+
  labs(title="Spore Abundance at 200 days",
       x="Thermal Optimum Difference (Pathogen-Host, °C)",
       y="Environmental Temperature (°C)",
       caption="Pathogen Thermal Optimum = 36°C\nHost Thermal Optimum Ranges from 26.2-33.2°C\nSpore Abundance = P",
       fill="Spore Abundance\nat 200 days")
