#mismatch between host and pathogen thermal optima

#host thermal optimum
#pathogen thermal optimum 
#how does increased host Topt (closer to pathogen Topt) affect model output?
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

#################### Section 1: within-host model #####################


#import the data and remove NAs
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Rate_Data") #Emily's data
# lookinga at Rate vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
Gene_data<-na.omit(Gene_data) #remove NAs first

#to shift Host thermal optimum:
shifted_gene_temp_data <- Gene_data[1]+8 #shift by 8C or temp of interest
shifted_gene_data<-cbind(shifted_gene_temp_data,Gene_data[2:5])
head(shifted_gene_data)
###################################
#Run the Flinn model by treatment group

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
head(fits) #check output of fits

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

plot(x=fraction_data$Temp,y=fraction_data$`Microbe Dependent Fraction`) #check to make sure fraction max is where topt is

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
##################################################################
#now, use a loop to run the model at the desired temps using T parameters calculated above:
output.df <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!

for (i in 1:length(fraction_data$Temp)){
  times <- seq(from=0,to=3/365,by=1/365/4) #original times
  xstart<-c(H=1,I=1000,B=10000)
  
  #select fractions as parameters for each temp
  #TMI <- fraction_data[i,2]
  TMD <- fraction_data[i,3]
  
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  
  
  #parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
  #         gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
  #        r=6000,d=0.003,c=0.005,sigma=500000,TMI=TMI,TMD=TMD,TB=TB) #original, no increase in I, with TMI
  
  #parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.00001,
  #gamma=1,KI=20000,alpha=0.00085,KB=2000000,w=0.0005,z=0.001,
  #r=6000,d=0.003,c=0.05,sigma=500000,TMI=TMI,TMD=TMD,TB=TB) #parms for high I mods - conservative changes, with TMI
 
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
  
  out<-mutate(out,temp=fraction_data[i,1]) #add temp to HIB vs. time output
  output.df<-rbind(output.df,out) #store in data frame
}

#next, organize data
output.df<-output.df %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(output.df)

as.data.frame(output.df)
##################
#save data frame (modify name of file as needed)
#write.xlsx(output.df, file = "HIB_within_host_model_output_OG+8C+c*TMD+conserv_better_I.xlsx",
#          sheetName = "1", append = FALSE,row.names = FALSE)

write.xlsx(output.df, file = "HIB_within_host_model_output_+No_shift+conserv_better_I+No_TMI.xlsx",
          sheetName = "1", append = FALSE,row.names = FALSE)
###########################################
#then plot models on graph
#modify names of files as needed
#un hashtag jpeg and dev.off to save as jpeg file

#jpeg(file="within_host_model_temps_colored.jpeg",width=1200,height=800)
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
#dev.off()
model_with_colors #view graph in Plots panel

##########################################################################
###########################################################################
#########################################################################
###################################################################
#Extract within host parameters for between host model
output.df.long <- data.frame() #Temp = Temp.vector, model = NA #must do this every time prior to running new models!

for (i in 1:length(fraction_data$Temp)){
  times <- seq(from=0,to=3/365,by=1/365/12) #****divide by 12 instead of 4 to add additional time intervals
  xstart<-c(H=1,I=1000,B=10000)
  
  #select fractions as parameters for each temp
  #TMI <- fraction_data[i,2]
  TMD <- fraction_data[i,3]
  
  Btoptrate<-Ratkowsky(temp=36,Tmin=6,Tmax=49,b=0.004,c=0.14)
  Btrate <- Ratkowsky(temp=fraction_data[i,1],Tmin=6,Tmax=49,b=0.004,c=0.14)
  TB <- (Btrate/Btoptrate)
  
  #parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.0001,
  #         gamma=0.1,KI=20000,alpha=0.0000085,KB=2000000,w=0.0005,z=0.001,
  #        r=6000,d=0.003,c=0.005,sigma=500000,TMI=TMI,TMD=TMD,TB=TB) #original
  
  #w/ TMI: parms <-c(psi=0.5, KH=1, beta=0.0005,p=0.0005,m=0.00001,
  #gamma=1,KI=20000,alpha=0.00085,KB=2000000,w=0.0005,z=0.001,
  #r=6000,d=0.003,c=0.05,sigma=500000,TMI=TMI,TMD=TMD,TB=TB) #parms for high I mods - conservative changes
  
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
  
  out<-mutate(out,temp=fraction_data[i,1]) #add temp to HIB vs. time output
  output.df.long<-rbind(output.df.long,out) #store in data frame
}

#next, organize data
output.df.long<-output.df.long %>%
  gather(variable,value,-temp,-time) %>% group_by(temp)
head(output.df.long)

as.data.frame(output.df.long)

##################
#save data frame
write.xlsx(output.df.long, file = "HIB_within_host_model_output_longversion.xlsx",
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
paramdf
write.xlsx(paramdf, file = "Temp_dependent_params_withinhost_Decchanges.xlsx",
           sheetName = "1", append = FALSE,row.names = FALSE)



###########################################################################
###########################################################################
###########################################################################


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
  #Kp<-params["Kp"]
  mu<-params["mu"]
  TB<-params["TB"] #bacterial growth rate fraction (from within host / Ratkowsky models)
  #f<-params["f"] #fraction of topt Median Survival Time / bacterial growth rate fraction?
  TR<-params["TR"]
  
  #model equations
  #### developing model: #####
  #dSdt <- r*(1-(S+I)/K)-(beta*S*P)+(gamma*I)-(d*S)
  #dIdt <- (beta*S*P)-(gamma*I)-((d+delta)*I)
  #dPdt <- g*((d+delta)*I)*(1-P/Kp)-(beta*P)
  
  #dSdt <- r*(1-(S+I)/K)-(beta*S*g*P*(1-P/Kp))+(gamma*I)-(d*S) #susceptible beetles
  #dIdt <- (beta*S*g*P*(1-P/Kp))-(gamma*I)-(f*(d+delta)*I) #infected beetles
  #dPdt <- (f*(d+delta)*I)-(mu*P) #compartment P = beetles that are dead and producing spores
  ##############
  
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