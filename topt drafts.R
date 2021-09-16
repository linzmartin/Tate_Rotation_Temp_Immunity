#must run these empty frames first to initialize:
appendtimes<-data.frame() 
outputdataframetotals<-data.frame()
HIB_at_times_alltemps<-data.frame()
paramdfall<-data.frame()
btdata<-data.frame()
#next run the model:
for (i in 1:length(temp_shifts$temp)){
  fits <- Gene_data %>% #Use "shifted_gene_data" if wanting to model shifted data
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
    
    Btoptrate<-Ratkowsky(temp=(36-i),Tmin=(6-i),Tmax=(49-i),b=0.004,c=0.14)
    #For Ratkowsky model parameters, refer to model from Auger et al. 2008
    Btrate <- Ratkowsky(temp=(fraction_data[q,1]),Tmin=(6-i),Tmax=(49-i),b=0.004,c=0.14)
    TB <- (Btrate/Btoptrate)
    paramsBt <- cbind(TB,tempshift=temp_shifts[i,1],temp=(fraction_data[q,1]))
    btdata<-rbind(btdata,paramsBt)

    paramcols<-cbind(TMI,TMD,TB)
    
    paramcols<-cbind(paramcols,temp=fraction_data[q,1]) #add temp to HIB vs. time output
    paramdf<-rbind(paramdf,paramcols) #store in data frame
  }
  paramdf<-cbind(paramdf,toptshift=i,Hosttopt=25.2,BT_topt=(36-i),Toptdiff=((36-i)-(25.2)))
  paramdfall<-rbind(paramdfall,paramdf)
  
}
head(paramdfall)
view(paramdfall)
