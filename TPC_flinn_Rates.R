#rTPC of Flinn model of Gene expression RATES (not ddCT2 values)
####################
# clear workspace
rm(list = ls())

###########################
#load packages
library(readxl)
library(ggplot2)

#remotes::install_github("padpadpadpad/rTPC")
#install.packages("nls.multstart")

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
###################################
#Run the Flinn model by treatment group

#use nls_multstart()
#this is the function:
#flinn_1991 <- function(temp, a, b, c){
# est <- 1/(1 + a + b * temp + c * temp^2)
#return(est)
#}

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

###############
#Gaussian 

fits <- Gene_data %>%
  group_by(., Rate_Type) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(Rate~gaussian_1987(temp=Temp,rmax,topt,a),
                                                data=.x,
                                                iter = 1000,
                                                start_lower=get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                start_upper=get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                lower=get_lower_lims(.x$Temp, .x$Rate, model_name = "gaussian_1987"),
                                                upper=get_upper_lims(.x$Temp, .x$Rate, model_name = "gaussian_1987"),
                                                supp_errors = "Y",
                                                convergence_count = FALSE)))
# 




# look at output object - should show a nls fit column & data tibble column by grouping
select(fits, data, fit) 

constit_topt <- calc_params(fits$fit[[1]]) %>% select(topt)
constit_topt <- as.numeric(constit_topt)
induced_topt <- calc_params(fits$fit[[2]]) %>% select(topt)
induced_topt <- as.numeric(induced_topt)


#check the first fit to see if it worked
summary(fits$fit[[1]])
summary(fits$fit[[2]])
fits$fit[[1]]
fits$fit[[2]]

glance(fits$fit[[1]])
glance(fits$fit[[2]])


## clean up:
# get summary
info <- fits %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary)
info
# get params
params <- fits %>%
  mutate(., p = map(fit, tidy)) %>%
  unnest(p)

# get confidence intervals
CI <- fits %>%
  mutate(., cis = map(fit, confint2),
         cis = map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., Rate_Type) #%>%
  #mutate(., term = c('rmax', 'topt', 'a')) #%>%
  #ungroup() %>%
  #select(., -data, -fit)

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  mutate(., p = map(fit, augment)) %>%
  unnest(p)

#check models
select(info, Rate_Type, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions - need more data for smooth curve
new_preds <- Gene_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Gene_data, Rate_Type) %>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()
#> `summarise()` ungrouping output (override with `.groups` argument)

# create new predictions
preds2 <- fits %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Rate_Type') %>%
  group_by(., Rate_Type) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# plot
jpeg(file="Flinn model of gene expression rates.jpeg",width=550,height=350)
Flinn_model_of_rates <- ggplot() +
  geom_point(aes(Temp, Rate), size = 2, Gene_data) +
  geom_line(aes(Temp, Rate, group = Rate_Type), alpha = 0.5, preds2) +
  facet_wrap(~ Rate_Type, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4')) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.9, 0.15)) +
  labs(title=expression(paste("Immune Gene Expression Rates in ",italic("T. castaneum"))),
       subtitle="Thermal Performance Curves: Flinn Model",
       y="Expression rates (per hour)",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,24,30,34), position="bottom")
Flinn_model_of_rates
dev.off()

Flinn_model_of_rates

#######################
#Calculate a, b, and c from Flinn model for diff treatment groups:
# calculate additional traits for BtU treatment group
calc_params(fits$fit[[3]]) %>%
  # round for easy viewing
  mutate_all(round, 2)
# calculate additional traits for heat killed treatment group
calc_params(fits$fit[[2]]) %>%
  # round for easy viewing
  mutate_all(round, 2)

################
#calculate ddCT2 "rate" at optimum T of 25.19
#then calculate rate at other temps, find fraction
#not sure how useful this is...

#flinn_1991 <- function(temp, a, b, c){
# est <- 1/(1 + a + b * temp + c * temp^2)
#return(est)
#}

##BtU Treatment:
#calculate ddCT2 "rate" at topt 25.19C:
Btu2519<-flinn_1991(temp=25.19,a=0.09064983,b=-0.1477239,c=0.0029320)

#24C:
Btu24<-flinn_1991(temp=24,a=0.09064983,b=-0.1477239,c=0.0029320)
#30C:
Btu30<-flinn_1991(temp=24,a=0.09064983,b=-0.1477239,c=0.0029320)
#34C:
Btu34<-flinn_1991(temp=24,a=0.09064983,b=-0.1477239,c=0.0029320)

#calculate fractions:
Btu24frac<-Btu24/Btu2519*100



#############################
#calculate ddCT2 based on model at temps to find max ddCT and %s
neededtemps <- Gene_data %>%
  do(., data.frame(Temp = c(25.19,24,30,34), stringsAsFactors = FALSE))

predsforcalc <- fits %>%
  mutate(., p = map(fit, augment, newdata = neededtemps)) %>%
  unnest(p) %>%
  group_by(., Rate_Type) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

ddCT2atTemp<-select(predsforcalc,-data,-fit)
ddCT2atTemp<-as.data.frame(ddCT2atTemp)

#install additionall packages
#install.packages("xlsx")
#install.packages("rJava")
#install.packages("Rtools")
library(xlsx)

#save as excel spreadsheet, then calculated fractions in excel
write.xlsx(ddCT2atTemp, file = "ddCT2_calcs_from_Flinnmodel.xlsx",
           sheetName = "1", append = FALSE)

######