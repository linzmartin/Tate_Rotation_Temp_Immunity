#Spore germination - Thermal performance

#Data set: Knaysi_1964_Fig1_Graph_grabber_data.csv

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
library(ggrepel)
library(MuMIn)  
library(broom)
#############################
#import the data and remove NAs
Germination_data <- read.csv("Knaysi_1964_Fig1_Graph_grabber_data.csv", stringsAsFactors = FALSE) 
###################################
#plotting and testing quadratic line of best fit:
plot(Germination_data$Temperature, Germination_data$Percent_Dark_Spores,main="Spore Germination vs. Temp",
     xlab="Temperature",ylab="Percent Dark Spores",
     sub="Data from Knaysi (1964) Figure 1")
###################################
#FITTING MULTIPLE MODELS AT ONCE: Gaussian and Quadratic Models #using package rTPC
Germ_fits <- nest(Germination_data,data=c(Temperature, Percent_Dark_Spores))%>%
  mutate(gaussian = purrr::map(data, ~nls_multstart(Percent_Dark_Spores~gaussian_1987(temp = Temperature, rmax, topt, a),
                                                    data = Germination_data,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(Germination_data$Temperature,Germination_data$Percent_Dark_Spores, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(Germination_data$Temperature, Germination_data$Percent_Dark_Spores, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Percent_Dark_Spores~quadratic_2008(temp = Temperature, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         flinn = purrr::map(data, ~nls_multstart(Percent_Dark_Spores~flinn_1991(temp = Temperature, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         weibull = purrr::map(data, ~nls_multstart(Percent_Dark_Spores~weibull_1995(temp = Temperature, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Percent_Dark_Spores~ratkowsky_1983(temp = Temperature, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$Temperature, .x$Percent_Dark_Spores, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)))


label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

head(Germ_fits)

## clean up:
# stack models
d_stack <- select(Germ_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', gaussian:ratkowsky)
d_stack
# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
new_preds <- Germination_data %>%
  do(., data.frame(Temperature = seq(min(.$Temperature), max(.$Temperature), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- summarise(Germination_data, min_Temp = min(Temperature), max_Temp = max(Temperature))

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min) %>%
  filter(., Temperature > unique(min_Temp) & Temperature < unique(max_Temp)) %>%
  rename(., Percent_Dark_Spores = .fitted)

##################
# plot both models on the same graph:
ggplot() +
  geom_point(aes(Temperature, Percent_Dark_Spores), size=2, Germination_data)+
  #geom_line(aes(Temp,ddCT2),alpha=0.5,preds2)+
  geom_line(aes(Temperature,Percent_Dark_Spores,col = model_name),preds2) +
  #geom_label_repel(aes(Temp, preds2$.fitted, label = model_name, col = model_name), 
  #                fill = 'white', nudge_y = 20, segment.size = 0.2, 
  #               segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Percent_Dark_Spores',
       title = 'Spore germination rate across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

################################################
#Summary of different model fits (e.g. AIC, BIC) - compare AIC values
d_ic <- d_stack %>%
  mutate(d_stack, info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic

#Ratkowsky model has lower AIC & BIC - preferred for spore germination TPC

##########################
#extract params for Ratkowsky model:
# rate = ((a.(temp - tmin)).(1 - exp(b.(temp - tmax))))^2
Ratkowsky <- function(temp, Tmin, Tmax, a,b){
  rate <- (b*(temp-Tmin)*(1-exp(c*(temp-Tmax))))^2
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
#############################################