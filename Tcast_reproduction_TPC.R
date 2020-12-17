#calculating fraction of optimal reproduction rate
#using Park 1948 data

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
Repro_data <- read_xlsx("Tcast_reproduction_Park1948data.xlsx") 
###################################
#plotting and testing quadratic line of best fit:
plot(Repro_data$Temp, Repro_data$Mean_eggs, main="Reproductive rate (eggs/30days) vs. Temp", 
     xlab="Temperature", ylab="Eggs/30days", pch=19)
par(new=TRUE)
lines(Repro_data$Temp,quadratic,col="red")

ggplot(Repro_data, aes(x=Temp,y=Mean_eggs))+
  geom_point() + geom_smooth(method="lm",
                             formula= y ~ poly(x, 2),se=FALSE)
######################################
#FITTING MULTIPLE MODELS AT ONCE: Gaussian and Quadratic Models 
Repro_fits <- nest(Repro_data,data=c(Temp, Mean_eggs))%>%
  mutate(gaussian = purrr::map(data, ~nls_multstart(Mean_eggs~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = Repro_data,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(Repro_data$Temp, Repro_data$Mean_eggs, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Mean_eggs~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Mean_eggs, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Mean_eggs, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Mean_eggs, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Mean_eggs, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))


label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

## clean up:
# stack models
d_stack <- select(Repro_fits, -data,-Time_days,-Source) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', gaussian:quadratic)
d_stack
# get parameters using tidy
params <- d_stack %>%
 mutate(., est = map(fit, tidy)) %>%
select(-fit) %>%
unnest(est)

# get predictions using augment
new_preds <- Repro_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- summarise(Repro_data, min_Temp = min(Temp), max_Temp = max(Temp))

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Mean_eggs = .fitted)

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()
#d_labs

##################
# plot both models on the same graph:
ggplot() +
  geom_point(aes(Temp, Mean_eggs), size=2, Repro_data)+
  #geom_line(aes(Temp,ddCT2),alpha=0.5,preds2)+
  geom_line(aes(Temp,Mean_eggs,col = model_name),preds2) +
  #geom_label_repel(aes(Temp, preds2$.fitted, label = model_name, col = model_name), 
   #                fill = 'white', nudge_y = 20, segment.size = 0.2, 
    #               segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Mean Eggs per 30 days',
       title = 'Reproduction rate across temperatures') +
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
gaussian_AIC <- d_ic$AIC[[1]]
quadratic_AIC <- d_ic$AIC[[2]]
#Gaussian model has lower AIC - preferred

##########################
#extract params for Gaussian model:
# rate = rmax.exp(-0.5.(abs(temp - topt)/a)^2) #Guassian model
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

coeff_list<-coeffs(gausfit)

Repro_rmax<-coeff_list[[1]]
Repro_topt<-coeff_list[[2]]
Repro_a<-coeff_list[[3]]
##########################################################