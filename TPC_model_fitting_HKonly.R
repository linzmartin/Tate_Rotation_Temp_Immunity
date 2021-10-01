#Immune Gene Expression and Temperature
#Model selection process 

#fitting thermal performance curves to the data (TPCs)

#plotting multiple models to data (all on one graph)

# clear workspace
rm(list = ls())

##################################################
#load packages
library(readxl)
library(rTPC) #remotes::install_github("padpadpadpad/rTPC")
library(nls.multstart) #install.packages("nls.multstart")
library(broom)
library(tidyverse)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(ggplot2)
library(ggrepel)
library(MuMIn)
#############################
#import the data
Gene_data <- read_xlsx("Gene_exp_rates_Sadiq_081221(new).xlsx",sheet="Immune_Rate_Data")#using Sadiq's new data
Gene_data <- Gene_data[complete.cases(Gene_data),]
# looking at rate of gene expression change (delta ddct2 / time) vs. temp vs. treatment
##################################
#calculate means:
Gene_exp_means <- group_by(Gene_data,Gene,Treatment,Temp) %>%
  summarise(
    count = n(),
    mean = mean(Rate, na.rm = TRUE),
    sd = sd(Rate, na.rm = TRUE),
    median = median(Rate, na.rm = TRUE),
    IQR = IQR(Rate, na.rm = TRUE)
  )
head(Gene_exp_means) #check output

#plot the means
ggplot(Gene_exp_means, aes(Temp, mean, shape=factor(Gene))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  facet_wrap(Treatment~Gene,labeller=label_both,scales="free_y")+
  geom_point(colour = "black",size=1.5) +
  labs(title="Mean Immune Gene Expression in T. castaneum",
       y="Rate of Expression (ddCT2/time))",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))

##
#BtUInjected <- filter(Gene_data, Treatment == "btu")
HeatKilled <- filter(Gene_data, Treatment == "hk_btu")

#########################
###########################
#############################
#FITTING MULTIPLE MODELS AT ONCE:
#Run diff models on Heat Killed group
Gene_fits <- HeatKilled %>% 
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
        modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         #rezende = purrr::map(data, ~nls_multstart(Rate~rezende_2019(temp = Temp, q10, a,b,c),
          #                                  data = .x,
           #                                 iter = c(4,4,4,4),
            #                                start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') - 10,
             #                               start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') + 10,
              #                              lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
               #                             upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                #                            supp_errors = 'Y',
                 #                           convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                          lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                          upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                          supp_errors = 'Y', convergence_count = FALSE)))

glimpse(select(Gene_fits, 1:7))
#Gene_fits$gaussian[[1]]


label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}


## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack
# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

preds <- d_stack %>%
 mutate(., p = map(fit, augment)) %>%
  unnest(p)
#select(info, Treatment, logLik, AIC, BIC, deviance, df.residual)

#new_preds <- BtUInjected %>%
new_preds <- Gene_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Gene_data,Treatment,Gene)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()
#d_labs

# plot
library(ggplot2)
library(ggrepel)
ggplot() +
  geom_point(aes(Temp, Rate),size=2, Gene_data)+
  geom_line(aes(Temp,Rate,col = model_name),alpha=0.05,preds2) +
  #facet_wrap(~Treatment,labeller = labeller(.multi_line = FALSE))+
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'Immune Gene Expression across Temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type='qual', palette = 2)


################################################
#Summary of different model fits (e.g. AIC, BIC)
#install.packages("MuMIn")

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic
 
######################################################################
######################################################################
#separate HeatKilled Data by gene:
Cactus1_exp <- filter(HeatKilled, Gene == "Cactus1")
Relish2_exp <- filter(HeatKilled, Gene == "Relish2")
PGRPSC2_exp <- filter(HeatKilled, Gene == "PGRPSC2")
Defensin1dg2_exp <- filter(HeatKilled, Gene == "Def1dg2")
Hsp27_exp <- filter(HeatKilled, Gene == "hsp27")

#now, fit models for each gene within the Heat Killed treatment group:

#Defensin fitting models:
Gene_fits <- Defensin1dg2_exp %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(Rate~rezende_2019(temp = Temp, q10, a,b,c),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 supp_errors = 'Y', convergence_count = FALSE)))

#glimpse(select(Gene_fits, 1:7))
#Gene_fits$gaussian[[1]]


label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}


## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack

new_preds <- Defensin1dg2_exp %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Defensin1dg2_exp,Treatment)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot() +
  geom_point(aes(Temp, Rate), size=2,Defensin1dg2_exp)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'Defensin1dg2 Expression vs. Temp \nin Heat Killed BtU Treatment') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic
###########################
#Cactus model fitting:
Gene_fits <- Cactus1_exp %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(Rate~rezende_2019(temp = Temp, q10, a,b,c),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 supp_errors = 'Y', convergence_count = FALSE)))

## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack

new_preds <- Cactus1_exp %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Cactus1_exp,Treatment)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot() +
  geom_point(aes(Temp, Rate), size=2,Cactus1_exp)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'Cactus1 Expression vs. Temp \nin Heat Killed BtU Treatment') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic

###########################
#Relish model fitting:
Gene_fits <- Relish2_exp %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(Rate~rezende_2019(temp = Temp, q10, a,b,c),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 supp_errors = 'Y', convergence_count = FALSE)))

## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack

new_preds <- Relish2_exp %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Relish2_exp,Treatment)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

#plot:
ggplot() +
  geom_point(aes(Temp, Rate), size=2,Relish2_exp)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'Relish2 Expression vs. Temp \nin Heat Killed BtU Treatment') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic

###########################
#PGRPSC2_exp model fitting:
Gene_fits <- PGRPSC2_exp %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(Rate~rezende_2019(temp = Temp, q10, a,b,c),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'rezende_2019') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'rezende_2019'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 supp_errors = 'Y', convergence_count = FALSE)))


## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack

new_preds <- PGRPSC2_exp %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(PGRPSC2_exp,Treatment)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

#plot:
ggplot() +
  geom_point(aes(Temp, Rate), size=2,PGRPSC2_exp)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'PGRPSC2 Expression vs. Temp \nin Heat Killed Treatment') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic


##################
#HSP27 model fitting:
Gene_fits <- Hsp27_exp %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(Rate~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'briere2_1999') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'briere2_1999'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(Rate~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(Rate~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(Rate~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(Rate~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(Rate~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$Rate, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$Rate, model_name = 'spain_1982'),
                                                 supp_errors = 'Y', convergence_count = FALSE)))

## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', flinn:spain)
d_stack

new_preds <- Hsp27_exp %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Hsp27_exp,Treatment)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

#plot:
ggplot() +
  geom_point(aes(Temp, Rate), size=2,Hsp27_exp)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'HSP27 Expression vs. Temp \nin Heat Killed Treatment') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic
