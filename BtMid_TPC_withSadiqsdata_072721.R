#BtMid TPC with Sadiq's new data:

# clear workspace
rm(list = ls())

#load libraries:

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
#import the data
BtMid_data <- read_xlsx("Gene_exp_rates_Sadiq_072721.xlsx",sheet="Btmid_Rate_Data")#using Sadiq's new data

##################################
#Plot:
ggplot(BtMid_data, aes(Temp, Rate)) +
  #geom_point(aes(colour = factor(Rate_Type)),size=4) +
  geom_point(colour = "black",size=1.5) +
  labs(title="BtMid Gene Expression in T. castaneum",
       y="Rate of Expression (ddCT2/time))",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))



meanfoldchange<-
  aggregate(x=BtMid_data$Rate,
            by=list(BtMid_data$Temp,BtMid_data$Rate_Type),
            FUN=mean)
colnames(meanfoldchange)<- c("Temp", "Rate_Type", "Rate") #restore column names
meanfoldchange

ggplot(meanfoldchange, aes(Temp, Rate, shape=factor(Rate_Type))) +
  geom_point(aes(colour = factor(Rate_Type)),size=4) +
  geom_smooth(mapping=aes(x=Temp,y=Rate)) +
  geom_point(colour = "grey90",size=1.5) +
  labs(title="Mean BtMid Expression")


#####################################

startgausall<-get_start_vals(meanfoldchange$Temp,meanfoldchange$ddCT2, model_name = "gaussian_1987")
fits <- meanfoldchange %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ddCT2~gaussian_1987(temp=Temp,rmax,topt,a),
                                                data=.x,
                                                iter = 1000,
                                                start_lower=startgausall -10,
                                                start_upper=startgausall+10,
                                                lower=get_lower_lims(meanfoldchange$Temp,meanfoldchange$ddCT2, model_name = "gaussian_1987"),
                                                upper=get_upper_lims(meanfoldchange$Temp,meanfoldchange$ddCT2, model_name = "gaussian_1987"),
                                                supp_errors = "Y",
                                                convergence_count = FALSE)))
# look at output object - should show a nls fit column & data tibble column by grouping
select(fits, data, fit) 
#check the first fit to see if it worked
summary(fits$fit[[1]])
glance(fits$fit[[1]])

## clean up:
# get summary
info <- fits %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary)

# get params
params <- fits %>%
  mutate(., p = map(fit, tidy)) %>%
  unnest(p)

# get confidence intervals - doesn't work for mean b/c don't have residuals
CI <- fits %>%
  mutate(., cis = map(fit, confint2),
         cis = map(cis, data.frame)) %>%
  unnest(cis) %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..) %>%
  group_by(., Treatment) %>%
  mutate(., term = c('rmax', 'topt', 'a')) %>%
  ungroup() %>%
  select(., -data, -fit)

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  mutate(., p = map(fit, augment)) %>%
  unnest(p)

#check models
select(info, Treatment, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions - need more data for smooth curve
new_preds <- meanfoldchange %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(meanfoldchange, Treatment) %>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()
#> `summarise()` ungrouping output (override with `.groups` argument)

# create new predictions
preds2 <- fits %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., ddCT2 = .fitted) %>%
  ungroup()


# plot
ggplot() +
  geom_point(aes(Temp, ddCT2), size = 2, meanfoldchange) +
  geom_line(aes(Temp, ddCT2, group = Treatment), alpha = 0.5, preds2) +
  facet_wrap(~ Treatment, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4')) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.9, 0.15)) +
  labs(title="BtMid Gene Expression in T. castaneum",
       subtitle="Gaussian thermal performance curves",
       y="Mean Fold change (ddCT2)",
       x="Temperature (ºC)")
#########################################################

##################
#Btmid data model fitting:
Gene_fits <- BtMid_data %>%
  group_by(., Rate_Type) %>%
  nest() %>%
  mutate(briere2 = purrr::map(data, ~nls_multstart(Rate~briere2_1999(temp = Temp, tmin, tmax, a,b),
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

label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

## clean up:
# stack models
d_stack <- select(Gene_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere2:spain)
d_stack

new_preds <- BtMid_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(BtMid_data,Rate_Type)%>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Rate_Type') %>%
  group_by(., Rate_Type) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., Rate = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(preds2, Temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
library(ggplot2)
library(ggrepel)
ggplot() +
  geom_point(aes(Temp, Rate), size=2,BtMid_data)+
  #geom_line(aes(Temp,ddCT2),alpha=0.5,preds2)+
  geom_line(aes(Temp,Rate,col = model_name),preds2) +
  geom_label_repel(aes(Temp, Rate, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'Rate of Gene Expression',
       title = 'BtMid Expression vs. Temp') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic
