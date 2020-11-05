#Immune Gene Expression and Temperature
#Model selection process 
#tested models w/ BtuInjected treatment group
#plotting multiple models to data (all on one graph)

##################################################
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
#install.packages("nlstools")
library(nlstools)
#############################
#import the data
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Gene_Data_Modified") #Emily's data
# lookinga at ddCT2 (fold change) vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
##################################
#Plot the data to view the general trends for Treatment across temperatures
ggplot(Gene_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5) +
  labs(title="Immune Gene Expression in T. castaneum",
       y="Relative Expression (ddCT2)",
       x="Temperature (ºC)") +
  scale_x_discrete(limits=c(20,24,30,34), position="bottom")

ggplot(Gene_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) +
  geom_point(colour = "grey90",size=1.5)

# filter Gene_data to separate by Treatment and by Temperature
Gene_data<-na.omit(Gene_data) #remove NAs first

##
BtUInjected <- filter(Gene_data, Treatment == "BtU")
BtUInjected
NonInjected <- filter(Gene_data, Treatment == "non-injected")
HeatKilled <- filter(Gene_data, Treatment == "heat killed")
#Temp24 <- filter(Gene_data, Temp == 24)
#Temp24Btu <- filter(Gene_data, Temp == 24, Treatment == "BtU")

##########
#if you want means only (instead of individual data points) use this:
meanfoldchange<-
  aggregate(x=Gene_data$ddCT2,
            by=list(Gene_data$Temp,Gene_data$Treatment),
            FUN=mean)
colnames(meanfoldchange)<- c("Temp", "Treatment", "ddCT2") #restore column names
meanfoldchange


Btumean<-filter(meanfoldchange, Treatment=="BtU")

#map mean fold change across temperatures
ggplot(meanfoldchange, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5) +
  labs(x = 'Temperature (ºC)',
       y = 'Mean Fold Change in Immune Gene Expression',
       title = 'Immune Gene Expression Across Temperatures') +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) 
#######################
#Fit Gaussian model to Gene data
startgausall<-get_start_vals(Gene_data$Temp,Gene_data$ddCT2, model_name = "gaussian_1987")
fits <- Gene_data %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ddCT2~gaussian_1987(temp=Temp,rmax,topt,a),
                                                data=.x,
                                                iter = 1000,
                                                start_lower=startgausall -10,
                                                start_upper=startgausall+10,
                                                lower=get_lower_lims(Gene_data$Temp,Gene_data$ddCT2, model_name = "gaussian_1987"),
                                                upper=get_upper_lims(Gene_data$Temp,Gene_data$ddCT2, model_name = "gaussian_1987"),
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

# get confidence intervals
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
new_preds <- Gene_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(Gene_data, Treatment) %>%
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
  geom_point(aes(Temp, ddCT2), size = 2, Gene_data) +
  geom_line(aes(Temp, ddCT2, group = Treatment), alpha = 0.5, preds2) +
  facet_wrap(~ Treatment, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4')) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.9, 0.15)) +
  labs(title="Immune Gene Expression in T. castaneum",
       subtitle="Gaussian thermal performance curves",
       y="Relative expression (ddCT2)",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,24,30,34), position="bottom")


#########################
###########################
#############################
#FITTING MULTIPLE MODELS AT ONCE:

#Run diff models on BtUInjected group
Gene_fits <- BtUInjected %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(flinn = purrr::map(data, ~nls_multstart(ddCT2~flinn_1991(temp = Temp, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = purrr::map(data, ~nls_multstart(ddCT2~gaussian_1987(temp = Temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
        modifiedgaussian = purrr::map(data, ~nls_multstart(ddCT2~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'modifiedgaussian_2006') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'modifiedgaussian_2006') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(ddCT2~quadratic_2008(temp = Temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(ddCT2~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(ddCT2~rezende_2019(temp = Temp, q10, a,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'rezende_2019') - 10,
                                            start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'rezende_2019') + 10,
                                            lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'rezende_2019'),
                                            upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'rezende_2019'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(ddCT2~spain_1982(temp = Temp, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'spain_1982') + 1,
                                          lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'spain_1982'),
                                          upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'spain_1982'),
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
# get parameters using tidy
#params <- d_stack %>%
 # mutate(., est = map(fit, tidy)) %>%
  #select(-fit) %>%
  #unnest(est)

#preds <- d_stack %>%
 # mutate(., p = map(fit, augment)) %>%
#  unnest(p)
#select(info, Treatment, logLik, AIC, BIC, deviance, df.residual)

new_preds <- BtUInjected %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(BtUInjected, Treatment) %>%
  summarise(., min_Temp = min(Temp), max_Temp = max(Temp)) %>%
  ungroup()

#MAKE PREDS W NEW TEMPS:
preds2 <- d_stack %>%
  mutate(., p = map(fit, augment, newdata = new_preds)) %>%
  unnest(p) %>%
  merge(., max_min, by = 'Treatment') %>%
  group_by(., Treatment) %>%
  filter(., Temp > unique(min_Temp) & Temp < unique(max_Temp)) %>%
  rename(., ddCT2 = .fitted) %>%
  ungroup()

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 34) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()
#d_labs

# plot
library(ggplot2)
library(ggrepel)
ggplot() +
  geom_point(aes(Temp, ddCT2), size=2, BtUInjected)+
  #geom_line(aes(Temp,ddCT2),alpha=0.5,preds2)+
  geom_line(aes(Temp,ddCT2,col = model_name),preds2) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), 
                   fill = 'white', nudge_y = 20, segment.size = 0.2, 
                   segment.colour = "grey50",d_labs) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'left') +
  labs(x = 'Temperature (ºC)',
       y = 'ddCT2',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

################################################
#Summary of different model fits (e.g. AIC, BIC)
#install.packages("MuMIn")
library(MuMIn)  
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)
d_ic
 
