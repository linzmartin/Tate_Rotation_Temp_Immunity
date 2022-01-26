#Immune Gene Expression and Temperature
#Model selection process - All genes
###################################################
#fitting thermal performance curves to the data (TPCs)

#plotting multiple models to data (all on one graph)

# clear workspace
rm(list = ls())
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
#Gene_data <- read_xlsx("Gene_exp_rates_Sadiq_081221(new).xlsx",sheet="Immune_Rate_Data")#using Sadiq's new data

Gene_data <- read_xlsx("Gene_exp_rates_Sadiq_011022.xlsx",sheet="Sheet1")#using Sadiq's new data



Gene_data <- read_xlsx("Gene_exp_rates_Sadiq_011722.xlsx",sheet="Sheet1")#using Sadiq's new data

Gene_data <- read_xlsx("Gene_exp_rates_Sadiq_011722.xlsx",sheet="Sheet2")#using Sadiq's new data, rates averaged by sex


# looking at rate of gene expression change (delta ddct2 / time) vs. temp vs. treatment
##################################
#Plot the data to view the general trends for Treatment across temperatures
ggplot(Gene_data, aes(Temp, Rate, shape=factor(Gene))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "black",size=1.5) +
  labs(title="Immune Gene Expression in T. castaneum",
       y="Rate of Expression (ddCT2/time))",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))+
  facet_grid(Gene~Treatment,labeller=label_both,scales="free")


#################
#plot individual replicates:
ggplot(Gene_data, aes(Temp, Rate, shape=factor(Gene))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  facet_wrap(~Gene,labeller=label_both,scales="free")+
  geom_point(colour = "black",size=1.5) +
  labs(title="Immune Gene Expression in T. castaneum",
       y="Rate of Expression (ddCT2/time))",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))

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

ggplot(Gene_exp_means, aes(Temp, mean, shape=factor(Gene))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  facet_wrap(~Gene,labeller=label_both,scales="free_y")+
  geom_point(colour = "black",size=1.5) +
  labs(title="Mean Immune Gene Expression in T. castaneum",
       y="Rate of Expression (ddCT2/time))",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))

#separate data by treatment group:
BtUInjected <- filter(Gene_data, Treatment == "btu")
HeatKilled <- filter(Gene_data, Treatment == "hk_btu")


##########
#calculate means again (instead of individual data points), use this for following code:
meanRate<-
  aggregate(x=Gene_data$Rate,
            by=list(Gene_data$Temp,Gene_data$Treatment),
            FUN=mean)
colnames(meanRate)<- c("Temp", "Treatment", "Mean_Rate") #restore column names
meanRate

#map mean fold change across temperatures for both hk and btu
ggplot(meanRate, aes(Temp, Mean_Rate,shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5) +
  labs(x = 'Temperature (ºC)',
       y = 'Mean Change of Immune Gene Expression',
       title = 'Immune Gene Expression Rates Across Temperatures') +
  geom_smooth(mapping=aes(x=Temp,y=Mean_Rate))  +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34))

#separate live and heat-killed btu treatments:
Btumean<-filter(meanRate, Treatment=="btu")
HKmean <- filter(meanRate, Treatment == "hk_btu")
#######################
#Fit Gaussian model to Gene data
startgausall<-get_start_vals(Gene_data$Temp,Gene_data$Rate, model_name = "gaussian_1987")
fits <- Gene_data %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(Rate~gaussian_1987(temp=Temp,rmax,topt,a),
                                                data=.x,
                                                iter = 1000,
                                                start_lower=startgausall -10,
                                                start_upper=startgausall+10,
                                                lower=get_lower_lims(Gene_data$Temp,Gene_data$Rate, model_name = "gaussian_1987"),
                                                upper=get_upper_lims(Gene_data$Temp,Gene_data$Rate, model_name = "gaussian_1987"),
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
  group_by(., Rate_Type) %>%
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
  rename(., Rate = .fitted) %>%
  ungroup()


# plot
ggplot() +
  geom_point(aes(Temp, Rate), size = 2, Gene_data) +
  geom_line(aes(Temp, Rate, group = Treatment), alpha = 0.5, preds2) +
  facet_wrap(~ Treatment, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4')) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.9, 0.15)) +
  labs(title="Immune Gene Expression in T. castaneum",
       subtitle="Gaussian thermal performance curves",
       y="Gene expression rate (ddCT2/time)",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,22,24,26,28,30,32,34), position="bottom")


