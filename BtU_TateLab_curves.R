#BtU data - from Justin C.
##################
# clear workspace
rm(list = ls())
#set wd to your project folder
setwd("C:/Users/linzm/Documents/Tate_Lab_Rotation/Tate_Rotation_Temp_Immunity")
###########################
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
############################
#d <- read.table("Btu_in_vitro_growthData_30_24_Lindsay.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#bacterial logisitc growth: K / (1 + ((K - N0) / N0) * exp(-r * t))
BtU_temp_data <- read_xlsx("Btu_in_vitro_growthData_30_24_Lindsay.xlsx")
#############


ggplot(BtU_temp_data, aes(temperature, r)) +
  geom_point(aes(colour = factor(temperature)),size=4) +
  geom_smooth(mapping=aes(x=temperature,y=r)) +
  geom_point(colour = "grey90",size=1.5) +
  labs(title="BtMid Expression")

meanfoldchange<-
  aggregate(x=BtU_temp_data$r,
            by=list(BtU_temp_data$temperature),
            FUN=mean)
colnames(meanfoldchange)<- c("temperature", "r") #restore column names
meanfoldchange
################
#ratkowsky_1983()
#rate = ((a.(temp - tmin)).(1 - exp(b.(temp - tmax))))^2

fit <- nls_multstart(r~ratkowsky_1983(temp = temperature, tmin, tmax, a, b),
                     data = meanfoldchange,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(BtU_temp_data$temperature, BtU_temp_data$r, model_name = 'ratkowsky_1983') - 10,
                     start_upper = get_start_vals(BtU_temp_data$temperature, BtU_temp_data$r, model_name = 'ratkowsky_1983') + 10,
                     lower = get_lower_lims(BtU_temp_data$temperature, BtU_temp_data$r, model_name = 'ratkowsky_1983'),
                     upper = get_upper_lims(BtU_temp_data$temperature, BtU_temp_data$r, model_name = 'ratkowsky_1983'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)

fit

# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
new_data <- data.frame(temp = seq(min(meanfoldchange$temp), max(meanfoldchange$temp), 0.5))
preds <- augment(fit, newdata = new_data)

# plot data and model fit
ggplot(meanfoldchange, aes(temperature, r)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

predictions<-augment(fit)


summary(gausmod3)
info3<-glance(gausmod3)
info3
params<-tidy(fit)

CI <- confint2(fit) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..)

# bind params and confidence intervals
params <- bind_cols(params, CI)
select(params, -c(statistic, p.value))

ggplot(predictions) +
  geom_point(aes(temperature, r)) +
  geom_line(aes(temperature, .fitted), col = 'blue',) +
  theme_bw() +
  labs(title="Immune Gene Expression in BtU Infected T. castaneum",
       subtitle="Gaussian model for thermal performance curve",
       y="Fold change",
       x="Temperature")
