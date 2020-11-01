#BtMid TPC

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
BtMid_data <- read_xlsx("Gene_temp_data.xlsx", sheet="BtMid_GeneData") #Emily's data
# lookinga at ddCT2 (fold change) vs. temp vs. treatment
# this is condensed data from Emily's work - find in my box folder
##################################



ggplot(BtMid_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) +
  geom_point(colour = "grey90",size=1.5) +
  labs(title="BtMid Expression")


BtMid_data<-na.omit(Gene_data) #remove NAs first


meanfoldchange<-
  aggregate(x=BtMid_data$ddCT2,
            by=list(BtMid_data$Temp,BtMid_data$Treatment),
            FUN=mean)
colnames(meanfoldchange)<- c("Temp", "Treatment", "ddCT2") #restore column names
meanfoldchange

ggplot(meanfoldchange, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) +
  geom_point(colour = "grey90",size=1.5) +
  labs(title="Mean BtMid Expression")

BtUInjected <- filter(BtMid_data, Treatment == "BtU")
BtUInjected
NonInjected <- filter(BtMid_data, Treatment == "non-injected")
HeatKilled <- filter(BtMid_data, Treatment == "heat killed")
#############################
plot(x=BtUInjected$Temp,y=BtUInjected$ddCT2)
plot(NonInjected$Temp, NonInjected$ddCT2)

linearmod<-lm(ddCT2~Temp, data=BtUInjected)
summary(linearmod)

plot(x=BtUInjected$Temp,y=BtUInjected$ddCT2)
abline(linearmod)
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