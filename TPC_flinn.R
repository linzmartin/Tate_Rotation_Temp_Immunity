#Immune Gene Expression and Temperature
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
Gene_data<-na.omit(Gene_data) #remove NAs first

get_model_names()
?gaussian_1987()

fits <- Gene_data %>%
  group_by(., Treatment) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~nls_multstart(ddCT2~flinn_1991(temp = Temp, a, b, c),
                                                 data = .x,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(.x$Temp, .x$ddCT2, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(.x$Temp, .x$ddCT2, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)))

head(fits)

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
       subtitle="Flinn model thermal performance curves",
       y="Relative expression (ddCT2)",
       x="Temperature (ºC)") +
  scale_x_continuous(breaks=c(20,24,30,34), position="bottom")
