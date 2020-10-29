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
Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Gene_Data") #Emily's data
  # lookinga at ddCT2 (fold change) vs. temp vs. treatment
  # this is condensed data from Emily's work - find in my box folder
##################################

#Plot the data to view the general trends for Treatment across temperatures
ggplot(Gene_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5)

ggplot(Gene_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) +
  geom_point(colour = "grey90",size=1.5)

# filter Gene_data to separate by Treatment and by Temperature
Gene_data<-na.omit(Gene_data) #remove NAs first

BtUInjected <- filter(Gene_data, Treatment == "BtU")
BtUInjected
NonInjected <- filter(Gene_data, Treatment == "non-injected")
HeatKilled <- filter(Gene_data, Treatment == "heat killed")
#Temp24 <- filter(Gene_data, Temp == 24)
#Temp24Btu <- filter(Gene_data, Temp == 24, Treatment == "BtU")

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
  

ggplot(Btumean, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5) +
  labs(x = 'Temperature (ºC)',
       y = 'Mean Fold Change in Immune Gene Expression',
       title = 'Immune Gene Expression Across Temperatures') +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) 

ggplot(BtUInjected, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5) +
  labs(x = 'Temperature (ºC)',
       y = 'Mean Fold Change in Immune Gene Expression',
       title = 'Immune Gene Expression Across Temperatures') +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2)) 

get_model_names()
# "briere2_1999" - doesn't work
# "gaussian_1987"  
# "sharpeschoolfull_1981" and  "sharpeschoolhigh_1981" - don't work
# "modifiedgaussian_2006" - doesn't work
# "weibull_1995"
?get_model_names()
??rTPC
#BtU only but w/ individ points
ggplot(BtUInjected, aes(Temp, ddCT2)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Fold Change in Immune Gene Expression',
       title = 'Immune Gene Expression Across Temperatures')
#mean
ggplot(Btumean, aes(Temp, ddCT2), shape=factor(Treatment)) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Fold Change in Immune Gene Expression',
       title = 'Immune Gene Expression Across Temperatures')


######
?modifiedgaussian_2006() #rate = rmax.exp(-0.5.(abs(temp - topt)/a)^b)
?gaussian_1987() #rate = rmax.exp(-0.5.(abs(temp - topt)/a)^2)
startgaus1<-get_start_vals(Btumean$Temp,Btumean$ddCT2, model_name = "modifiedgaussian_2006")
gausmod1<-nls.multstart::nls_multstart(ddCT2~modifiedgaussian_2006(temp=Temp,rmax,topt,a,b),
                                       data=Btumean,
                                       iter = c(3,3,3,3),
                                       start_lower=startgaus1 -10,
                                       start_upper=startgaus1+10,
                                       lower=get_lower_lims(Btumean$Temp,Btumean$ddCT2, model_name = "modifiedgaussian_2006"),
                                       upper=get_upper_lims(Btumean$Temp,Btumean$ddCT2, model_name = "modifiedgaussian_2006"),
                                       supp_errors = "Y",
                                       convergence_count = FALSE)

gausmod1
#mod1 doesn't work
startgaus2<-get_start_vals(Btumean$Temp,Btumean$ddCT2, model_name = "gaussian_1987")
gausmod2<-nls.multstart::nls_multstart(ddCT2~gaussian_1987(temp=Temp,rmax,topt,a),
                                       data=Btumean,
                                       iter = c(4,4,4),
                                       start_lower=startgaus2 -10,
                                       start_upper=startgaus2+10,
                                       lower=get_lower_lims(Btumean$Temp,Btumean$ddCT2, model_name = "gaussian_1987"),
                                       upper=get_upper_lims(Btumean$Temp,Btumean$ddCT2, model_name = "gaussian_1987"),
                                       supp_errors = "Y",
                                       convergence_count = FALSE)

summary(gausmod2)
#doesn't work - 46 iterations, but zero stats b/c topt is 24
preds <- data.frame(temp = seq(min(Btumean$Temp), max(Btumean$Temp), length.out = 100))
preds <- broom::augment(gausmod2, newdata = preds)

info<-glance(gausmod2)
info
preds


######################

#mod3 works! - use BtuInjected, Gaussian model
#rate = rmax.exp(-0.5.(abs(temp - topt)/a)^2)
startgaus3<-get_start_vals(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "gaussian_1987")
gausmod3<-nls.multstart::nls_multstart(ddCT2~gaussian_1987(temp=Temp,rmax,topt,a),
                                       data=BtUInjected,
                                       iter = c(4,4,4),
                                       start_lower=startgaus3 -10,
                                       start_upper=startgaus3+10,
                                       lower=get_lower_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "gaussian_1987"),
                                       upper=get_upper_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "gaussian_1987"),
                                       supp_errors = "Y",
                                       convergence_count = FALSE)

#model fit
summary(gausmod3)
info3<-glance(gausmod3)
info3
params<-tidy(gausmod3)

CI <- confint2(gausmod3) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..)

# bind params and confidence intervals
params <- bind_cols(params, CI)
select(params, -c(statistic, p.value))



#predictions
preds3 <- data.frame(temp = seq(min(BtUInjected$Temp), max(BtUInjected$Temp), length.out = 52))
preds3 <- broom::augment(gausmod3, newdata = preds3)
preds3

preds4 <- augment(gausmod3)
preds4
summary(preds3)
summary(preds4)

# plot
ggplot(preds3) +
  geom_point(aes(Temp, ddCT2), BtUInjected) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw()

ggplot() +
  geom_point(aes(Temp, ddCT2), BtUInjected) +
  geom_line(aes(Temp, .fitted), preds4)

ggplot(preds4) +
  geom_point(aes(Temp, ddCT2), BtUInjected) +
  geom_line(aes(Temp, .fitted), col = 'blue',) +
  theme_bw() +
  labs(title="Immune Gene Expression in BtU Infected T. castaneum",
       subtitle="Gaussian model for thermal performance curve",
       y="Fold change",
       x="Temperature")

#plot(preds3$.fitted~preds4$.fitted)


#########################################################
#try Weibull
#weibull_1995
# rate = ((a.(((c-1)/c)^((1-c)/c)).((((temp-topt)/b)+(((c-1)/c)^(1/c)))^(c-1)).(exp(-((((temp-topt)/b)+(((c-1)/c)^(1/c)))^c)+((c-1)/c)))))
?weibull_1995()
startW<-get_start_vals(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "weibull_1995")
weibullmod<-nls.multstart::nls_multstart(ddCT2~weibull_1995(temp=Temp,a,topt,b,c),
                                       data=BtUInjected,
                                       iter = c(4),
                                       start_lower=startW -10,
                                       start_upper=startW+10,
                                       lower=get_lower_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "weibull_1995"),
                                       upper=get_upper_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "weibull_1995"),
                                       supp_errors = "Y",
                                       convergence_count = FALSE)

##error in parameters

#model fit
summary(gausmod3)
info3<-glance(gausmod3)

#######################################
#try flinn_1991() 
#rate = 1 / (1 + a + b.temp + c.temp^2)
#a controls height of curve
#b controls slope of initial increase
#c controls position and steepness of decline of curve

?flinn_1991()

startF<-get_start_vals(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "flinn_1991")
flinnmod<-nls.multstart::nls_multstart(ddCT2~flinn_1991(temp=Temp,a,b,c),
                                         data=BtUInjected,
                                         iter = c(4),
                                         start_lower=startF -10,
                                         start_upper=startF +10,
                                         lower=get_lower_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "flinn_1991"),
                                         upper=get_upper_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "flinn_1991"),
                                         supp_errors = "Y",
                                         convergence_count = FALSE)
summary(flinnmod)
infoflinn<-glance(flinnmod)
infoflinn


params<-tidy(flinnmod)

CI <- confint2(flinnmod) %>%
  data.frame() %>%
  rename(., conf.low = X2.5.., conf.high = X97.5..)

# bind params and confidence intervals
params <- bind_cols(params, CI)
select(params, -c(statistic, p.value))

predsflinn <- augment(flinnmod)
predsflinn

ggplot(predsflinn) +
  geom_point(aes(Temp, ddCT2), BtUInjected) +
  geom_line(aes(Temp, .fitted), col = 'blue',) +
  theme_bw() +
  labs(title="Immune Gene Expression in BtU Infected T. castaneum",
       subtitle="Flinn model for thermal performance curve",
       y="Fold change",
       x="Temperature")
#slightly better AIC than Gaussian
###################################################

#boatman_2017(temp, rmax, tmin, tmax, a, b)
# rate = rmax.(sin(pi.((temp - tmin)/(tmax - tmin))^a))^b
?boatman_2017

startboat<-get_start_vals(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "boatman_2017")
boatmod<-nls.multstart::nls_multstart(ddCT2~boatman_2017(temp=Temp,rmax,tmin,tmax,a,b),
                                       data=BtUInjected,
                                       iter = c(4),
                                       start_lower=startboat -10,
                                       start_upper=startboat +10,
                                       lower=get_lower_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "boatman_2017"),
                                       upper=get_upper_lims(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "boatman_2017"),
                                       supp_errors = "Y",
                                       convergence_count = FALSE)

#doesn't work

########################################
#multiple curves at once
# https://github.com/padpadpadpad/nls.multstart

fits <- Chlorella_TRC %>%
  group_by(., flux, growth.temp, process, curve_id) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(ln.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc=-1000, E=0.1, Eh=0.5, Th=285),
                                                start_upper = c(lnc=1000, E=2, Eh=10, Th=330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))
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
       y="Fold change (ddCT2)",
       x="Temperature (ºC)")