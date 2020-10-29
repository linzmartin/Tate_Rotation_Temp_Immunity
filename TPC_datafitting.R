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
startgaus3<-get_start_vals(BtUInjected$Temp,BtUInjected$ddCT2, model_name = "gaussian_1987")
gausmod3<-nls.multstart::nls_multstart(ddCT2~gaussian_1987(temp=Temp,rmax,topt,a),
                                       data=BtUInjected,
                                       iter = c(4,4,4),
                                       start_lower=startgaus2 -10,
                                       start_upper=startgaus2+10,
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

plot(preds3$.fitted~preds4$.fitted)
