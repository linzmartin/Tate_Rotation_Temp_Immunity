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
  geom_point(colour = "grey90",size=1.5)


BtMid_data<-na.omit(Gene_data) #remove NAs first

BtUInjected <- filter(BtMid_data, Treatment == "BtU")
BtUInjected
NonInjected <- filter(BtMid_data, Treatment == "non-injected")
HeatKilled <- filter(BtMid_data, Treatment == "heat killed")

plot(x=BtUInjected$Temp,y=BtUInjected$ddCT2)
plot(NonInjected$Temp, NonInjected$ddCT2)

linearmod<-lm(ddCT2~Temp, data=BtUInjected)
summary(linearmod)

plot(x=BtUInjected$Temp,y=BtUInjected$ddCT2)
abline(linearmod)
