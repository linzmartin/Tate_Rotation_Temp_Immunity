#modeling bacterial growth for BtU (with Sadiq's data)

##################################################
# clear workspace
rm(list = ls())
#set wd to your project folder
setwd("C:/Users/linzm/Documents/Tate_Lab_Rotation/Tate_Rotation_Temp_Immunity")
##################################################
#load packages
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

#install.packages("growthcurver")
library(growthcurver)

#install.packages("growthrates")
library(growthrates)
#############################
#("Btu_IVG_combined.xlsx",sheet = )
###########################################################

#import growth data:
BtU_growth_data <- read_xlsx("Btu_IVG_combined.xlsx",sheet = "R")
#plot the growth rate data:
ggplot(BtU_growth_data, aes(time, OD)) +
  geom_point(aes(colour = factor(temp)),size=1.5) +
  labs(title="BtU Growth Over Time",
       y="OD 600",
       x="Time (minutes)")

#plot the natural log of the OD:
ggplot(BtU_growth_data, aes(time, log(OD))) +
  geom_point(aes(colour = factor(temp)),size=1.5) +
  labs(title="BtU Growth Over Time",
       y="ln(OD 600)",
       x="Time (minutes)")

###########################################
#find the growth rate for each temperature:

#########################################
#using growthcurver package:
#https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html

###
#subset 20 deg C temp for narrow view:
BtU_20_only <- read_xlsx("BtU_IVG_20only.xlsx")

gc_fit_20 <- SummarizeGrowth(BtU_20_only$time, BtU_20_only$sample1a)
gc_fit_20

plot(gc_fit_20)
gc_fit_20$vals
gc_fit_20$vals$r

#####
#find growth rates for each temp mean on IVG combined sheet
#BtU_growth_data <- read_xlsx("Btu_IVG_combined.xlsx",sheet = "Sheet1")
BtU_growth_data <- read_xlsx("Btu_IVG_combined.xlsx",sheet = "all_replicates")

r_values <- data.frame()
for (i in 2:dim(BtU_growth_data)[2]){
  gc_x <- SummarizeGrowth(BtU_growth_data$time, BtU_growth_data[,i])
  rval <- gc_x$vals$r
  

  temporary_r_values <- data.frame("Temp"=colnames(BtU_growth_data[,i]),"r"=rval)
  r_values <- rbind(r_values,temporary_r_values)
  if(!is.null(dev.list())) dev.off()
  plotx <- plot(gc_x,main=colnames(BtU_growth_data[,i]))
}
r_values

#?plot()

#correct data frame to remove replicate number added to each time
r_values$Temp[1:9] <- rep(20,times=9)
r_values$Temp[10:18] <- rep(22,times=9)
r_values$Temp[19:27] <- rep(24,times=9)
r_values$Temp[28:36] <- rep(26,times=9)
r_values$Temp[37:45] <- rep(28,times=9)
r_values$Temp[46:54] <- rep(30,times=9)
r_values$Temp[55:63] <- rep(32,times=9)
r_values$Temp[64:72] <- rep(34,times=9)

ggplot(r_values, aes(Temp,r)) +
  #geom_point(aes(colour = factor(Rate_Type)),size=4) +
  geom_point(colour = "black",size=1.5) +
  labs(title="BtU Growth Rates",
       y="r",
       x="Temperature (ºC)")


r_means <- group_by(r_values,Temp) %>%
  summarise(
    count = n(),
    mean = mean(r, na.rm = TRUE),
    sd = sd(r, na.rm = TRUE),
    median = median(r, na.rm = TRUE),
    IQR = IQR(r, na.rm = TRUE)
  )
r_means


ggplot(r_means, aes(Temp, mean)) +
  #facet_wrap(~Gene,labeller=label_both,scales="free_y")+
  geom_point(colour = "black",size=1.5) +
  labs(title="Mean BtU Growth Rates",
       y="r (min-1)",
       x="Temperature (ºC)")
############################
#alternative method to calculating r with rolling regression:
#https://padpadpadpad.github.io/post/calculating-microbial-growth-rates-from-od-using-rolling-regression/
# load packages
library(tidyverse) #install.packages(tidyverse)
library(zoo) #install.packages(zoo)
library(broom) #install.packages(broom)
library(growthcurver) # install.packages(growthcurver)
library(nls.multstart) # install.packages(nls.multstart)
# remotes::install_github('padpadpadpad/MicrobioUoE)

# load example data
d <- growthcurver::growthdata %>%
  gather(., well, od, -time) %>%
  mutate(ln_od = log(od),
         log10_od = log10(od))

# have a look at the data
glimpse(d)


# filter for just a single well
d_a1 <- filter(d, well == 'A1')
d_a1
# define gompertz growth model
gompertz <- function(log10_nmax, log10_n0, mumax, t, lag){
  log10_n0 + (log10_nmax - log10_n0) * exp(-exp(mumax * exp(1) * (lag - t)/((log10_nmax - log10_n0) * log(10)) + 1))
}

# fit gompertz model
fit_gomp <- nls.multstart::nls_multstart(log10_od ~ gompertz(log10_nmax, log10_n0, mumax, t = time, lag),
                                         data = d_a1,
                                         start_lower = c(log10_nmax = -0.75, log10_n0 = -3, mumax = 0, lag = 0),
                                         start_upper = c(log10_nmax = 0.5, log10_n0 = -1, mumax = 10, lag = 25),
                                         lower = c(log10_nmax = -0.6, log10_n0 = -2, mumax = 0, lag = 0),
                                         iter = 500,
                                         supp_errors = 'Y')

# get predictions
gomp_preds <- augment(fit_gomp)

# plot on original scale
ggplot(d_a1, aes(time, od)) +
  geom_line(aes(time, 10^.fitted), gomp_preds, col = 'red') +
  geom_point() +
  theme_bw(base_size = 16) +
  labs(x = 'time (hours)',
       y = 'OD') +
  annotate(geom = 'text', x = 0, y = 0.37, label = paste('µ = ', round(coef(fit_gomp)[3], 2), ' hr-1', sep = ''), hjust = 0, size = (8))

##############################################
#####################
#model fitting:
r_fits <- r_values %>%
  #group_by(., Temp)%>%
  #nest()# %>%
  mutate(flinn = purrr::map(r_values, ~nls_multstart(r~flinn_1991(temp = Temp, a, b, c),
                                                 data = r_values,
                                                 iter = c(5,5,5),
                                                 start_lower = get_start_vals(r_values$Temp, r_values$r, model_name = 'flinn_1991') - 10,
                                                 start_upper = get_start_vals(r_values$Temp, r_values$r, model_name = 'flinn_1991') + 10,
                                                 lower = get_lower_lims(r_values$Temp, r_values$r, model_name = 'flinn_1991'),
                                                 upper = get_upper_lims(r_values$Temp, r_values$r, model_name = 'flinn_1991'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)))
#,
         gaussian = purrr::map(, ~nls_multstart(r~gaussian_1987(temp = Temp, rmax, topt, a),
                                                    data = .x,
                                                    iter = c(4,4,4),
                                                    start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'gaussian_1987') - 10,
                                                    start_upper = get_start_vals(.x$Temp, .x$r, model_name = 'gaussian_1987') + 10,
                                                    lower = get_lower_lims(.x$Temp, .x$r, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$Temp, .x$r, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         modifiedgaussian = purrr::map(data, ~nls_multstart(r~modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                                                            data = .x,
                                                            iter = c(4,4,4,4),
                                                            start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'modifiedgaussian_2006') - 10,
                                                            start_upper = get_start_vals(.x$Temp, .x$rRate, model_name = 'modifiedgaussian_2006') + 10,
                                                            lower = get_lower_lims(.x$Temp, .x$r, model_name = 'modifiedgaussian_2006'),
                                                            upper = get_upper_lims(.x$Temp, .x$r, model_name = 'modifiedgaussian_2006'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)),
         quadratic = purrr::map(data, ~nls_multstart(r~quadratic_2008(temp = Temp, a, b, c),
                                                     data = .x,
                                                     iter = c(4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'quadratic_2008') - 0.5,
                                                     start_upper = get_start_vals(.x$Temp, .x$r, model_name = 'quadratic_2008') + 0.5,
                                                     lower = get_lower_lims(.x$Temp, .x$r, model_name = 'quadratic_2008'),
                                                     upper = get_upper_lims(.x$Temp, .x$r, model_name = 'quadratic_2008'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         ratkowsky = purrr::map(data, ~nls_multstart(r~ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'ratkowsky_1983') - 10,
                                                     start_upper = get_start_vals(.x$Temp, .x$r, model_name = 'ratkowsky_1983') + 10,
                                                     lower = get_lower_lims(.x$Temp, .x$r, model_name = 'ratkowsky_1983'),
                                                     upper = get_upper_lims(.x$Temp, .x$r, model_name = 'ratkowsky_1983'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         rezende = purrr::map(data, ~nls_multstart(r~rezende_2019(temp = Temp, q10, a,b,c),
                                                   data = .x,
                                                   iter = c(4,4,4,4),
                                                   start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'rezende_2019') - 10,
                                                   start_upper = get_start_vals(.x$Temp, .x$r, model_name = 'rezende_2019') + 10,
                                                   lower = get_lower_lims(.x$Temp, .x$r, model_name = 'rezende_2019'),
                                                   upper = get_upper_lims(.x$Temp, .x$r, model_name = 'rezende_2019'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         spain = purrr::map(data, ~nls_multstart(r~spain_1982(temp = Temp, a,b,c,r0),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$Temp, .x$r, model_name = 'spain_1982') - 1,
                                                 start_upper = get_start_vals(.x$Temp, .x$r, model_name = 'spain_1982') + 1,
                                                 lower = get_lower_lims(.x$Temp, .x$r, model_name = 'spain_1982'),
                                                 upper = get_upper_lims(.x$Temp, .x$r, model_name = 'spain_1982'),
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

#new_preds <- BtUInjected %>%
new_preds <- Gene_data %>%
  do(., data.frame(Temp = seq(min(.$Temp), max(.$Temp), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
#max_min <- group_by(BtUInjected, Treatment) %>%
max_min <- group_by(Gene_data,Rate_Type)%>%
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
#d_labs
