#Btu curves by plate
#calculating growth rates using rolling regression
#https://padpadpadpad.github.io/post/calculating-microbial-growth-rates-from-od-using-rolling-regression/
#https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html#a-simple-first-example

# load packages
library(tidyverse) #install.packages(tidyverse)
library(zoo) #install.packages(zoo)
library(broom) #install.packages(broom)
library(growthcurver) # install.packages(growthcurver)
library(nls.multstart) # install.packages(nls.multstart)

#import the data:

#background corrected:
Btu20growth <- read.delim(file = 'BtU_growth/2021_06_03_btu_IVG_20_MSR.tab', header = TRUE,
                        stringsAsFactors = FALSE)

#not background corrected but temp and blank included in columns:
Btu20growth_withblank <- read.delim(file = 'BtU_growth/2021_06_03_btu_IVG_20_MSR_withblank.tab', header = TRUE,
                          stringsAsFactors = FALSE)

Btu20growth_withblank <- subset(Btu20growth_withblank, select=-c(2)) #remove Temp column


time_list <- seq(30,1440,by=30)
Btu20growth_withblank$Time <- time_list
Btu20growth$Time <- time_list

#generate growth data
mydata <- Btu20growth %>%
  gather(., well, od, -Time) %>%
  mutate(ln_od = log(od),
         log10_od = log10(od))

glimpse(mydata)

#######################################################
#Fit modified gompertz:

# filter for just a single well
mydata_D6 <- filter(mydata, well == 'D6')
mydata_D6 <- filter(mydata_D6, ln_od != "-Inf" & ln_od != "NaN")

# define gompertz growth model
gompertz <- function(log10_nmax, log10_n0, mumax, t, lag){
  log10_n0 + (log10_nmax - log10_n0) * exp(-exp(mumax * exp(1) * (lag - t)/((log10_nmax - log10_n0) * log(10)) + 1))
}

# fit gompertz model
fit_gomp <- nls.multstart::nls_multstart(log10_od ~ gompertz(log10_nmax, log10_n0, mumax, t = Time, lag),
                                         data = mydata_D6,
                                         start_lower = c(log10_nmax = -0.75, log10_n0 = -3, mumax = 0, lag = 0),
                                         start_upper = c(log10_nmax = 0.5, log10_n0 = -1, mumax = 10, lag = 25),
                                         lower = c(log10_nmax = -0.6, log10_n0 = -2, mumax = 0, lag = 270),
                                         iter = 500,
                                         supp_errors = 'Y')
summary(fit_gomp)
#plot(mydata_D6$Time,mydata_D6$log10_od)
# get predictions
gomp_preds <- augment(fit_gomp)

coef(fit_gomp)

# plot on original scale
ggplot(mydata_D6, aes(Time, od)) +
  geom_line(aes(Time, 10^.fitted), gomp_preds, col = 'red') +
  geom_point() +
  theme_bw(base_size = 16) +
  labs(x = 'Time (minutes)',
       y = 'OD') +
  annotate(geom = 'text', x = 0, y = 0.37, label = paste('µ = ', round(coef(fit_gomp)[3], 2), ' hr-1', sep = ''), hjust = 0, size = 6)
#################################################################
# create the rolling regression function
roll_regress <- function(x){
  temp <- data.frame(x)
  mod <- lm(temp)
  temp <- data.frame(slope = coef(mod)[[2]],
                     slope_lwr = confint(mod)[2, ][[1]],
                     slope_upr = confint(mod)[2, ][[2]],
                     intercept = coef(mod)[[1]],
                     rsq = summary(mod)$r.squared, stringsAsFactors = FALSE)
  return(temp)
}

# define window - here every ~1.5 hours, measurements taken every 1/2 hour
num_points = ceiling(1.5*60/(60*0.5)) 

# run rolling regression on ln od ~ time
models <- mydata_D6 %>%
  do(cbind(model = select(., ln_od, Time) %>% 
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = select(., Time),
           ln_od = select(., ln_od))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')

# create predictions
preds <- models %>%
  filter(., !is.na(slope)) %>%
  group_by(Time) %>%
  do(data.frame(time2 = c(.$Time - 2, .$Time + 2))) %>%
  left_join(., models) %>%
  mutate(pred = (slope*time2) + intercept)


# calculate the exponential growth rate
growth_rate <- filter(models, slope == max(slope, na.rm = TRUE))

# plot rolling regression on original scale
ggplot(mydata_D6, aes(Time, ln_od)) +
  geom_point() +
  geom_line(aes(time2, pred, group = Time), col = 'red', preds, alpha = 0.5,size=2) +
  theme_bw(base_size = 16) +
  geom_segment(aes(x = Time, y = -3, xend = Time, yend = ln_od), growth_rate) +
  geom_segment(aes(x = 0, y = ln_od, xend = Time, yend = ln_od), growth_rate) +
  annotate(geom = 'text', x = 0, y = -1, 
           label = paste('µ = ', round(growth_rate$slope, 2), ' hr-1\n95%CI:(',round(growth_rate$slope_lwr, 2), '-', 
                         round(growth_rate$slope_upr, 2), ')', 
                         sep = ''), hjust = 0,size=5) +
  labs(x = 'Time (minutes)',
       y = 'OD')
#########################################################
#run rolling reg on all wells using group_by():
mydata <- filter(mydata, ln_od != "-Inf" & ln_od != "NaN")
# run rolling regression on ln od_cor ~ time
models <- mydata %>%
  group_by(well) %>%
  do(cbind(model = select(., ln_od, Time) %>% 
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = select(., Time),
           ln_od = select(., ln_od))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')

# calculate growth rate for each one
growth_rates <- models %>%
  filter(slope == max(slope, na.rm = TRUE)) %>%
  ungroup()

glimpse(growth_rates)


growth_rates$temp <- 20
sample_list <- c(1,1,1,2,2,2,3,3,3)
replicate_list <- c("a","a","a","b","b","b","c","c","c")
growth_rates$sample <- sample_list
growth_rates$replicate <- replicate_list
glimpse(growth_rates)

Bt20 <- growth_rates


#########################################
growth_by_plate_data <- SummarizeGrowthByPlate(Btu20growth_withblank,bg_correct = "blank")

growth_by_plate_data <- SummarizeGrowthByPlate(Btu20growth_withblank,bg_correct = "blank",plot_fit = TRUE,
                                 plot_file = "growthcurves_20C_plots.pdf")

head(growth_by_plate_data)
?SummarizeGrowthByPlate()
Bt20_nonrr <- growth_by_plate_data
Growth_at20 <- subset(growth_by_plate_data, select=c(sample,k,n0,r,sigma))

#############################################
