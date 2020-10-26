
library(readxl)
MosquitoLifespan <- read_xlsx("Survival Data.xlsx", sheet="AeD7L1 and WT Survival") %>% # import data, change file name/location
  mutate(survival.days = Date.of.death - Emerge.Date) #add column to table for number of days survived
my_data <- read_excel("my_file.xlsx", sheet = "data")

Gene_data <- read_xlsx("Gene_temp_data.xlsx", sheet="Gene_Data")

library(ggplot2)
?plot()
?ggplot2()


ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))

ggplot(data=Gene_data) +
  geom_point(mapping=aes(x=Temp,y=ddCT2))
ggplot(data=Gene_data) +
  geom_smooth(mapping=aes(x=Temp,y=ddCT2))


ggplot(Gene_data, aes(Temp, ddCT2, shape=factor(Treatment))) +
  geom_point(aes(colour = factor(Treatment)),size=4) +
  geom_point(colour = "grey90",size=1.5)
