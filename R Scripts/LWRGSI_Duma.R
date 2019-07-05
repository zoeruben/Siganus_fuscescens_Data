###### LWR and GSI By Location ######
#LWR Duma

#### Set-Up ####
remove(list=ls())

#Set working directory
getwd()

#Load libraries
library("tidyverse")
library("ggplot2")

#Read in data
Duma <- read.csv("Duma_Population.csv", header = TRUE)

#### Calculate LWR ####

#Log transform length and weight
Duma$logL <- log(Duma$Total.Length..cm.)
Duma$logW <- log(Duma$Total.Weight..g.)

#Run linear model with log-transformed weight and length
lm_lLlW <- lm(logW ~ logL, data = Duma)
lm_lLlW
summary(lm_lLlW) #show intercept, error, and r2


#Plot of model
log_LW_plot <- ggplot(data = Duma, aes(x = logL, y = logW)) +
  geom_point() + geom_smooth(method = 'lm')
LW_plot <- ggplot(data = Duma, aes(x = Total.Length..cm., y = Total.Weight..g.)) +
  geom_point()
log_LW_plot
LW_plot
