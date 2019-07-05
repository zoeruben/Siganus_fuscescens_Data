###### LWR and GSI by Location #####
#####LWR Bais######
#### Set-Up ####
remove(list=ls())

#Set working directory
getwd()

#Load libraries
library("tidyverse")
library("ggplot2")

#Read in data
Bais <- read.csv("Bais_Population.csv", header = TRUE)

#### Calculate LWR ####

#Log transform length and weight
Bais$logL <- log(Bais$Total.Length..cm.)
Bais$logW <- log(Bais$Total.Weight..g.)

#Run linear model with log-transformed weight and length
lm_lLlW <- lm(logW ~ logL, data = Bais)
lm_lLlW
summary(lm_lLlW) #show intercept, error, and r2


#Plot of model
log_LW_plot <- ggplot(data = Bais, aes(x = logL, y = logW)) +
  geom_point() + geom_smooth(method = 'lm')
LW_plot <- ggplot(data = Bais, aes(x = Total.Length..cm., y = Total.Weight..g.)) +
  geom_point()
log_LW_plot
LW_plot

#GSI Bais

#Calculate GSI
Bais$GSI <- (Bais$Gonad.Weight..g./Bais$Total.Weight..g.) * 100

                        