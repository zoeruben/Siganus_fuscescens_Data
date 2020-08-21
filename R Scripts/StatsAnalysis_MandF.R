############ANALYSIS & GRAPH BUILDING OF MALE AND FEMALE SIGANUS FUSCESCENS#############
#### Set-Up ####
rm(list=ls())

#set working directory to location of this script
setwd("C:/Users/Zoe/Texas A&M University - Corpus Christi/Bird, Chris - Ruben-Bird_PIRE-REU/Individual Project/Analysis/Spreadsheets")
#Set working directory
getwd()

#Load libraries
library(tidyverse)
library(ggplot2)
library(Rmisc)
library(ellipse)
library(nlme)
library(geomorph)
library(tidyr)
library(RColorBrewer)
library(car)
library(boot)
library(lme4)
library(lmerTest)
library(dplyr)
library(nortest)
library(Hmisc)
library(reshape2)
library(emmeans) #very good vignettes

#Read in Data
All <- read.csv("Data_Sheet_ZR_Sfuscescens_MF.csv", header = TRUE, na.strings=c(""," ", "NA"))

#fix data types
All$GSI <- as.numeric(as.character(All$GSI))
All$Gonad_Weight_g <- as.numeric(as.character(All$Gonad_Weight_g))
All$Sampling_LocationSex <- paste(All$Sampling_Location, All$Sex, sep="_")

#Subset data by location
All_Duma <- All[All$Sampling_Location == "Dumaguete",]
All_Bais <- All[All$Sampling_Location == "Bais City",]
All_Kalibo <- All[All$Sampling_Location == "Kalibo",]
All_Ayungon <- All[All$Sampling_Location == "Ayungon",]
All_Sipalay <- All[All$Sampling_Location == "Sipalay",]


#Visualize Data
#Boxplot of Length v Sex and Location
ggplot(data = All, aes(y=Total_Length_cm, x=Sex, fill=(Sex) ))+
  geom_boxplot() +
  facet_grid(. ~ Sampling_Location) +
  xlab("Sex") +
  ylab("Total Length (cm)") +
  ylim(6,20) +
  labs(fill = "Sex") +
  labs(shape="Sex", color="Sampling_Location") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5))

#Check for normality in order to run ANOVAs
str(All)
ad.test(All$logTotal_Length_cm)


#two-way ANOVA on TW/TL for males v females
twoWayANOVA_TW_TL <- aov(TW_TL ~ Sampling_Location * Sex, data = All)
summary(twoWayANOVA_TW_TL)

###TOTAL LENGTH 
full_model_continuous <- lmer(Total_Length_cm ~ Sex * Num_RegisteredFishers_Div_HabArea + (1 | Sampling_Location), data = All)
full_model_continuous_lme <- lme(Total_Length_cm ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)
boxCox(lm(Total_Length_cm ~ Sex * Num_RegisteredFishers_Div_HabArea, data = All))

full_model_continuous_lme_log <- lme(log(Total_Length_cm) ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Sampling_Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'Num_RegisteredFishers_Div_HabArea')

prediction_data <- expand.grid(Sex = c('F', 'M'), Num_RegisteredFishers_Div_HabArea = seq(min(All$Num_RegisteredFishers_Div_HabArea), max(All$Num_RegisteredFishers_Div_HabArea), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#Boxplot of Total Length v Num Fishers
ggplot(data = All, aes(y=Total_Length_cm, x=Num_RegisteredFishers_Div_HabArea, 
                       group=interaction(Sampling_Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Number of Registered Fishers per km") +
  ylab("Total Length (cm)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

####TOTAL WEIGHT V NUM FISHERS
full_model_continuous <- lmer(Total_Weight_g ~ Sex * Num_RegisteredFishers_Div_HabArea + (1 | Sampling_Location), data = All)
full_model_continuous_lme <- lme(Total_Weight_g ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)
boxCox(lm(Total_Weight_g ~ Sex * Num_RegisteredFishers_Div_HabArea, data = All))

full_model_continuous_lme_log <- lme(log(Total_Weight_g) ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Sampling_Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'Num_RegisteredFishers_Div_HabArea')



prediction_data <- expand.grid(Sex = c('F', 'M'), Num_RegisteredFishers_Div_HabArea = seq(min(All$Num_RegisteredFishers_Div_HabArea), max(All$Num_RegisteredFishers_Div_HabArea), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#Graph of Total Weight v Num Fishers
ggplot(data = All, aes(y=Total_Weight_g, x=Num_RegisteredFishers_Div_HabArea, 
                       group=interaction(Sampling_Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Number of Registered Fishers per km") +
  ylab("Total Weight (g)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###TOTAL WEIGHT/TOTAL LENGTH
full_model_continuous <- lmer(TW_TL ~ Sex * Num_RegisteredFishers_Div_HabArea + (1 | Sampling_Location), data = All)
full_model_continuous_lme <- lme(TW_TL ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)
boxCox(lm(TW_TL ~ Sex * Num_RegisteredFishers_Div_HabArea, data = All))

full_model_continuous_lme_log <- lme(log(TW_TL) ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Sampling_Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'Num_RegisteredFishers_Div_HabArea')



prediction_data <- expand.grid(Sex = c('F', 'M'), Num_RegisteredFishers_Div_HabArea = seq(min(All$Num_RegisteredFishers_Div_HabArea), max(All$Num_RegisteredFishers_Div_HabArea), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#Total Weight/Total Length v Num Fishers
ggplot(data = All, aes(y=TW_TL, x=Num_RegisteredFishers_Div_HabArea, 
                       group=interaction(Sampling_Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Number of Registered Fishers per km") +
  ylab("Total Weight/Total Length (g/cm)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))


###B-Value
full_model_continuous <- lmer(B_value ~ Sex * Num_RegisteredFishers_Div_HabArea + (1 | Sampling_Location), data = All)
full_model_continuous_lme <- lme(B_value ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)
boxCox(lm(Total_Length_cm ~ Sex * Num_RegisteredFishers_Div_HabArea, data = All))

full_model_continuous_lme_log <- lme(log(B_value) ~ Sex * Num_RegisteredFishers_Div_HabArea, random = ~1|Sampling_Location, data = All)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Sampling_Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'Num_RegisteredFishers_Div_HabArea', "Sampling_Location"), All)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'Num_RegisteredFishers_Div_HabArea')


prediction_data <- expand.grid(Sex = c('F', 'M'), Num_RegisteredFishers_Div_HabArea = seq(min(All$Num_RegisteredFishers_Div_HabArea), max(All$Num_RegisteredFishers_Div_HabArea), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#B-value v Num Fishers
ggplot(data = All, aes(y=B_value, x=Num_RegisteredFishers_Div_HabArea, 
                       group=interaction(Sampling_Location,Sex), colour = Sex)) +
  geom_boxplot() +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Number of Registered Fishers per km^2") +
  ylab("B-value") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5))

#line plot with multiple groups
ggplot(data=All) +
  aes(x=Num_RegisteredFishers_Div_HabArea, y=Total_Length_cm, group_by(Sex), fill=Sex)+
  geom_boxplot()

#two-way ANOVA on TL for males v females and by fishing pressure category
twoWayANOVA_TL <- aov(Total_Length_cm ~ Num_RegisteredFishers_Div_HabArea * Sex, data = All)
summary(twoWayANOVA_TL)         

1.34-1.16 #high TW/TL between sexes
1.49-1.19 #low TW/TL between sexes


#LWR plot for all locations, male and female
ggplot(data = All, aes(y=Total_Weight_g, x=Total_Length_cm, color=Sex)) +
  facet_grid(Num_RegisteredFishers_Div_HabArea ~ Sex) +
  geom_point() +
  stat_smooth(method = "lm", aes(color=Sex)) +
  xlab("Total Weight (g)") +
  ylab("Total Length (cm)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))


####B-Value v Number of Registered Fishers (GOOD GRAPH)
ggplot(data=All, aes(y=B_value, x=Num_RegisteredFishers, color = Sex, shape=Sex)) +
  geom_point(shape=1, size=15, stroke=5) +
  stat_smooth(method = "lm", aes(fill=Sex, colour=Sex)) +
  xlab("Number of Registered Fishers") +
  ylab("Fish Growth (B-value)") +
  labs(color="Sex") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

#Get averages of total length, total weight, and total weight div total length
#install.packages('doBy')
library(doBy)
library(psych)           
describeBy(All$B_value, group=All$Num_RegisteredFishers)
summaryBy(data=All, B_value + Total_Length_cm + Total_Weight_g + TW_TL ~ Sex * Sampling_Location)

Average_Values <- as.data.frame(summaryBy(data=All, B_value + Total_Length_cm + Total_Weight_g + TW_TL + GSI ~ Sex * Sampling_Location * Num_RegisteredFishers))           

library(tidyverse)

Average_Values %>% spread(Sex , B_value.mean)

unstack(Average_Values, Sex ~ B_value.mean + Total_Length_cm.mean + Total_Weight_g.mean + TW_TL.mean + GSI.mean + Sampling_Location)

unstack(Average_Values, Sex ~ B_value.mean)


###Average Total Length v Num Fishers
ggplot(data=Average_Values, aes(y=Total_Length_cm.mean, x=Num_RegisteredFishers, color = Sex, shape=Sex)) +
  geom_point(shape=1, size=15, stroke=5, aes(fill=Sex)) +
  stat_smooth(method = "lm", aes(color=Sex, fill=Sex)) +
  xlab("Number of Registered Fishers") +
  ylab("Average Total Length (cm)") +
  labs(color="Sex") +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###Average Total Weight v Num Fishers
ggplot(data=Average_Values, aes(y=Total_Weight_g.mean, x=Num_RegisteredFishers, color = Sex, shape=Sex)) +
  geom_point(shape=1, size=15, stroke=5, aes(fill=Sex)) +
  stat_smooth(method = "lm", aes(color=Sex, fill=Sex)) +
  xlab("Number of Registered Fishers") +
  ylab("Average Total Weight (g)") +
  labs(color="Sex") +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###Average TW/TL v Num Fishers
ggplot(data=Average_Values, aes(y=TW_TL.mean, x=Num_RegisteredFishers, color = Sex, shape=Sex)) +
  geom_point(shape=1, size=15, stroke=5, aes(fill=Sex)) +
  stat_smooth(method = "lm", aes(color=Sex, fill=Sex)) +
  xlab("Number of Registered Fishers") +
  ylab("Average Total Weight/Total Length (g/cm)") +
  labs(color="Sex") +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

#Visualize TW/TL by Sex and Number of Registered Fishers
ggplot(data = All, aes(y=TW_TL, x=Sex, fill=(Sex) ))+
  geom_boxplot() +
  facet_grid(. ~ Num_RegisteredFishers) +
  xlab("Sex") +
  ylab("Total Weight/Total Length (g/cm)") +
  ylim(0,12) +
  labs(fill = "Sex") +
  labs(shape="Sex", color="Num_RegisteredFishers") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))


DFM <- Average_Values[Average_Values$Sex == "M", 2:8]
DFF <- Average_Values[Average_Values$Sex == "F", 2:8]
DF <- cbind(DFF, DFM)
FdivM<-cbind(Average_Values[Average_Values$Sex == "M", 2:3],DFF[3:7]/DFM[3:7])
colnames(FdivM) <- c("Sampling Location", "Num_RegisteredFishers", "B_value_FdivM", "Total_Length_cm_FdivM", "Total_Weight_g_FdivM", "TW_TL_FdivM", "GSI_FdivM")
FminM<-cbind(Average_Values[Average_Values$Sex == "M", 2:3],DFF[3:7]-DFM[3:7])
colnames(FminM) <- c("Sampling Location", "Num_RegisteredFishers", "B_value_FminM", "Total_Length_cm_FminM", "Total_Weight_g_FminM", "TW_TL_FminM", "GSI_FminM")
FminM <- FminM[,3:7]
FvM_data <- cbind(FdivM, FminM)

FvM_data$Num_RegisteredFishers <- as.numeric(as.character(FvM_data$Num_RegisteredFishers))

#FUNCTION FOR PLOTTING FEMALE V MALE MORPHOMETRICS
Plot_FvM_data <- function(DF, Y, X) {
  ggplot(data=DF, aes(y=Y, x=X)) +
    geom_point(shape=1 , size=10 , stroke=5) +
    xlab("Number of Registered Fishers") +
    ylab("Female v Male") +
    theme(legend.position="bottom") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.text = element_text(size=18)) +
    theme(text = element_text(size=25))
}


#function for plotting female v male morphometrics incl string
Plot_FvM_data2 <- function(DF=FvM_data, Y, X, Ylab="") {
  ggplot(data=DF, aes_string(y=Y, x=X)) +
    geom_point(shape=1 , size=10 , stroke=5) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x, 2, raw=TRUE), color="red") +
    xlab("Number of Registered Fishers") +
    ylab(Ylab) +
    theme(legend.position="bottom") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.text = element_text(size=18)) +
    theme(text = element_text(size=25))
}
Plot_FvM_data2(Y='1/(Total_Weight_g_FdivM)', X='Num_RegisteredFishers', Ylab="Ratio Female/Male Total Weight")

Plot_FvM_data2(Y='log(Total_Weight_g_FminM)', X='Num_RegisteredFishers', Ylab="Absolute Difference in Female/Male Total Weight")

Plot_FvM_data2(DF=DFF, Y="GSI.mean", X="(Num_RegisteredFishers)") +
  ylab("Mean GSI")

###Ratio of B-value v Num Fishers
ggplot(data=FdivM, aes(y=B_value_FdivM, x=Num_RegisteredFishers)) +
  geom_point(shape=1 , size=10 , stroke=5) +
  xlab("Number of Registered Fishers") +
  ylab("Ratio of Female/Male B-value") +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###Diff of B-value v Num Fishers
ggplot(data=FvM_data, aes(y=B_value_FminM, x=Num_RegisteredFishers)) +
  geom_point(shape=1 , size=10 , stroke=5) +
  xlab("Number of Registered Fishers") +
  ylab("Difference b/w Female/Male B-value") +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))
