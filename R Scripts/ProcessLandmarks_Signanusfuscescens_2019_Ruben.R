rm(list=ls())

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

#CEB forked geomorph to modify its functions.  Those mods are sourced here
#source("C:/Users/cbird/Documents/GCL/scripts/geomorph/R/plotTangentSpace.r")
# lapply(list.files(path="C:/Users/cbird/Documents/GCL/scripts/geomorph/R", pattern = "[.][Rr]$", recursive = FALSE, full.names=TRUE), source)

#pdf("regression_shape_fishingpressure_stringent_PLOTS.pdf")

#setwd("C:/Users/cbird/OneDrive - Texas A&M University - Corpus Christi/SharePoint/Bachner, Micah - PIRE/Landmarks_&_Data/!Final_MorphJ_With_Chris_Nov2018")
setwd("C:/Users/Zoe/Texas A&M University - Corpus Christi/Bird, Chris - Ruben-Bird_PIRE-REU/Individual Project/Analysis/geomorph_analysis")

############################LANDMARKS DATA HACK: TO REMOVE NA's AND TAIL#############################
landmarks_notfixed1<-readland.tps("Siganusfuscescens_landmarks_2019_Ruben_Play.TPS", specID="ID")
class(landmarks_notfixed1)
dim(landmarks_notfixed1)
#this is to remove all specimens with "NA" for Sex AND Gonadal Stage
#omit <- c(65:98, 101, 109:110, 113, 120, 124, 128, 130, 159, 163, 193:203, 264, 266, 270, 272:273,
          275, 277, 280, 283, 286, 288:289, 209, 216, 223, 228, 233:263)
#this is to remove all specimens with "NA" for Sex ONLY
omit <- c(65:96, 159, 163, 193, 196, 216, 223, 228, 237:239, 245:246, 249:250, 254, 257, 259, 263)
landmarks_notfixed2 <- landmarks_notfixed1[,,-omit]
landmarks <- landmarks_notfixed2[-(25:27),,]
dim(landmarks)


#read in data
#landmarks<-readland.tps("TotalComplete_ForAnalysis.TPS", specID="ID")
#edited landmarks in a few fish using tpsdig2 and saved to new tps file
metadata_notfixed<-read.csv("metadata_ruben.csv")
names(metadata_notfixed)
#fix metadata to exclude samples that have "NA" for Sex and Gonadal_Stage
metadata <- metadata_notfixed[complete.cases(metadata_notfixed$Sex),]
#metadata <- metadata_notfixed[complete.cases(metadata_notfixed$Gonadal_Stage),]
names(metadata)
length(metadata)
#################################################################################################
#let's make aggregate index of number of fishers relative to Siganus habitat
data_fishingpressure <- metadata[,c(24,28)]
data_fishingpressure$distanceFromPopHub<-metadata[,17]+1


INV_data_fishingpressure<-1/(data_fishingpressure)

#pairs(data_fishingpressure)
invfishingpressure_PCA<-princomp(INV_data_fishingpressure, cor = TRUE, scores=T)
summary(invfishingpressure_PCA)
plot(invfishingpressure_PCA)
plot(invfishingpressure_PCA$scores[,'Comp.1'] ~ invfishingpressure_PCA$scores[,'Comp.2'])

fishlist <- c(min(which(metadata$Location_Code == "BAIS")),
              min(which(metadata$Location_Code == "SIP")),
              min(which(metadata$Location_Code == "AYU")),
              min(which(metadata$Location_Code == "DUMA")),
              min(which(metadata$Location_Code == "KAL")))
for(i in fishlist){
  print(metadata$Location_Code[i])
  print(invfishingpressure_PCA$scores[i,'Comp.1'])
  print(invfishingpressure_PCA$scores[i,'Comp.2'])
}

metadata$invfishingpressure_PC1<-princomp(INV_data_fishingpressure, cor = TRUE, scores=T)$scores[,'Comp.1']
metadata$invfishingpressure_PC2<-princomp(INV_data_fishingpressure, cor = TRUE, scores=T)$scores[,'Comp.2']

invfishingpressure_PC1<-princomp(INV_data_fishingpressure, cor = TRUE, scores=T)$scores[,'Comp.1']
invfishingpressure_PC2<-princomp(INV_data_fishingpressure, cor = TRUE, scores=T)$scores[,'Comp.2']


#################################################################################################

#################################################################################################
#visualize raw data, make wire frame, do procrustes 
mean.shape<-mshape(landmarks)
#wireframe<-define.links(mean.shape)
#write.csv(wireframe, file= "wireframe.csv", row.names = F)
links<-as.matrix(read.csv("wireframe.csv"))
plotAllSpecimens(landmarks, links=links, label=F, plot.param = list(labels=metadata$MorphoJ_ID, txt.cex=2.5,txt.col="red",pt.cex=2, pt.bg="white", mean.cex=3, link.col="black", link.lwd=5))
procrustes<-gpagen(landmarks)
mean.shape.proc<-mshape(procrustes$coords)
plotAllSpecimens(procrustes$coords,links=links, label=F, plot.param = list(labels=metadata$MorphoJ_ID,txt.cex=2.5,txt.col="red",pt.cex=2, pt.bg="white", mean.cex=3, link.col="black", link.lwd=5))

#identify specimen closest to mean shape
findMeanSpec(procrustes$coords)

#conduct PCA on Procrustes coordinates, note: this creates a pca graph, but I prefer the ggplot version below
par(cex=1, cex.axis=2, cex.lab=2, font.lab=2, mar=c(5,5,4,2)+0.1)
pca_procrustes<-plotTangentSpace(procrustes$coords,warpgrids=F, groups=metadata$Location, legend = TRUE, pt.size=3, method="points")
summary(pca_procrustes)
#pca_procrustes<-prcomp(two.d.array(procrustes$coords))

#visualize fish shapes at PC extremes
pdf(file="PC_wireframes.pdf")
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC1max, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC1min, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC2max, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC2min, links=links, method ="points", mag=1)
# plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC3max, links=links, method ="points", mag=1)
# plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC3min, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC3max, links=links, method ="points", mag=1, gridPars=gridPar(pt.bg = "white", pt.size = 1))
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC3min, links=links, method ="points", mag=1, gridPars=gridPar(pt.bg = "white", pt.size = 1))
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC4max, links=links, method ="points", mag=1, gridPars=gridPar(pt.bg = "white", pt.size = 1))
plotRefToTarget(mean.shape.proc, pca_procrustes$pc.shapes$PC4min, links=links, method ="points", mag=1, gridPars=gridPar(pt.bg = "white", pt.size = 1))
dev.off()

#ggplot version of pca
pc.scores<-as.data.frame(pca_procrustes$pc.scores)
#ggplot version of pca
pdf(file="Rimages/PCA_sex.pdf")
ggplot(pc.scores, aes(PC1,PC2, color=metadata$Sex))+
  geom_point(size=5)+
  stat_ellipse(type="norm", level=0.5, size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  guides(color=guide_legend(title="Sex"))+
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept=0))+
  scale_color_brewer(palette="Set1")

ggplot(pc.scores, aes(PC3,PC4, color=metadata$Sex))+
  geom_point(size=5)+
  stat_ellipse(type="norm", level=0.5, size=3)+
  theme_bw()+
  theme(axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  guides(color=guide_legend(title="Sex"))+
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept=0))+
  scale_color_brewer(palette="Set1")
dev.off()

#################################################################################################

#################################################################################################
##look for outlier fish
#note, this makes the landmark id Var1, the X,Y labels Var2, the sites Var3; etc
#also note that the landmark id and xy labels come out as letters from A-Z
# coords_stacked<-as.data.frame.table(procrustes$coords, responseName = "coordinate", stringsAsFactors = F)
# coords_stacked$landmarkID<-match(coords_stacked$Var1, LETTERS)
# coords_stacked$Var1<-NULL
# coords_stacked$Var2[coords_stacked$Var2 == "X"]<- "x"
# coords_stacked$Var2[coords_stacked$Var2 == "Y"]<- "y"
# coords_stacked<-rename(coords_stacked, c("Var3" = "FishID"))
# 
# coords_unstacked<-spread(coords_stacked, Var2, coordinate)
# ggplot(coords_unstacked, aes(x=x, y=y, label=FishID))+
#   geom_point()

#plot the spread of points for each landmark and label by fishID
# for (i in seq(1:24)) {
#   print(i)
#   sub_coords<-subset(coords_unstacked,landmarkID == i)
#   print(ggplot(sub_coords, aes(x=x, y=y, label=FishID))+
#     geom_point() +
#     geom_text(aes(label=FishID),hjust=0, vjust=0) + ggtitle(paste("landmark ", i)))
# }

outliers_procrustes<-plotOutliers(procrustes$coords, groups = metadata$Location, inspect.outliers = TRUE)

#################################ALL OF THIS NEEDS ADJUSTED FOR MY DATA###################################
#visualize outliers
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"dDU091 "], links=links, method ="TPS")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"dDU091 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"dDU091 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY069 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY069 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY066 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY066 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY068 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"aAY068 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM134 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM134 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM121 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM121 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM117 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM117 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM128 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM128 "], links=links, method ="points") 
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM115 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM115 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM124 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"cAM124 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"bBA151 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"bBA151 "], links=links, method ="points")
# 
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"bBA147 "], method ="vector")
# plotRefToTarget(mean.shape.proc, procrustes$coords[,,"bBA147 "], links=links, method ="points")

#plotAllSpecimens(procrustes$coords,links=links, plot.param = list(labels=metadata$MorphoJ_ID))
#################################################################################################

#################################################################################################
#quantify morphological integration (see manual and Bookstein 2015)
globalIntegration(procrustes$coords)
#################################################################################################

#################################################################################################
#test for allometry
#library(geomorph)
allometry.geomorph<-geomorph.data.frame(procrustes,Location=metadata$Location, Sex=metadata$Sex,
                                        Mature=metadata$Mature, Gonad_Stage=metadata$Gonadal_Stage,
                                        Gonadal_Stage_Index=metadata$Gonadal_Stage_Index, 
                                        invFishingPressure_PC1=metadata$invfishingpressure_PC1)
#allometry.geomorph<-geomorph.data.frame(procrustes,Location=metadata$Location, Sex=metadata$Sex)
length(metadata$Location)
length(metadata$Sex)

#allometry.geomorph$location<-metadata$Location
#allometry.geomorph$invCsize <- 1/allometry.geomorph$Csize
allometry.geomorph$lnCsize <- log(allometry.geomorph$Csize)

#group<-allometry.geomorph$Location
group<-allometry.geomorph$Sex
#group<-allometry.geomorph$Mature
#group<-allometry.geomorph$Gonad_Stage

#LNallometry<-procD.allometry(coords~lnCsize, ~group, logsz=F, iter=9999, alpha=0.01, data=allometry.geomorph)
#plot(LNallometry)
#plot(LNallometry, method="CAC", warpgrids=F)
# plot(LNallometry, method="RegScore", warpgrids=F)
# plot(LNallometry, method="PredLine", warpgrids=F)
# LNallometry
# 
# #prepare data for ggplot
# LNallometry.df<-cbind(as.data.frame(LNallometry$gps),as.data.frame(LNallometry$CAC),as.data.frame(allometry.geomorph$lnCsize))
# LNallometry.df<-rename(LNallometry.df, c("V1" = "CAC"))
# LNallometry.df<-rename(LNallometry.df, c("allometry.geomorph$lnCsize" = "lnCsize"))
# LNallometry.df<-rename(LNallometry.df, c("LNallometry$gps" = "Group"))
# LNallometry.df<-cbind(LNallometry.df, as.data.frame(LNallometry$Reg.proj))
# LNallometry.df<-rename(LNallometry.df, c("V1" = "Reg.proj"))

# pdf(file="Rimages/allometry_figs.pdf")
# ggplot(LNallometry.df, aes(x=lnCsize, y=CAC, color=Group)) + 
#   geom_smooth(method=lm, size=2, se=TRUE, level=0.95) +
#   geom_smooth(method=loess, color="black",size=2, se=TRUE, level=0.95) +
#   geom_point(size=1, shape = 21, stroke = 1.5) +
#   theme_bw()+
#   theme(axis.text=element_text(size=18), 
#         axis.title=element_text(size=18,face="bold"),
#         legend.text=element_text(size=18),
#         legend.title=element_text(size=18,face="bold"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+
#   #panel.border = element_blank())+
#   guides(color=guide_legend(title="Sex"))+
#   scale_color_brewer(palette="PuBuGn")
# 
# ggplot(LNallometry.df, aes(x=lnCsize, y=Reg.proj, color=Group)) + 
#   geom_smooth(method=lm, size=2, se=TRUE, level=0.95) +
#   geom_smooth(method=loess, color="black",size=2, se=TRUE, level=0.95) +
#   geom_point(size=1, shape = 1, stroke = 1.5) +
#   theme_bw()+
#   theme(axis.text=element_text(size=18), 
#         axis.title=element_text(size=18,face="bold"),
#         legend.text=element_text(size=18),
#         legend.title=element_text(size=18,face="bold"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+
#   #panel.border = element_blank())+
#   guides(color=guide_legend(title="Sex"))+
#   scale_color_brewer(palette="Set1")
# dev.off()

#obtain size-adjusted resids
LNallometry.lm <- procD.lm(coords~lnCsize*Sex, data=allometry.geomorph, iter = 9999, RRPP=TRUE)
summary(LNallometry.lm) # same ANOVA Table

########################CHECK GONADAL STAGE################################################
#group<-allometry.geomorph$Gonadal_Stage_Index
#LNallometry.lm <- procD.lm(coords~lnCsize*Gonadal_Stage_Index*Sex, data=allometry.geomorph, iter = 9999, RRPP=TRUE)
#summary(LNallometry.lm) # same ANOVA Table
##################################################################################################

################################HACK######################################################
#LNallometry.lm <- procD.lm(coords~lnCsize, data=allometry.geomorph, iter = 9999, RRPP=TRUE)
#summary(LNallometry.lm) # same ANOVA Table

plot(LNallometry.lm, type = "PC", pch = 19, col = "blue")
plotLNallometry.lm <- plot(LNallometry.lm, type = "regression",
     predictor = allometry.geomorph$lnCsize, reg.type = "RegScore",
     pch = 19, col = "green")
variablesforLNallometry <- as.data.frame(cbind(allometry.geomorph$lnCsize,  
                                               plotLNallometry.lm$RegScore,
                                               as.character(allometry.geomorph$Sex),
                                               allometry.geomorph$Gonad_Stage))
colnames(variablesforLNallometry) <- c("lnCsize", "RegScore", "Sex", "GonadStage")
variablesforLNallometry$lnCsize <- as.numeric(as.character(variablesforLNallometry$lnCsize))
variablesforLNallometry$RegScore <- as.numeric(as.character(variablesforLNallometry$RegScore))
ggplot(variablesforLNallometry, aes(y=RegScore, x=lnCsize, color=Sex)) +
  geom_point(size=3) + geom_smooth(method="lm", aes(fill=Sex, colour=Sex)) +
  xlab("ln(Centroid Size)") +
  ylab("Regression Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))



ggplot(variablesforLNallometry, aes(y=RegScore, x=lnCsize, color=GonadStage)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)



#obtain size-adjusted resids
LNallometry.lm <- procD.lm(coords~lnCsize, data=allometry.geomorph, iter = 9999, RRPP=TRUE)
LNallometry.resid <- arrayspecs(LNallometry.lm$residuals, p=dim(procrustes$coords)[1], k=dim(procrustes$coords)[2]) # allometry-adjusted residuals
adj.shape <- LNallometry.resid + array(procrustes$consensus, dim(LNallometry.resid)) # allometry-free shapes

mean.shape.proc.adj<-mshape(adj.shape)


#PLOT GONADAL STAGE WITH SIZE ADJUSTED RESIDUALS
# #group<-allometry.geomorph$Gonadal_Stage_Index
# adjshapeGonads.lm <- procD.lm(adj.shape~Gonadal_Stage_Index*lnCsize, data=allometry.geomorph, iter = 9999, RRPP=TRUE)
# summary(adjshapeGonads.lm) # same ANOVA Table
# plot(adjshapeGonads.lm, type = "PC", pch = 19, col = "blue")
# plotadjshapeGonads.lm <- plot(LNallometry.lm, type = "regression",
#                            predictor = allometry.geomorph$Gonadal_Stage_Index, reg.type = "RegScore",
#                            pch = 19, col = "green")
# variablesforLNallometry$lnCsize <- allometry.geomorph$lnCsize
# variablesforLNallometry$RegScore <- plotLNallometry.lm$RegScore
# variablesforLNallometry$Sex <- allometry.geomorph$Sex
# variablesforLNallometry$GonadStage <- allometry.geomorph$Gonad_Stage
# variablesforLNallometry$AdjRegScore <- plotadjshapeGonads.lm$RegScore
# 
# ggplot(variablesforLNallometry, aes(y=AdjRegScore, x=Gonadal_Stage_Index)) +
#   geom_point()
# ggplot(variablesforLNallometry, aes(y=AdjRegScore, x=lnCsize, color=Sex)) +
#   geom_point() + geom_smooth(method="lm")
# ggplot(variablesforLNallometry, aes(y=AdjRegScore, x=lnCsize, color=GonadStage)) +
#   geom_point() + geom_smooth(method="lm", se=FALSE)

#this is the shape of a fish before adjusting for size and group
plotRefToTarget(mean.shape.proc, procrustes$coords[,,"30 "], links=links, method ="points")

#this is the shape of a fish after adjusting for size and group
plotRefToTarget(mean.shape.proc, adj.shape[,,"30 "], links=links, method ="points")

# PCA of allometry-free shape
pca_adj.shape<-plotTangentSpace(adj.shape) 
summary(pca_adj.shape)

#visualize fish shapes at PC extremes
pdf(file="Rimages/PCA_subtracted_allometry_sex_wireframes.pdf")
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC1max, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC1min, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC2max, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC2min, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC3max, links=links, method ="points", mag=1)
plotRefToTarget(mean.shape.proc, pca_adj.shape$pc.shapes$PC3min, links=links, method ="points", mag=1)
dev.off()

#ggplot version of pca
pc.scores<-as.data.frame(pca_adj.shape$pc.scores)
pdf(file="Rimages/PCA_subtracted_allometry_sex_.pdf")
ggplot(pc.scores, aes(PC1,PC2, color=metadata$Sex))+
   geom_point(size=5)+
   stat_ellipse(type="norm", level=0.5, size=3)+
   theme_bw()+
   theme(axis.text=element_text(size=18), 
         axis.title=element_text(size=18,face="bold"),
         legend.text=element_text(size=18),
         legend.title=element_text(size=18,face="bold"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_blank())+
   guides(color=guide_legend(title="Sex"))+
   geom_hline(aes(yintercept = 0)) +
   geom_vline(aes(xintercept=0))+
   scale_color_brewer(palette="Set1")

ggplot(pc.scores, aes(PC1,PC2, color=metadata$Location))+
  geom_point(size=3)+
  stat_ellipse(type="norm", level=0.5, size=1)+
  theme_bw()+
  theme(axis.text=element_text(size=18), 
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  guides(color=guide_legend(title="Location"))+
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept=0))+
  scale_color_brewer(palette="Set1")

# ggplot(pc.scores, aes(PC1,PC3, color=metadata$Location))+
#   geom_point(size=5)+
#   stat_ellipse(type="norm", level=0.5, size=3)+
#   theme_bw()+
#   theme(axis.text=element_text(size=18), 
#         axis.title=element_text(size=18,face="bold"),
#         legend.text=element_text(size=18),
#         legend.title=element_text(size=18,face="bold"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank())+
#   guides(color=guide_legend(title="Location"))+
#   geom_hline(aes(yintercept = 0)) +
#   geom_vline(aes(xintercept=0))+
#   scale_color_brewer(palette="Set1")
dev.off()
###############################################################################################
#test for effect of site
metadata_noU <- metadata[metadata$Sex != "U",]
invfishingpressure_PC1_noU <- invfishingpressure_PC1[metadata$Sex != "U"]
adj.shape_noU <- adj.shape[,,metadata$Sex != "U"]

Location<-metadata_noU$Location
Sex <- metadata_noU$Sex
adj.shape.lm.location <- procD.lm(adj.shape_noU~Location*Sex, logsz=F, iter=9999, RRPP = TRUE)
summary(adj.shape.lm.location)
plot(adj.shape.lm.location, type = "diagnostics",outliers = TRUE)

# PC plot rotated to major axis of fitted values
plot(adj.shape.lm.location, type = "PC", pch = 19, col = "blue")

# Uses residuals from model to find the commonom regression component
# for a predictor from the model
plot(adj.shape.lm.location, type = "regression", predictor = invfishingpressure_PC1_noU, reg.type = "RegScore",
     pch = 19, col = "green")
# Uses residuals from model to find the projected regression scores
adj.shape.lm.location.RegScore.plot <- plot(adj.shape.lm.location, type = "regression", predictor = invfishingpressure_PC1_noU, reg.type = "RegScore", pch = 21, bg = "black")
#############################################################################


###############################################################################################
#test for effect of fishing pressure
adj.shape.lm <- procD.lm(adj.shape_noU~invfishingpressure_PC1_noU*Sex , logsz=F, iter=9999, RRPP = TRUE)
summary(adj.shape.lm)
plot(adj.shape.lm, type = "diagnostics",outliers = TRUE)

# PC plot rotated to major axis of fitted values
plot(adj.shape.lm, type = "PC", pch = 19, col = "blue")

plot(adj.shape.lm, method="PredLine")

# Uses residuals from model to find the commonom regression component
# for a predictor from the model
plot(adj.shape.lm, type = "regression", predictor = invfishingpressure_PC1_noU, reg.type = "RegScore",
     pch = 19, col = "green")
# Uses residuals from model to find the projected regression scores
adj.shape.lm.RegScore.plot <- plot(adj.shape.lm, type = "regression", predictor = invfishingpressure_PC1_noU, reg.type = "RegScore", pch = 21, bg = "black")
#############################################################################

##############################################################################
#test for effect of fishing pressure, stringent

#CRC
# adj.shape.lm.stringent<-as.data.frame(adj.shape.lm.RegScore.plot$CRC)
# adj.shape.lm.stringent$Location_Code<-metadata$Location_Code
# adj.shape.lm.stringent$invfishingpressure_PC1<-metadata$invfishingpressure_PC1
# 
# model.stringent<-lme(V1 ~ invfishingpressure_PC1,random=~1|Location_Code,data=adj.shape.lm.stringent)
# summary(model.stringent)
# AIC(model.stringent)
# dataf_summary<-summarySE(adj.shape.lm.stringent,measurevar="V1", groupvars=c("invfishingpressure_PC1","Location_Code"), conf.interval=0.99)
# ggplot(dataf_summary, aes(x=invfishingpressure_PC1, y=V1, color=Location_Code)) + 
#   geom_smooth(method=lm, color="black",size=2, se=TRUE, level=0.95) +
#   geom_errorbar(aes(ymin=V1-se, ymax=V1+se), size=2, width=0.5) +
#   geom_point(size=4) +
#   theme(text = element_text(size=18))


#RegScore
adj.shape.lm.stringent.regscore<-as.data.frame(adj.shape.lm.RegScore.plot$RegScore)
colnames(adj.shape.lm.stringent.regscore)[1] <- "RegScore"
adj.shape.lm.stringent.regscore$Location_Code<-metadata_noU$Location_Code
adj.shape.lm.stringent.regscore$invfishingpressure_PC1_noU<-invfishingpressure_PC1_noU
adj.shape.lm.stringent.regscore$Sex <- metadata_noU$Sex
#adj.shape.lm.stringent.regscore$invfishingpressure_PC2<-metadata$invfishingpressure_PC2

#remove Sex = "U" from dataset
#metadata_noU <- subset(metadata, Sex == "F" | Sex == "M")

model.stringent.regscore<-lme(RegScore ~ invfishingpressure_PC1_noU*Sex,random=~1|Location_Code,data=adj.shape.lm.stringent.regscore)
summary(model.stringent.regscore)
anova(model.stringent.regscore)
AIC(model.stringent.regscore)
dataf_summary.regscore<-summarySE(adj.shape.lm.stringent.regscore,measurevar="RegScore", groupvars=c("invfishingpressure_PC1_noU","Location_Code", "Sex"), conf.interval=0.99)
colnames(dataf_summary.regscore)[2] <- "Location"
dataf_summary.regscore$lnRegScore <- log(dataf_summary.regscore$RegScore + (1-min(dataf_summary.regscore$RegScore)))
pdf("Rimages/RegScoreVsFishingPressure.pdf")
ggplot(dataf_summary.regscore, aes(x=invfishingpressure_PC1_noU, y=RegScore, color=Sex, group=Sex)) + 
  geom_smooth(method=lm, color="black",size=3, se=FALSE, level=0.95) +
  geom_errorbar(aes(ymin=RegScore-se, ymax=RegScore+se), size=3, width=0.5) +
  geom_point(size=5) +
  xlab("Fishing Pressure (PC1)") +
  ylab("Regression Score") +
  theme_bw()+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_brewer(palette="Set1")

dev.off()

pdf("Rimages/RegScoreVsFishingPressure_wireframes_eachfish.pdf")
for(i in order(adj.shape.lm.stringent.regscore$RegScore)){
  i
  print(plotRefToTarget(mean.shape.proc, adj.shape[,,i], links=links, method ="points", mag=1))
}
dev.off()


#make list with index pos of first fish for each location
fishlist <- c(min(which(metadata$Location_Code == "BAIS")),
              min(which(metadata$Location_Code == "SIP")),
              min(which(metadata$Location_Code == "AYU")),
              min(which(metadata$Location_Code == "KAL")),
              min(which(metadata$Location_Code == "DUMA")))

pdf("Rimages/RegScoreVsFishingPressure_wireframes_4sites.pdf")
for(i in fishlist){
  a<-t(array(adj.shape.lm$fitted[i,], dim=c(2,24)))
  #plot(a[,2]~a[,1])
  print(plotRefToTarget(mean.shape.proc, a, links=links, method ="points", mag=2))
}
dev.off()

#############################################################################################

#plot lnCsize across locations
ggplot(allometry.geomorph, aes(lnCsize, color=metadata$Location)) +
  geom_point(size=3)

plot(allometry.geomorph$lnCsize, type = "PC", pch = 19, col = "blue")

#write Csize and lnCsize to CSV for inclusion in metadata
#Csize <- (allometry.geomorph$Csize)
#write.csv(Csize, "Csize.csv")


#Subset metadata to include only males and females
#NOTE: may need to go back and re-run metadata_notfixed and metadata without removing gonadal_stage = NA
metadata_notfixed<-read.csv("metadata_ruben.csv")
metadata <- metadata_notfixed[complete.cases(metadata_notfixed$Sex),]
Metadata_SexMF <- metadata[metadata$Sex %in% c("M", "F"), ]

#LWR plot for all locations, male and female
ggplot(data=Metadata_SexMF, aes(y=Weight_g, x=Total_Length_cm, color=Sex)) +
  facet_grid(PC1 ~ Sex) +
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
ggplot(data=Metadata_SexMF, aes(y=B.value, x=PC1, color = Sex, shape=Sex)) +
  geom_point(shape=1, size=15, stroke=5) +
  stat_smooth(method = "lm", aes(fill=Sex, colour=Sex)) +
  xlab("Fishing Pressure (PC1)") +
  ylab("Fish Growth (B-value)") +
  labs(color="Sex") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

####TOTAL WEIGHT V FISHING PRESSURE (PC1)
full_model_continuous <- lmer(Weight_g ~ Sex * PC1 + (1 | Location), data = Metadata_SexMF)
full_model_continuous_lme <- lme(Weight_g ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)
boxCox(lm(Weight_g ~ Sex * PC1, data = Metadata_SexMF))

full_model_continuous_lme_log <- lme(log(Weight_g) ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'PC1')



prediction_data <- expand.grid(Sex = c('F', 'M'), PC1 = seq(min(Metadata_SexMF$PC1), max(Metadata_SexMF$PC1), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)



#Graph of Total Weight v PC1
ggplot(data = Metadata_SexMF, aes(y=Weight_g, x=PC1, 
                       group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Fishing Pressure (PC1)") +
  ylab("Total Weight (g)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###TOTAL LENGTH 
full_model_continuous <- lmer(Total_Length_cm ~ Sex * PC1 + (1 | Location), data = Metadata_SexMF)
full_model_continuous_lme <- lme(Total_Length_cm ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)
boxCox(lm(Total_Length_cm ~ Sex * PC1, data = Metadata_SexMF))

full_model_continuous_lme_log <- lme(log(Total_Length_cm) ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'PC1')

prediction_data <- expand.grid(Sex = c('F', 'M'), PC1 = seq(min(Metadata_SexMF$PC1), max(Metadata_SexMF$PC1), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#Boxplot of Total Length v PC1
ggplot(data = Metadata_SexMF, aes(y=Total_Length_cm, x=PC1, 
                       group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Fishing Pressure (PC1)") +
  ylab("Total Length (cm)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

###TOTAL WEIGHT/TOTAL LENGTH
full_model_continuous <- lmer(TW_TL ~ Sex * PC1 + (1 | Location), data = Metadata_SexMF)
full_model_continuous_lme <- lme(TW_TL ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)
boxCox(lm(TW_TL ~ Sex * PC1, data = Metadata_SexMF))

full_model_continuous_lme_log <- lme(log(TW_TL) ~ Sex * PC1, random = ~1|Location, data = Metadata_SexMF)
diag.plots(full_model_continuous_lme_log, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)

full_model_continuous_lme_log_varIdent <- update(full_model_continuous_lme_log, weights = varIdent(form=~1|Location))
diag.plots(full_model_continuous_lme_log_varIdent, col.nos = c('Sex', 'PC1', "Location"), Metadata_SexMF)


ranova(full_model_continuous)
anova(full_model_continuous_lme_log_varIdent)
summary(full_model_continuous_lme_log_varIdent)
emtrends(full_model_continuous_lme_log_varIdent, pairwise ~ Sex, var = 'PC1')

prediction_data <- expand.grid(Sex = c('F', 'M'), PC1 = seq(min(Metadata_SexMF$PC1), max(Metadata_SexMF$PC1), length.out = 100))
prediction_data$fit <- predict(full_model_continuous_lme_log_varIdent, newdata = prediction_data, level = 0)
prediction_data$fit <- exp(prediction_data$fit)


#Boxplot of Total Weight/Total Length v PC1
ggplot(data = Metadata_SexMF, aes(y=TW_TL, x=PC1, 
                                  group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_line(data = prediction_data, aes(y = fit, group = Sex)) + 
  xlab("Fishing Pressure (PC1)") +
  ylab("Total Weight/Total Length (g/cm)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))


#Boxplot of Csize v PC1
ggplot(data = Metadata_SexMF, aes(y=Csize, x=PC1, 
                                  group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_smooth(method = "lm", fill = NA, aes(group=Metadata_SexMF$Sex)) +  
  xlab("Fishing Pressure (PC1)") +
  ylab("Centroid Size") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

#Boxplot of lnCsize v PC1
ggplot(data = Metadata_SexMF, aes(y=lnCsize, x=PC1, 
                                  group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_smooth(method = "lm", fill = NA, aes(group=Metadata_SexMF$Sex)) +  
  xlab("Fishing Pressure (PC1)") +
  ylab("ln(Centroid Size)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

#Boxplot of TW/Csize v PC1
ggplot(data = Metadata_SexMF, aes(y=TW_Csize, x=PC1, 
                                  group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_smooth(method = "lm", fill = NA, aes(group=Metadata_SexMF$Sex)) +  
  xlab("Fishing Pressure (PC1)") +
  ylab("Total Weight/Centroid Size") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25))

#Boxplot of lnTW/lnCsize v PC1
ggplot(data = Metadata_SexMF, aes(y=lnTW_lnCsize, x=PC1, 
                                  group=interaction(Location,Sex), colour = Sex)) +
  geom_boxplot(size=2) +
  geom_smooth(method = "lm", fill = NA, aes(group=Metadata_SexMF$Sex)) +  
  xlab("Fishing Pressure (PC1)") +
  ylab("ln(Total Weight)/ln(Centroid Size)") +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
