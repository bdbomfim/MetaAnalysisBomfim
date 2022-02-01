##Random Forest Code for Litterfall Resilience

#Required packages
library(forestmodel)
library(metaforest)
library(caret)
library(ranger)
library(dmetar)
library(MuMIn)
library(leaps)
library(tidyverse)
library(ggpubr)

####Uploading and saving data in metadat dataframe####
metadat<-read.csv(file.choose())#Litterfall_Mass.csv in 'data' folder
str(metadat)
#transforming variables to numeric
metadat$TSD_months=as.numeric(metadat$TSD_months)
metadat$HURRECON_wind_ms=as.numeric(metadat$HURRECON_wind_ms)
metadat$Gale_wind_duration_minutes=as.numeric(metadat$Gale_wind_duration_minutes)

##Resilience Data wrangling####
#subset including observations 1 month and on post-cyclone#
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
#subset inculding 1 to 36 months post-cyclone and excluding duplicated studies and obs. in Bisley
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)
#subset including ambient conditions only
res_amb<-res_all %>% filter(Treatment=="Ambient")#|Treatment=="TrimDeb")
str(res_amb)#945 obs with 417 from Puerto Rico (44%) #1125obs with TrimDeb and 597 from Puerto Rico(53%)
summary(res_amb$Treatment)
#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36####
data_esall <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_amb, measure = "ROM")
str(data_esall)#945 obs. including all litterfall mass fractions

#Creating new columns yi_new as the result of yi divided by duration####
names(data_esall)
data_esall$yi_new <- data_esall$yi / data_esall$TSD_months
#checking new column yi_new
summary(data_esall$yi_new)
summary(data_esall$yi)
#Same for the variance vi_new
data_esall$vi_new <- data_esall$vi / data_esall$TSD_months
#checking new column vi_new
summary(data_esall$vi_new)
summary(data_esall$vi)

#Filtering data frame to include only total litterfall mass flux between 1 and 21 months post-cyclone
tot_lit_amb<-data_esall %>% filter(Fraction=="TotLitfall")#%>% filter(Treatment=="Ambient")
#Data 1 to 21 months - Ambient only
tot_lit_amb_1to21_final<-tot_lit_amb %>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
summary(tot_lit_amb_1to21_final$Treatment)#192 obs of 66 variables

###Random Forest Total litterfall mass flux Resilience 1 to 21 months Ambient####
#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
tot_lit_amb_1to21_final$tsd<-z.trans(tot_lit_amb_1to21_final$TSD_months)
tot_lit_amb_1to21_final$long<-z.trans(tot_lit_amb_1to21_final$Longitude)
tot_lit_amb_1to21_final$elev<-z.trans(tot_lit_amb_1to21_final$Elevation_m)
tot_lit_amb_1to21_final$mat_map<-z.trans(tot_lit_amb_1to21_final$MAT_MAP_x100)
tot_lit_amb_1to21_final$soilP<-z.trans(tot_lit_amb_1to21_final$Other_soil_P)
tot_lit_amb_1to21_final$stormfreq<-z.trans(tot_lit_amb_1to21_final$StormFrequencyNorm)
tot_lit_amb_1to21_final$timesincestorm<-z.trans(tot_lit_amb_1to21_final$YearsSinceLastStorm)
tot_lit_amb_1to21_final$distrain<-z.trans(tot_lit_amb_1to21_final$Disturb_Rainfall_mm)
tot_lit_amb_1to21_final$hurrwind<-z.trans(tot_lit_amb_1to21_final$HURRECON_wind_ms)
tot_lit_amb_1to21_final$windur<-z.trans(tot_lit_amb_1to21_final$Gale_wind_duration_minutes)

names(tot_lit_amb_1to21_final)
summary(tot_lit_amb_1to21_final$Treatment)

#Data frame for random forest analysis with effect sizes, variance and desirable moderators
datametaforest_restot<-tot_lit_amb_1to21_final[,c(3,8,15,23,25,27,29,63:74)]%>% filter(hurrwind!="NA")
str(datametaforest_restot)#192 obs 17 variables, 213 with CTE
# Rename column where name is "yi_new" and "vi_new" to allow the MetaForest function to run
names(datametaforest_restot)[names(datametaforest_restot) == "yi_new"] <- "yi"
names(datametaforest_restot)[names(datametaforest_restot) == "vi_new"] <- "vi"

#Run model with many trees to check convergence
#For any random forest model, it is important to check whether the model converges.
check_conv_restot <- MetaForest(yi~.,
                         data = datametaforest_restot,
                         study = "Effectsize_ID",
                         whichweights = "random",
                         num.trees = 20000)
plot(check_conv_restot) #model converged at about 10000 trees, which will be used in the next step

#Model with 5000trees for replication
mf_rep_restot <- MetaForest(yi~.,
                     data = datametaforest_restot,
                     study = "Effectsize_ID",
                     whichweights = "random",
                     num.trees = 7500)
results_res<-summary(mf_rep_restot, digits=2)
results_res

#Now apply recursive pre-selection using the preselect function.

# Run recursive preselection, store results in object 'preselect'
#There are three possible algorithms for variable selection : replicate, recursive, and bootstrap (but latter does not work with this data)
#applying the recursive algorithm
preselected_restot <- preselect(mf_rep_restot,
                         replications = 100,
                         algorithm = "recursive")
plot(preselected_restot)# of variables

#Applying the replicate algorithm
preselected_restot1 <- preselect(mf_rep_restot,
                                replications = 100,
                                algorithm = "replicate")
plot(preselected_restot1)# of variables

#Applying the bootstrap algorithm to circumvent data dependency
preselected2_restot <- preselect(mf_rep_restot,
                          replications = 100,
                          algorithm = "bootstrap")#super weird
plot(preselected2_restot)
#Both generated similar results, so sticking with the 'recursive' algorithm for variable selection

#Using preselect_vars, we retain only those moderators for which a 50% percentile interval 
#of the variable importance metrics does not include zero
#Retain only moderators with positive variable importance in more than 50% of replications
retain_mods_restot <- preselect_vars(preselected_restot, cutoff = .5)

#Tuning the model with package caret
# Set up 5-fold grouped (=clustered) cross-validation (CV)
grouped_cv_restot <- trainControl(method = "cv", 
                           index = groupKFold(datametaforest_restot$Effectsize_ID, k = 5))#5-fold clustered cross-validation
grouped_boot_restot <- trainControl(method = "cboot", 
                             index = groupKFold(datametaforest_restot$Effectsize_ID, k = 5))

# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid_restot <- expand.grid(whichweights = c("random", "fixed", "unif"),
                           mtry = 2:6,
                           min.node.size = 2:6)

# X should contain only retained moderators, clustering variable, and vi
X_restot <- datametaforest_restot[, c("Effectsize_ID", "vi", retain_mods_restot)]

# Train the model
mf_cv_restot <- train(y = datametaforest_restot$yi,
               x = X_restot,
               study = "Effectsize_ID", # Name of the clustering variable
               method = ModelInfo_mf(), 
               trControl = grouped_cv_restot,
               tuneGrid = tuning_grid_restot,
               num.trees = 10000)

# Examine optimal tuning parameters
mf_cv_restot$results[which.min(mf_cv_restot$results$RMSE), ]
#Based on the root mean squared error, the best combination of tuning parameters
#were (random) weights - this changes if group of moderators change
#The object returned by train already contains the final model, 
#estimated with the best combination of tuning parameters

#Inspecting the results
r2_cv_restot<-mf_cv_restot$results$Rsquared[which.min(mf_cv_restot$results$RMSE)]
r2_cv_restot #0.39

# Extract the estimate of predictive performance R^2_{oob} from the final model
#This is an estimate of how much variance the model would explain in a new data set 
r2_oob_restot <- final_restot$forest$r.squared
r2_oob_restot #predictive performance of 31%

# For convenience, extract final model
final_restot <- mf_cv_restot$finalModel
final_restot
#Checking convergence
plot(final_restot)

#We can conclude that the model has converged, and has a positive estimate of explained variance in new data. 
#Now, we proceed to interpreting the model findings. 

# Plot variable importance
#MetaForest’s variable importance metrics reflect each moderator’s contribution 
#to the predictive power of the final model across all linear-, non-linear-, and interaction effects
VarImpPlot(final_restot)

# Sort the variable names by importance, so that the
# partial dependence plots will be ranked by importance
ordered_vars_restot <- names(final_restot$forest$variable.importance)[
  order(final_restot$forest$variable.importance, decreasing = TRUE)]
ordered_vars_restot

##Preparing Figs 8a and 8b

#Preparein dataframe with moderators
Vimp_restot <- data.frame(final_restot$forest$variable.importance)
str(Vimp_restot)
final_restot$forest$variable.importance

Vimp3_restot <- cbind(data.frame(varimp=final_restot$forest$variable.importance,predictors=retain_mods_restot))
Vimp3_restot
names(Vimp3_restot)
(final_restot$forest$variable.importance)


#### Figure 8a Resilience META FOREST

#Data
var_importance_restot <- data.frame(variable=setdiff(colnames(Vimp3_restot), "Predictors"),
                             importance=as.vector(final_restot$forest$variable.importance))
var_importance_restot
names(datametaforest_restot)
predictor_names_restot<-c("Holdridge zone","Geological group","Parent material P","Parent material","Soil order","Time since cyclone","Longitude","Elevation","MAT/MAP","Soil P","Storm frequency","Time since last storm", "Cyclone rainfall", "Wind speed","Wind duration")
predictor_names_restot
Vimp3_restot <- cbind(data.frame(varimp=final_restot$forest$variable.importance,predictors=predictor_names_restot))
str(Vimp3_restot)
Vimp3_restot
Vimp3_restot$predictors<- factor(Vimp3_restot$predictors, levels=Vimp3_restot$predictors)

#Figure8a####
Fig8a <- ggplot(Vimp3_restot, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig8a <- Fig8a + geom_bar(col="black",fill="darkgray")+coord_flip()+theme_pubr()
Fig8a <- Fig8a + ylab("Moderator (permutation) importance") + xlab("")+ theme(
  axis.text.y=element_text(size=24),
  axis.text.x=element_text(size=24),
  axis.title.x=element_text(size=24),
  legend.title=element_blank(),
  legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)+ annotate("text", y = 0.11, x = 1, label = "a", size=12,hjust=0,fontface="bold",colour="black")#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig8a

#Saving high res
ggsave(filename = "Fig8a_Random-forest_Resilience-TL-region.png",
       plot = Fig8a, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 600)

####META FOREST Leaf fall####
leaf_all<-data_esall %>% filter(Fraction=="Leaf fall")
str(leaf_all)#358 obs

#1 to 36 months - Ambient only
leaf_amb<-leaf_all%>% filter(Treatment=="Ambient")#|Treatment=="TrimDeb")
str(leaf_amb)#231 obs no CTE, 267 obs with CTE
summary(leaf_amb$Treatment)
leaf_amb$Case_study= paste(leaf_amb$Site, leaf_amb$DisturbanceName, sep="|")

##Subseting the data 1 - 21 months by deleting cyclones Ivor, Jova, Patricia and Gilbert and Study #4
leaf_amb_1to21<- leaf_amb  %>% filter (TSD_months<22) %>% filter(Study_ID!="4")%>% filter(DisturbanceName!="Keith")%>% filter(DisturbanceName!="Ivor")%>% filter(DisturbanceName!="Jova")%>% filter(DisturbanceName!="Patricia")%>% filter(DisturbanceName!="Gilbert")
str(leaf_amb_1to21)#193 obs | 172 obs Ambient only
#Filtering to include only data with HURRECON wind speed
Data_leaf_amb_1to21<-leaf_amb_1to21 %>% filter(HURRECON_wind_ms!="NA")
str(Data_leaf_amb_1to21)#193 obs | 172 obs Ambient only

#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
Data_leaf_amb_1to21$tsd<-z.trans(Data_leaf_amb_1to21$TSD_months)
Data_leaf_amb_1to21$long<-z.trans(Data_leaf_amb_1to21$Longitude)
Data_leaf_amb_1to21$elev<-z.trans(Data_leaf_amb_1to21$Elevation_m)
Data_leaf_amb_1to21$mat_map<-z.trans(Data_leaf_amb_1to21$MAT_MAP_x100)
Data_leaf_amb_1to21$soilP<-z.trans(Data_leaf_amb_1to21$Other_soil_P)
Data_leaf_amb_1to21$stormfreq<-z.trans(Data_leaf_amb_1to21$StormFrequencyNorm)
Data_leaf_amb_1to21$timesincestorm<-z.trans(Data_leaf_amb_1to21$YearsSinceLastStorm)
Data_leaf_amb_1to21$distrain<-z.trans(Data_leaf_amb_1to21$Disturb_Rainfall_mm)
Data_leaf_amb_1to21$hurrwind<-z.trans(Data_leaf_amb_1to21$HURRECON_wind_ms)
Data_leaf_amb_1to21$windur<-z.trans(Data_leaf_amb_1to21$Gale_wind_duration_minutes)

names(Data_leaf_amb_1to21)
summary(Data_leaf_amb_1to21$yi_new)
summary(Data_leaf_amb_1to21$yi)

datametaforest_reslf<-Data_leaf_amb_1to21[,c(3,8,15,23,25,27,29,63:64,66:75)]%>% filter(hurrwind!="NA")
str(datametaforest_reslf)#172 obs (amb), 193 (with CTE) of 17 variables
# Rename column where name is "yi_new" and "vi_new" to allow the MetaForest function to run
names(datametaforest_reslf)[names(datametaforest_reslf) == "yi_new"] <- "yi"
names(datametaforest_reslf)[names(datametaforest_reslf) == "vi_new"] <- "vi"
names(datametaforest_reslf)

# Run model with many trees to check convergence
check_conv_reslf <- MetaForest(yi~.,
                            data = datametaforest_reslf,
                            study = "Effectsize_ID",
                            whichweights = "random",
                            num.trees = 20000)

plot(check_conv_reslf) #model converged at about 8000 trees, which will be used in the next step

#Model with 10000 trees for replication
mf_rep_reslf<- MetaForest(yi~.,
                       data = datametaforest_reslf,
                       study = "Effectsize_ID",
                       whichweights = "random",
                       num.trees = 8000)
results_reslf<-summary(mf_rep_reslf, digits=2)
results_reslf
#preselecting variables
preselected_reslf <- preselect(mf_rep_reslf,
                            replications = 100,
                            algorithm = "recursive")
retain_mods_reslf <- preselect_vars(preselected_reslf, cutoff = .5)
grouped_cv_reslf <- trainControl(method = "cv", index = groupKFold(datametaforest_reslf$Effectsize_ID, k = 5))
grouped_boot_reslf <- trainControl(method = "cboot", index = groupKFold(datametaforest_reslf$Effectsize_ID, k = 5))

# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid_reslf <- expand.grid(whichweights = c("random", "fixed", "unif"),
                              mtry = 2:6,
                              min.node.size = 2:6)
# X should contain only retained moderators, clustering variable, and vi
X2_reslf <- datametaforest_reslf[, c("Effectsize_ID", "vi", retain_mods_reslf)]
# Train the model
mf_cv_reslf <- train(y = datametaforest_reslf$yi,
                  x = X2_reslf,
                  study = "Effectsize_ID", # Name of the clustering variable
                  method = ModelInfo_mf(), 
                  trControl = grouped_cv_reslf,
                  tuneGrid = tuning_grid_reslf,
                  num.trees = 8000)
#Examine optimal tuning parameters
mf_cv_reslf$results[which.min(mf_cv_reslf$results$RMSE), ]
r2_cv_reslf<-mf_cv_reslf$results$Rsquared[which.min(mf_cv_reslf$results$RMSE)]
r2_cv_reslf#0.54

#Predictive performance
r2_oob_reslf <- final_reslf$forest$r.squared
r2_oob_reslf #predictive performance of 48%

# For convenience, extract final model
final_reslf <- mf_cv_reslf$finalModel
final_reslf

#Residual heterogeneity
Vimp_reslf <- data.frame(final_reslf$forest$variable.importance)
plot(Vimp_reslf)
final_reslf$forest$variable.importance
#Variable importance plot
VarImpPlot(final_reslf)
#Data frame with moderator names and importance metric
Vimp_reslf <- cbind(data.frame(varimp=final_reslf$forest$variable.importance,predictors=retain_mods_reslf))
str(Vimp_reslf)
Vimp_reslf

#### Figure8c Random Forest Resilience Leaf litterfall####

#Prepare dataframe with moderators
Vimp_reslf <- data.frame(final_reslf$forest$variable.importance)
str(Vimp_reslf)
final_reslf$forest$variable.importance

Vimp_reslf <- cbind(data.frame(varimp=final_reslf$forest$variable.importance,predictors=retain_mods_reslf))
Vimp_reslf
str(Vimp_reslf)
head(Vimp_reslf)
Vimp_reslf
(final_reslf$forest$variable.importance)

#### Figure 8a Resilience META FOREST

#Data
#Data frame
var_importance_reslf <- data.frame(variable=setdiff(colnames(Vimp_reslf), "Predictors"),
                               importance=as.vector(final_reslf$forest$variable.importance))

names(var_importance_reslf)
#Assigning new names
predictor_names_reslf<-c("Holdridge zone","Geological group","Parent material P","Parent material","Soil order","Time since cyclone","Longitude","Elevation","MAT/MAP","Soil P","Storm frequency","Time since last storm","Cyclone rainfall","Wind speed","Wind duration")
predictor_names_reslf
Vimp3_reslf <- cbind(data.frame(varimp=final_reslf$forest$variable.importance,predictors=predictor_names_reslf))
Vimp3_reslf$predictors<- factor(Vimp3_reslf$predictors, levels=Vimp3_reslf$predictors)

#Figure
Fig8c <- ggplot(Vimp3_reslf, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig8c <- Fig8c + geom_bar(col="black",fill="#197D32")+coord_flip()+theme_pubr()
Fig8c <- Fig8c + ylab("Moderator (permutation) importance") + xlab("")+ theme(
  axis.text.y=element_text(size=26), axis.text.x=element_text(size=26),
  axis.title.x=element_text(size=28), legend.title=element_blank(), legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)+ annotate("text", y = 0.29, x = 1, label = "b", size=12,hjust=0,fontface="bold",colour="black")#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig8c

#Saving in high res
ggsave(filename = "Fig8c_VarImp_Leaf_Resilience.final-noCTE.png",
       plot = Fig8c, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 800)

##Figure8ab####
Fig8ac<-Fig8a+Fig8c
Fig8ac
ggsave(filename = "Fig8ab_VarImp_TotLeaf_Resilience_fig8nocte.png",
       plot = Fig8ac, width = 26, height = 12, units = 'cm',
       scale = 2, dpi = 1000)

###END###
