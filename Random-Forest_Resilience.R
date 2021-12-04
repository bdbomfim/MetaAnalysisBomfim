##Response Meta Forest code

library(forestmodel)
library(metaforest)
library(caret)
library(ranger)
library(dmetar)
library(MuMIn)
library(leaps)

###META FOREST Total Litterfall Resilience 1 to 21 months Ambient####

str(tot_lit_amb_1to21_final)#213 obs of 87 variables - including CTE with Hugo data

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

datametaforest_restot<-tot_lit_amb_1to21_final[,c(3,15,23,25,27,29,78:89)]%>% filter(hurrwind!="NA")
str(datametaforest_restot)#213obs 18 variables

# Run model with many trees to check convergence
check_conv_restot <- MetaForest(yi~.,
                         data = datametaforest_restot,
                         study = "Effectsize_ID",
                         whichweights = "random",
                         num.trees = 20000)

plot(check_conv_restot) #model converged at about 10000 trees, which will be used in the next step

#Model with 10000 trees for replication
mf_rep_restot <- MetaForest(yi~.,
                     data = datametaforest_restot,
                     study = "Effectsize_ID",
                     whichweights = "random",
                     num.trees = 10000)
results_res<-summary(mf_rep_restot, digits=2)
results_res

#now apply recursive pre-selection using the preselect function. 

# Run recursive preselection, store results in object 'preselect'
#There are three possible algorithms for variable selection : replicate, recursive, and bootstrap (but latter does not work with this data)

preselected_restot <- preselect(mf_rep_restot,
                         replications = 100,
                         algorithm = "recursive")
plot(preselected_restot)# of variables

preselected2_restot <- preselect(mf_rep_restot,
                          replications = 100,
                          algorithm = "replicate")
plot(preselected2_restot)

#Both generated similar results

#Using preselect_vars, we retain only those moderators for which a 50% percentile interval 
#of the variable importance metrics does not include zero

#Retain only moderators with positive variable importance
#in more than 50% of replications
retain_mods_restot <- preselect_vars(preselected_restot, cutoff = .5)

plot(retain_mods_restot)#number of variables retained

# Set up 3-fold grouped (=clustered) CV

grouped_cv_restot <- trainControl(method = "cv", 
                           index = groupKFold(datametaforest_restot$Effectsize_ID, k = 3))
grouped_boot_restot <- trainControl(method = "cboot", 
                             index = groupKFold(datametaforest_restot$Effectsize_ID, k = 3))

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
#Examine optimal tuning parameters
mf_cv_restot$results[which.min(mf_cv_restot$results$RMSE), ]
#Based on the root mean squared error, the best combination of tuning parameters
#were RANDOM weights
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

#REsidual heterogeneity
Vimp_restot <- data.frame(final_restot$forest$variable.importance)
str(Vimp_restot)
final_restot$forest$variable.importance

Vimp3_restot <- cbind(data.frame(varimp=final_restot$forest$variable.importance,predictors=retain_mods_restot))
Vimp3_restot
str(Vimp3_restot)
head(Vimp3_restot)
Vimp3_restot
(final_restot$forest$variable.importance)

Fig4<-VarImpPlot(final)
Fig4

#### Figure 8a Resilience META FOREST

#Data
var_importance_restot <- data.frame(variable=setdiff(colnames(Vimp3_restot), "Predictors"),
                             importance=as.vector(final_restot$forest$variable.importance))
var_importance_restot
names(datametaforest_restot)
predictor_names_restot<-c("Holdridge zone","Geological group","Parent material P","Parent material","Soil order","Soil P","Wind speed","Wind duration","Time since cyclone","Longitude","Elevation","MAT/MAP","Storm frequency","Time since last storm","Cyclone rainfall")
predictor_names_restot
Vimp3_restot <- cbind(data.frame(varimp=final_restot$forest$variable.importance,predictors=predictor_names_restot))
str(Vimp3_restot)
Vimp3_restot
Vimp3_restot$predictors<- factor(Vimp3_restot$predictors, levels=Vimp3_restot$predictors)

#Figure
Fig8a <- ggplot(Vimp3_restot, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig8a
Fig8a <- Fig8a + geom_bar(col="black",fill="darkgray")+coord_flip()+theme_pubr()
Fig8a
Fig8a <- Fig8a + ylab("Relative importance of predictors") + xlab("")+ theme(
  axis.text.y=element_text(size=26),
  axis.text.x=element_text(size=26),
  axis.title.x=element_text(size=28),
  legend.title=element_blank(),
  legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)+ annotate("text", y = 0.12, x = 1, label = "a", size=12,hjust=0,fontface="bold",colour="black")#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig8a#+annotate INCLUDE ANNOTATION letter and Total litterfall mass

#Saving high res
ggsave(filename = "Fig8a_VarImp_TotLit_Resilience.png",
       plot = Fig8a, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 600)

####META FOREST Leaf fall####

str(leaf_amb_1to21)
Data_leaf_amb_1to21<-leaf_amb_1to21 %>% filter(HURRECON_wind_ms!="NA")
str(Data_leaf_amb_1to21)
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

datametaforest_reslf<-Data_leaf_amb_1to21[,c(3,15,23,25,27,29,79:84,91:96)]%>% filter(hurrwind!="NA")
str(datametaforest_reslf)#198 obs of 18 variables

# Run model with many trees to check convergence
check_conv_reslf <- MetaForest(yi~.,
                            data = datametaforest_reslf,
                            study = "Effectsize_ID",
                            whichweights = "random",
                            num.trees = 20000)

plot(check_conv_reslf) #m odel converged at about 8000 trees, which will be used in the next step

#Model with 10000 trees for replication

mf_rep_reslf<- MetaForest(yi~.,
                       data = datametaforest_reslf,
                       study = "Effectsize_ID",
                       whichweights = "random",
                       num.trees = 8000)
results_reslf<-summary(mf_rep_reslf, digits=2)
results_reslf

preselected_reslf <- preselect(mf_rep_reslf,
                            replications = 100,
                            algorithm = "recursive")
retain_mods_reslf <- preselect_vars(preselected_reslf, cutoff = .5)

grouped_cv_reslf <- trainControl(method = "cv", 
                              index = groupKFold(datametaforest_reslf$Effectsize_ID, k = 3))
grouped_boot_reslf <- trainControl(method = "cboot", 
                                index = groupKFold(datametaforest_reslf$Effectsize_ID, k = 3))

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

# For convenience, extract final model
final_reslf <- mf_cv_reslf$finalModel
final_reslf

r2_oob_reslf <- final_reslf$forest$r.squared
r2_oob_reslf #predictive performance of 48%

#REsidual heterogeneity
Vimp_reslf <- data.frame(final_reslf$forest$variable.importance)
plot(Vimp_reslf)
final_reslf$forest$variable.importance

VarImpPlot(final_reslf)

Vimp_reslf <- cbind(data.frame(varimp=final_reslf$forest$variable.importance,predictors=retain_mods_reslf))
str(Vimp_reslf)
Vimp_reslf

#### Figure 8c Resilience META FOREST Leaf fall####

#Data frame
var_importance_reslf <- data.frame(variable=setdiff(colnames(Vimp_reslf), "Predictors"),
                               importance=as.vector(final_reslf$forest$variable.importance))

names(var_importance_reslf)
var_importance_reslf
predictor_names_reslf<-c("Holdridge zone","Geological group","Parent material P","Parent material","Soil order","Time since cyclone","Longitude","Elevation","MAT/MAP","Soil P","Storm frequency","Time since last storm","Cyclone rainfall","Wind speed","Wind duration")
predictor_names_reslf
Vimp3_reslf <- cbind(data.frame(varimp=final_reslf$forest$variable.importance,predictors=predictor_names_reslf))
str(Vimp3_reslf)

Vimp3_reslf$predictors<- factor(Vimp3_reslf$predictors, levels=Vimp3_reslf$predictors)

#Figure
Fig8c <- ggplot(Vimp3_reslf, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig8c <- Fig8c + geom_bar(col="black",fill="#197D32")+coord_flip()+theme_pubr()
Fig8c
Fig8c <- Fig8c + ylab("Relative importance of predictors") + xlab("")+ theme(
  axis.text.y=element_text(size=26),
  axis.text.x=element_text(size=26),
  axis.title.x=element_text(size=28),
  legend.title=element_blank(),
  legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)+ annotate("text", y = 0.25, x = 1, label = "b", size=12,hjust=0,fontface="bold",colour="black")#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig8c

#Saving in high res
ggsave(filename = "Fig8c_VarImp_Leaf_Resilience.png",
       plot = Fig8c, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 800)

Fig8ac<-Fig8a+Fig8c
Fig8ac
ggsave(filename = "Fig8ac_VarImp_TotLeaf_Resilience.png",
       plot = Fig8ac, width = 24, height = 12, units = 'cm',
       scale = 2, dpi = 800)

###END###
