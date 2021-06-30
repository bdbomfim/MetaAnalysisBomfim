##Response Meta Forest code

library(forestmodel)
library(metaforest)
library(caret)
library(ranger)
library(dmetar)
library(MuMIn)
library(rJava)
library(leaps)

###META FOREST Total Litterfall####
str(data_es0ia)

#Standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ia$long<-z.trans(data_es0ia$Longitude)
data_es0ia$elev<-z.trans(data_es0ia$Elevation_m)
data_es0ia$mat_map<-z.trans(data_es0ia$MAT_MAP_x100)
data_es0ia$soilP<-z.trans(data_es0ia$Other_soil_P)
data_es0ia$stormfreq<-z.trans(data_es0ia$StormFrequencyNorm)
data_es0ia$timesincestorm<-z.trans(data_es0ia$YearsSinceLastStorm)
data_es0ia$distrain<-z.trans(data_es0ia$Disturb_Rainfall_mm)
data_es0ia$hurrwind<-z.trans(data_es0ia$HURRECON_wind_ms)
data_es0ia$windur<-z.trans(data_es0ia$Gale_wind_duration_minutes)

names(data_es0ia)

#Data frame
datametaforest<-data_es0ia[,c(3,15,23,25,27,29,79:89)]%>% filter(hurrwind!="NA")
str(datametaforest)#46obs 17 variables

# Run model with many trees to check convergence
#MetaForest uses a weighted random forest to explore heterogeneity in meta-analytic data
check_conv <- MetaForest(yi~.,
                         data = datametaforest,
                         study = "Effectsize_ID",
                         whichweights = "random",
                         num.trees = 20000)
plot(check_conv) #model converged at about 7500 trees, which will be used in the next step

#Model with 9000 trees for replication
mf_rep <- MetaForest(yi~.,
                     data = datametaforest,
                     study = "Effectsize_ID",
                     whichweights = "random",#using random effects, as in classic meta-analysis
                     num.trees = 7500)
#Printing the results of the analysis
results<-summary(mf_rep, digits=2)
results
VarImpPlot(mf_rep)#can take a look at the importance value at this point, too

#Preselecting variables for MetaForest analysis
#Applying recursive pre-selection using the preselect function. 
# Run recursive preselection, store results in object 'preselect'
#There are three possible algorithms for variable selection : replicate, recursive, and bootstrap (but latter does not work with this data)
preselected <- preselect(mf_rep,
                         replications = 100,
                         algorithm = "recursive")
plot(preselected)# of variables

#Alternative 'replicate' algorithm
preselected2 <- preselect(mf_rep,
                          replications = 100,
                          algorithm = "replicate")
plot(preselected2)
#Both generated similar results

#Using preselect_vars to retain only those moderators for which a 50% percentile interval 
#of the variable importance metrics does not include zero
#Retain only moderators with positive variable importance in more than 50% of replications
retain_mods <- preselect_vars(preselected, cutoff = .5)

# Set up 3-fold grouped (=clustered) CV
grouped_cv <- trainControl(method = "cv", 
                           index = groupKFold(datametaforest$Effectsize_ID, k = 3))
grouped_boot <- trainControl(method = "cboot", 
                             index = groupKFold(datametaforest$Effectsize_ID, k = 3))

# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid <- expand.grid(whichweights = c("random", "fixed", "unif"),
                           mtry = 2:6,
                           min.node.size = 2:6)

# X should contain only retained moderators, clustering variable, and vi
X <- datametaforest[, c("Effectsize_ID", "vi", retain_mods)]

# Train the random forest model
mf_cv <- train(y = datametaforest$yi,
               x = X,
               study = "Effectsize_ID", # Name of the clustering variable
               method = ModelInfo_mf(), 
               trControl = grouped_cv,
               tuneGrid = tuning_grid,
               num.trees = 7500)

#Examine optimal tuning parameters
mf_cv$results[which.min(mf_cv$results$RMSE), ]
#Based on the root mean squared error, the best combination of tuning parameters
#were UNIF weights

#R2 cv
r2_cv<-mf_cv$results$Rsquared[which.min(mf_cv$results$RMSE)]
r2_cv#0.62

#The object returned by train already contains the final model, 
#estimated with the best combination of tuning parameters

#Inspecting the results
# Extracting the final model
final <- mf_cv$finalModel
final

# Extracting the estimate of predictive performance R^2_{oob} from the final model
#This is an estimate of how much variance the model would explain in a new data set 
r2_oob <- final$forest$r.squared
r2_oob #predictive performance of 51.4%

#Residual heterogeneity
Vimp <- data.frame(final$forest$variable.importance)
Vimp

# Plot convergence
plot(final)

#We can conclude that the model has converged, and has a positive estimate of explained variance in new data. 
#Proceed to interpreting the model findings. 

#Plotting the variable importance, and partial dependence plots.
#MetaForest’s variable importance metrics reflect each moderator’s contribution 
#to the predictive power of the final model across all linear-, non-linear-, and interaction effects
VarImpPlot(final)

# Sort the variable names by importance
ordered_vars <- names(final$forest$variable.importance)[
  order(final$forest$variable.importance, decreasing = TRUE)]
ordered_vars

##Preparing Figs 4a and 4b

Vimp <- data.frame(final$forest$variable.importance)
str(Vimp)
final$forest$variable.importance

Vimp3 <- cbind(data.frame(varimp=final$forest$variable.importance,predictors=retain_mods))
Vimp3
str(Vimp3)
head(Vimp3)
Vimp3

#### Figure 4a Relative importance of predictors - Response of Total Litterfall####

#Data frame
var_importance <- data.frame(variable=setdiff(colnames(Vimp3), "Predictors"),
                             importance=as.vector(final$forest$variable.importance))
var_importance
names(datametaforest)
predictor_names<-c("Holdridge zone","Geological group","Parent material P","Parent material","Soil order","Longitude","Elevation","MAT/MAP","Soil P","Storm frequency","Time since last storm","Cyclone rainfall","Wind speed","Wind duration")
predictor_names
Vimp3 <- cbind(data.frame(varimp=final$forest$variable.importance,predictors=predictor_names))
str(Vimp3)

Vimp3$predictors<- factor(Vimp3$predictors, levels=Vimp3$predictors)

#Figure
Fig4a <- ggplot(Vimp3, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig4a <- Fig4a + geom_bar(col="black",fill="darkgray")+coord_flip()+theme_pubr()
Fig4a <- Fig4a + ylab("Relative importance of predictors") + xlab("")+ theme(
  axis.text.y=element_text(size=26),
  axis.text.x=element_text(size=26),
  axis.title.x=element_text(size=28),
  legend.title=element_blank(),
  legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)#+ annotate("text", y = 1.1, x = 1, label = "a", size=12,hjust=0,fontface="bold",colour="black")#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig4a

#Saving in high res
ggsave(filename = "Fig4a_RF-Response-Tot.png",
       plot = Fig4a, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 1000)

#Correlation plot of final variables in "final2"

#Data
names(datametaforest)
CPFig4b_cor<-datametaforest[,c(2:6,9:17)]
names(CPFig4b_cor)## Final!
names(CPFig4b_cor)
names(CPFig4b_cor)[names(CPFig4b_cor) == "Holdridge_ID"] <- "Holdridge zone"
names(CPFig4b_cor)[names(CPFig4b_cor) == "long"] <- "Longitude"
names(CPFig4b_cor)[names(CPFig4b_cor) == "elev"] <- "Elevation"
names(CPFig4b_cor)[names(CPFig4b_cor) == "mat_map"] <- "MAT/MAP"
names(CPFig4b_cor)[names(CPFig4b_cor) == "soilP"] <- "Soil P"
names(CPFig4b_cor)[names(CPFig4b_cor) == "stormfreq"] <- "Storm frequency"
names(CPFig4b_cor)[names(CPFig4b_cor) == "timesincestorm"] <- "Time since last storm"
names(CPFig4b_cor)[names(CPFig4b_cor) == "distrain"] <- "Cyclone rainfall"
names(CPFig4b_cor)[names(CPFig4b_cor) == "hurrwind"] <- "Wind speed"
names(CPFig4b_cor)[names(CPFig4b_cor) == "windur"] <- "Wind duration"
names(CPFig4b_cor)[names(CPFig4b_cor) == "Rocktype_ID"] <- "Geological group"
names(CPFig4b_cor)[names(CPFig4b_cor) == "RockP_ID"] <- "Parent material P"
names(CPFig4b_cor)[names(CPFig4b_cor) == "Par_Mat_ID"] <- "Parent material"
names(CPFig4b_cor)[names(CPFig4b_cor) == "Soil_ID"] <- "Soil order"

#Calculating correlations and p values
corrf4b <- round(cor(CPFig4b_cor,method="pearson"),2)
p.matf4b <- cor_pmat(CPFig4b_cor)

#Figure
Fig4b_2<-ggcorrplot(corrf4b, hc.order = TRUE, type = "lower",hc.method = "ward.D2",sig.level = 0.05,
                    outline.col = "white", p.mat = p.matf4b,method="square",ggtheme=ggplot2::theme_classic(),show.legend=TRUE, 
                    legend.title="Pearson's r", lab=TRUE, lab_size=6, tl.cex=28,insig="blank",
                    colors = c("#ABA0A0", "white", "#ffa600",pch.cex=20,nbreaks = 8,legend.text.cex=26))+font("legend.text",size=18)+font("legend.title", size=22)#+theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),axis.text.y = element_text(margin=margin(0,-2,0,0)))
Fig4b_2

#Saving figure in high res
ggsave(filename = "Fig4b_Final.png",
       plot = Fig4b_2, width = 16, height = 18, units = 'cm',
       scale = 2, dpi = 1000)

Fig4a+Fig4b

####META FOREST Leaf fall####

#Data
str(data_es0ilf)
Data_es0ilf<-data_es0ilf %>% filter(HURRECON_wind_ms!="NA")
str(Data_es0ilf)
#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
Data_es0ilf$long<-z.trans(Data_es0ilf$Longitude)
Data_es0ilf$elev<-z.trans(Data_es0ilf$Elevation_m)
Data_es0ilf$mat_map<-z.trans(Data_es0ilf$MAT_MAP_x100)
Data_es0ilf$soilP<-z.trans(Data_es0ilf$Other_soil_P)
Data_es0ilf$stormfreq<-z.trans(Data_es0ilf$StormFrequencyNorm)
Data_es0ilf$timesincestorm<-z.trans(Data_es0ilf$YearsSinceLastStorm)
Data_es0ilf$distrain<-z.trans(Data_es0ilf$Disturb_Rainfall_mm)
Data_es0ilf$hurrwind<-z.trans(Data_es0ilf$HURRECON_wind_ms)
Data_es0ilf$windur<-z.trans(Data_es0ilf$Gale_wind_duration_minutes)
names(Data_es0ilf)

#Final data for meta forest analysis
datametaforestlf<-Data_es0ilf[,c(3,15,23,25,27,79:89)]%>% filter(hurrwind!="NA")
str(datametaforestlf)#30 obs and 16 variables

# Run model with many trees to check convergence
check_conv2lf <- MetaForest(yi~.,
                            data = datametaforestlf,
                            study = "Effectsize_ID",
                            whichweights = "random",
                            num.trees = 20000)

plot(check_conv2lf) #model converged at about 7500 trees, which will be used in the next step

#Model with 7500 trees for replication
mf_rep2lf<- MetaForest(yi~.,
                       data = datametaforestlf,
                       study = "Effectsize_ID",
                       whichweights = "random",
                       num.trees = 7500)

preselected2lf <- preselect(mf_rep2lf,replications = 100,algorithm = "recursive")
plot(preselected2lf)

retain_mods2lf <- preselect_vars(preselected2lf, cutoff = .5)

#Control the computational nuances of the train function
grouped_cv2lf <- trainControl(method = "cv", #resampling method
                              index = groupKFold(datametaforestlf$Effectsize_ID, k = 3))
grouped_boot2lf <- trainControl(method = "cboot", 
                                index = groupKFold(datametaforestlf$Effectsize_ID, k = 3))

# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid2lf <- expand.grid(whichweights = c("random", "fixed", "unif"),
                              mtry = 2:6,
                              min.node.size = 2:6)

# X should contain only retained moderators, clustering variable, and vi
X2lf <- datametaforestlf[, c("Effectsize_ID", "vi", retain_mods2lf)]

# Train the model
mf_cv2lf <- train(y = datametaforestlf$yi,
                  x = X2lf,
                  study = "Effectsize_ID", # Name of the clustering variable
                  method = ModelInfo_mf(), 
                  trControl = grouped_cv2lf,
                  tuneGrid = tuning_grid2lf,
                  num.trees = 7500)

#Examine optimal tuning parameters
mf_cv2lf$results[which.min(mf_cv2lf$results$RMSE), ]

r2_cvlf<-mf_cv2lf$results$Rsquared[which.min(mf_cv2lf$results$RMSE)]
r2_cvlf

# For convenience, extract final model
final2lf <- mf_cv2lf$finalModel
final2lf

r2_oob2lf <- final2lf$forest$r.squared
r2_oob2lf #predictive performance of 61%

#Plotting relative importance of predictors on package default
VarImpPlot(final2lf)

#### Figure 4c RESPONSE META FOREST Leaf fall####

#Data frame variables
Vimp2lf <- data.frame(final2lf$forest$variable.importance)
plot(Vimp2lf)
final2lf$forest$variable.importance

Vimp3lf <- cbind(data.frame(varimp=final2lf$forest$variable.importance,predictors=retain_mods2lf))
str(Vimp3lf)
Vimp3lf

var_importancelf <- data.frame(variable=setdiff(colnames(Vimp3lf), "Predictors"),
                               importance=as.vector(final2lf$forest$variable.importance))

names(var_importancelf)
var_importancelf

predictor_nameslf<-c("Holdridge zone","Geological group","Parent material P","Parent material","Longitude","Elevation","MAT/MAP","Soil P","Storm frequency","Time since last storm","Cyclone rainfall","Wind speed","Wind duration")
predictor_nameslf
Vimp3lf <- cbind(data.frame(varimp=final2lf$forest$variable.importance,predictors=predictor_nameslf))
str(Vimp3lf)

Vimp3lf$predictors<- factor(Vimp3lf$predictors, levels=Vimp3lf$predictors)

#Figure 4c Variable Importance
Fig4c <- ggplot(Vimp3lf, aes(x=reorder(predictors,varimp), weight=varimp, fill=varimp))
Fig4c <- Fig4c + geom_bar(col="black",fill="#197D32")+coord_flip()+theme_pubr()
Fig4c
Fig4c <- Fig4c + ylab("Relative importance of predictors") + xlab("")+ theme(
  axis.text.y=element_text(size=26),
  axis.text.x=element_text(size=26),
  axis.title.x=element_text(size=28),
  legend.title=element_blank(),
  legend.text=element_blank(),
  legend.key = element_blank())+guides(fill=FALSE)#+annotate INCLUDE ANNOTATION letter and Total litterfall mass
Fig4c

#Saving in high res
ggsave(filename = "Fig4c_RF_Leaf_Response.png",
       plot = Fig4c, width = 14, height = 12, units = 'cm',
       scale = 2, dpi = 1000)

#Correlation plot of final variables in "final2"

names(datametaforestlf)

CPFig4d_cor<-datametaforestlf[,c(2:5,8:16)]
names(CPFig4d_cor)
#Updating names
names(CPFig4d_cor)[names(CPFig4d_cor) == "Holdridge_ID"] <- "Holdridge zone"
names(CPFig4d_cor)[names(CPFig4d_cor) == "long"] <- "Longitude"
names(CPFig4d_cor)[names(CPFig4d_cor) == "elev"] <- "Elevation"
names(CPFig4d_cor)[names(CPFig4d_cor) == "mat_map"] <- "MAT/MAP"
names(CPFig4d_cor)[names(CPFig4d_cor) == "soilP"] <- "Soil P"
names(CPFig4d_cor)[names(CPFig4d_cor) == "stormfreq"] <- "Storm frequency"
names(CPFig4d_cor)[names(CPFig4d_cor) == "timesincestorm"] <- "Time since last storm"
names(CPFig4d_cor)[names(CPFig4d_cor) == "distrain"] <- "Cyclone rainfall"
names(CPFig4d_cor)[names(CPFig4d_cor) == "hurrwind"] <- "Wind speed"
names(CPFig4d_cor)[names(CPFig4d_cor) == "windur"] <- "Wind duration"
names(CPFig4d_cor)[names(CPFig4d_cor) == "Rocktype_ID"] <- "Geological group"
names(CPFig4d_cor)[names(CPFig4d_cor) == "RockP_ID"] <- "Parent material P"
names(CPFig4d_cor)[names(CPFig4d_cor) == "Par_Mat_ID"] <- "Parent material"
#names(CPFig4d_cor)[names(CPFig4d_cor) == "Soil_ID"] <- "Soil order"

names(CPFig4d_cor)

#Calculating correlation coefficients and p values
corrf4d <- round(cor(CPFig4d_cor,method="pearson"),2)
p.matf4d <- cor_pmat(CPFig4d_cor)

#Figure
Fig4d_2<-ggcorrplot(corrf4d, hc.order = TRUE, type = "lower",hc.method = "ward.D2",sig.level = 0.05,
                    outline.col = "white", p.mat = p.matf4d,method="square",ggtheme=ggplot2::theme_classic(),show.legend=TRUE, 
                    legend.title="Pearson's r", lab=TRUE, lab_size=6, tl.cex=28,insig="blank",
                    colors = c("#46A332", "white", "#ffa600",pch.cex=20,nbreaks = 8,legend.text.cex=26))+font("legend.text",size=18)+font("legend.title", size=22)#+theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),axis.text.y = element_text(margin=margin(0,-2,0,0)))
Fig4d_2

#Saving in high res
ggsave(filename = "Fig4d_Final.png",
       plot = Fig4d_2, width = 16, height = 18, units = 'cm',
       scale = 2, dpi = 1000)

##END####
