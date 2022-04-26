#Resilience GAMM predictions to create synthetic case

#Packages
library(mgcv)
library(nlme)
library(gamm4)
library(ggplot2)
library(dplyr)
library(tidymv)
library(metafor)
library(tidyverse)
library(dplyr)
library(patchwork)

##Upload data####
metadat<-read.csv(file.choose())#Litterfall_Mass in data: folder
str(metadat)#2370 obs of 62 variables

##Data wrangling####
#subset including observations 1 month and on post-cyclone#
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
rec

tot_lit_all<-data_esall %>% filter(Fraction=="TotLitfall")

#Filtering data to include Ambient and CTE
tot_lit_amb_cte<-data_esall %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")

tot_lit_amb_1to21_final<-tot_lit_amb_cte %>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")

z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
tot_lit_amb_1to21_final$soilP<-z.trans(tot_lit_amb_1to21_final$Other_soil_P)
tot_lit_amb_1to21_final$hurrwind<-z.trans(tot_lit_amb_1to21_final$HURRECON_wind_ms)
tot_lit_amb_1to21_final$windur<-z.trans(tot_lit_amb_1to21_final$Gale_wind_duration_minutes)
tot_lit_amb_1to21_final$tsd<-z.trans(tot_lit_amb_1to21_final$TSD_months)

tot_meta<- rma.mv(yi,vi,random = list(~1|Site,~1|DisturbanceName),
                  tdist = TRUE,
                  data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta)
#these are the sigma2 values used to calculate weight2
tot_meta$sigma2 # both values are included in the next step

#Calculate weights for GAMMs: adding the values obtained by running tot_meta_amb$sigma2
tot_lit_amb_1to21_final$weight2<-(1/(tot_lit_amb_1to21_final$vi+0.02803322+0.17168580))
tot_lit_amb_1to21_final$weight2

#Best model Table 3 Model 1a####
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP,tsd, k=20)                       ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=T)
summary(gamm_2y_mixed_1$gam)#R2=0.25

#Summary of each predictor variable
summary(tot_lit_amb_1to21_final$soilP)
summary(tot_lit_amb_1to21_final$tsd)

# Preparing new data frames based on existing soil P and tsd values
#transformed data
#For time since disturbance at 12 months
newdat2.2a <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.1853929,0.1853929,0.1853929,0.1853929))
newdat2.2a
#At 1 month (min)
newdat2.3 <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.7348944,-0.7348944,-0.7348944,-0.7348944))
newdat2.3
#At 5 months (1st quartile)
newdat2.3a <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.4002445,-0.4002445,-0.4002445,-0.4002445))
newdat2.3a
#At 9 months (median)
newdat2.3b <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.06559,-0.06559,-0.06559,-0.06559))
newdat2.3b
#At 15 months (3rd quartile)
newdat2.3c <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.43638,0.43638,0.43638,0.43638))
newdat2.3c

# Getting predicted outcomes for new data
# These include the splines but ignore other REs
preds <- predict(gamm_2y$gam, type="terms")[,"s(soilP)"]
plot(gamm_2y$gam)

predictions2.2a = predict(gamm_2y_mixed_1$gam, newdata=newdat2.2a, se.fit = TRUE)
predictions2.3 = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3, se.fit = TRUE)
predictions2.3a = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3a, se.fit = TRUE)
predictions2.3b = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3b, se.fit = TRUE)
predictions2.3c = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3c, se.fit = TRUE)

# Combining new data and predictions ####
newdat2.2a = cbind(newdat2.2a, predictions2.2a) 
newdat2.2a$month<- c(12,12,12,12)

newdat2.3 = cbind(newdat2.3, predictions2.3)
newdat2.3$month<- c(1,1,1,1)

newdat2.3a = cbind(newdat2.3a, predictions2.3a) # 5 months
newdat2.3a$month<- c(5,5,5,5)

newdat2.3b = cbind(newdat2.3b, predictions2.3b) # 9 months
newdat2.3b$month<- c(9,9,9,9)

newdat2.3c = cbind(newdat2.3c, predictions2.3c) # 15 months
newdat2.3c$month<- c(15,15,15,15)

# Calculating CIs
newdat2.2a <- within(newdat2.2a, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3 <- within(newdat2.3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3a <- within(newdat2.3a, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3b <- within(newdat2.3b, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3c <- within(newdat2.3c, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

#Combining dataframes by row to prepare single figure

#Combining dataframes by row
gam_pred_data<-rbind(data.frame(SoilP=newdat2.2a$soilP,Tsd=newdat2.2a$tsd, PredRes=newdat2.2a$fit,CIlow=newdat2.2a$lower,Ciup=newdat2.2a$upper,Month=newdat2.2a$month,MonthBin="12"),
                     data.frame(SoilP=newdat2.3$soilP,Tsd=newdat2.3$tsd, PredRes=newdat2.3$fit,CIlow=newdat2.3$lower,Ciup=newdat2.3$upper,Month=newdat2.3$month,MonthBin="1"),
                     data.frame(SoilP=newdat2.3a$soilP,Tsd=newdat2.3a$tsd, PredRes=newdat2.3a$fit,CIlow=newdat2.3a$lower,Ciup=newdat2.3a$upper,Month=newdat2.3a$month,MonthBin="5"),
                     data.frame(SoilP=newdat2.3b$soilP,Tsd=newdat2.3b$tsd, PredRes=newdat2.3b$fit,CIlow=newdat2.3b$lower,Ciup=newdat2.3b$upper,Month=newdat2.3b$month,MonthBin="9"),
                     data.frame(SoilP=newdat2.3c$soilP,Tsd=newdat2.3c$tsd, PredRes=newdat2.3c$fit,CIlow=newdat2.3c$lower,Ciup=newdat2.3c$upper,Month=newdat2.3c$month,MonthBin="15"))
gam_pred_data

# Plotting the predicted resilience as a function of soil P

#Figure all month bins
egplot_all <- ggplot(data = gam_pred_data, aes(x=SoilP, y=PredRes, group=MonthBin)) + 
   geom_point(aes(group=MonthBin,col=MonthBin))+geom_line(aes(group=MonthBin,col=MonthBin)) + labs(y="Predicted resilience", x="Standardized soil P")+
    theme(axis.title.x=element_text(vjust = 0.5,size=18),
        axis.text=element_text(vjust = 1,size=16),
        axis.title.y=element_text(vjust = 1,size=18))+ggtitle("Table 3 model 1a - wind duration not included")
egplot_all

##Figure S8a####
ggsave(filename = "FigS8a.png",
       plot = egplot_all, width = 8, height = 6, units = 'cm',
       scale = 2, dpi = 1000)

## Preparing figure by using predictive model including gale wind duration##

#Table 3 Model 2a (ambient + CTE)####
gamm_2y_mixed_1e <- gamm4(yi ~ s(soilP, tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1e$gam)#R2 = 0.4
summary(gamm_2y_mixed_1e$mer)

#Summary of wind data
summary(tot_lit_amb_1to21_final$Gale_wind_duration_minutes)
summary(tot_lit_amb_1to21_final$windur)
summary(tot_lit_amb_1to21_final$tsd)

# Data frames including soil P, tsd and wind dur range of values ####

## 1st quartile wind duration
#At 12 months
newdat2.2aw1 <- data.frame(windur = -0.5360, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.1853929,0.1853929,0.1853929,0.1853929))
newdat2.2aw1
#At 1 month (min)
newdat2.3w1 <- data.frame(windur = -0.5360, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.7348944,-0.7348944,-0.7348944,-0.7348944))
newdat2.3w1
#At 5 months (1st quartile)
newdat2.3aw1 <- data.frame(windur = -0.5360, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.4002445,-0.4002445,-0.4002445,-0.4002445))
newdat2.3aw1
#At 9 months (median)
newdat2.3bw1 <- data.frame(windur = -0.5360, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.06559,-0.06559,-0.06559,-0.06559))
newdat2.3bw1
#At 15 months (3rd quartile)
newdat2.3cw1 <- data.frame(windur = -0.5360, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.43638,0.43638,0.43638,0.43638))
newdat2.3cw1

## 3rd quartile wind duration
#At 12 months
newdat2.2aw3 <- data.frame(windur = 0.6407, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.1853929,0.1853929,0.1853929,0.1853929))
newdat2.2aw3
#At 1 month (min)
newdat2.3w3 <- data.frame(windur = 0.6407, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.7348944,-0.7348944,-0.7348944,-0.7348944))
newdat2.3w3
#At 5 months (1st quartile)
newdat2.3aw3 <- data.frame(windur = 0.6407, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.4002445,-0.4002445,-0.4002445,-0.4002445))
newdat2.3aw3
#At 9 months (median)
newdat2.3bw3 <- data.frame(windur = 0.6407, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.06559,-0.06559,-0.06559,-0.06559))
newdat2.3bw3
#At 15 months (3rd quartile)
newdat2.3cw3 <- data.frame(windur = 0.6407, soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(0.43638,0.43638,0.43638,0.43638))
newdat2.3cw3

# Predictions with wind duration ####

# 1st quartile wind duration
predictions2.2aw1 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.2aw1, se.fit = TRUE)
predictions2.2aw1
predictions2.3w1 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3w1, se.fit = TRUE)
predictions2.3w1
predictions2.3aw1 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3aw1, se.fit = TRUE)
predictions2.3aw1
predictions2.3bw1 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3bw1, se.fit = TRUE)
predictions2.3bw1
predictions2.3cw1 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3cw1, se.fit = TRUE)
predictions2.3cw1

# 3rd quartile wind duration
predictions2.2aw3 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.2aw3, se.fit = TRUE)
predictions2.2aw3
predictions2.3w3 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3w3, se.fit = TRUE)
predictions2.3w3
predictions2.3aw3 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3aw3, se.fit = TRUE)
predictions2.3aw3
predictions2.3bw3 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3bw3, se.fit = TRUE)
predictions2.3bw3
predictions2.3cw3 = predict(gamm_2y_mixed_1e$gam, newdata=newdat2.3cw3, se.fit = TRUE)
predictions2.3cw3

# Combining new data and predictions including wind duration####

#1st quartile wind duration
newdat2.2aw1 = cbind(newdat2.2aw1, predictions2.2aw1) 
newdat2.2aw1$month<- c(12,12,12,12)

newdat2.3w1 = cbind(newdat2.3w1, predictions2.3w1)
newdat2.3w1$month<- c(1,1,1,1)

newdat2.3aw1 = cbind(newdat2.3aw1, predictions2.3aw1) # 5 months
newdat2.3aw1$month<- c(5,5,5,5)

newdat2.3bw1= cbind(newdat2.3bw1, predictions2.3bw1) # 9 months
newdat2.3bw1$month<- c(9,9,9,9)

newdat2.3cw1 = cbind(newdat2.3cw1, predictions2.3cw1) # 15 months
newdat2.3cw1$month<- c(15,15,15,15)

#3rd quartile wind duration
newdat2.2aw3 = cbind(newdat2.2aw3, predictions2.2aw3) 
newdat2.2aw3$month<- c(12,12,12,12)

newdat2.3w3 = cbind(newdat2.3w3, predictions2.3w3)
newdat2.3w3$month<- c(1,1,1,1)

newdat2.3aw3 = cbind(newdat2.3aw3, predictions2.3aw3) # 5 months
newdat2.3aw3$month<- c(5,5,5,5)

newdat2.3bw3= cbind(newdat2.3bw3, predictions2.3bw3) # 9 months
newdat2.3bw3$month<- c(9,9,9,9)

newdat2.3cw3 = cbind(newdat2.3cw3, predictions2.3cw3) # 15 months
newdat2.3cw3$month<- c(15,15,15,15)

# Calculating CIs with wind duration ####

#1st quartile wind duration
newdat2.2aw1 <- within(newdat2.2aw1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3w1 <- within(newdat2.3w1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3aw1 <- within(newdat2.3aw1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3bw1 <- within(newdat2.3bw1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3cw1 <- within(newdat2.3cw1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

#3rd quartile wind duration
newdat2.2aw3 <- within(newdat2.2aw3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3w3 <- within(newdat2.3w3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3aw3 <- within(newdat2.3aw3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3bw3 <- within(newdat2.3bw3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.3cw3 <- within(newdat2.3cw3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

#Combining dataframes by row to prepare single figure ####

#Combining dataframes by row
#1st quartile
newdat2.2aw1
gam_pred_dataw1<-rbind(data.frame(Windur=newdat2.2aw1$windur, SoilP=newdat2.2aw1$soilP,Tsd=newdat2.2aw1$tsd, PredRes=newdat2.2aw1$fit,CIlow=newdat2.2aw1$lower,Ciup=newdat2.2aw1$upper,Month=newdat2.2aw1$month,MonthBin="12"),
                     data.frame(Windur=newdat2.3w1$windur, SoilP=newdat2.3w1$soilP,Tsd=newdat2.3w1$tsd, PredRes=newdat2.3w1$fit,CIlow=newdat2.3w1$lower,Ciup=newdat2.3w1$upper,Month=newdat2.3w1$month,MonthBin="1"),
                     data.frame(Windur=newdat2.3aw1$windur, SoilP=newdat2.3aw1$soilP,Tsd=newdat2.3aw1$tsd, PredRes=newdat2.3aw1$fit,CIlow=newdat2.3aw1$lower,Ciup=newdat2.3aw1$upper,Month=newdat2.3aw1$month,MonthBin="5"),
                     data.frame(Windur=newdat2.3bw1$windur, SoilP=newdat2.3bw1$soilP,Tsd=newdat2.3bw1$tsd, PredRes=newdat2.3bw1$fit,CIlow=newdat2.3bw1$lower,Ciup=newdat2.3bw1$upper,Month=newdat2.3bw1$month,MonthBin="9"),
                     data.frame(Windur=newdat2.3cw1$windur, SoilP=newdat2.3cw1$soilP,Tsd=newdat2.3cw1$tsd, PredRes=newdat2.3cw1$fit,CIlow=newdat2.3cw1$lower,Ciup=newdat2.3cw1$upper,Month=newdat2.3cw1$month,MonthBin="15"))
gam_pred_dataw1

# Plotting the predicted resilience as a function of soil P
summary(gam_pred_dataw1)

#Figure all month bins
egplot_allw1 <- ggplot(data = gam_pred_dataw1, aes(x=SoilP, y=PredRes, group=MonthBin)) + 
  geom_point(aes(group=MonthBin,col=MonthBin))+geom_line(aes(group=MonthBin,col=MonthBin)) + labs(y="Predicted resilience", x="Standardized soil P")+
  theme(axis.title.x=element_text(vjust = 0.5,size=18),
        axis.text=element_text(vjust = 1,size=16),
        axis.title.y=element_text(vjust = 1,size=18))+ ggtitle("Table 3 model 2a - 1st quartile wind duration")
egplot_allw1

#Combining dataframes by row
#3rd quartile
newdat2.2aw3
gam_pred_dataw3<-rbind(data.frame(Windur=newdat2.2aw3$windur, SoilP=newdat2.2aw3$soilP,Tsd=newdat2.2aw3$tsd, PredRes=newdat2.2aw3$fit,CIlow=newdat2.2aw3$lower,Ciup=newdat2.2aw3$upper,Month=newdat2.2aw3$month,MonthBin="12"),
                       data.frame(Windur=newdat2.3w3$windur, SoilP=newdat2.3w3$soilP,Tsd=newdat2.3w3$tsd, PredRes=newdat2.3w3$fit,CIlow=newdat2.3w3$lower,Ciup=newdat2.3w3$upper,Month=newdat2.3w3$month,MonthBin="1"),
                       data.frame(Windur=newdat2.3aw3$windur, SoilP=newdat2.3aw3$soilP,Tsd=newdat2.3aw3$tsd, PredRes=newdat2.3aw3$fit,CIlow=newdat2.3aw3$lower,Ciup=newdat2.3aw3$upper,Month=newdat2.3aw3$month,MonthBin="5"),
                       data.frame(Windur=newdat2.3bw3$windur, SoilP=newdat2.3bw3$soilP,Tsd=newdat2.3bw3$tsd, PredRes=newdat2.3bw3$fit,CIlow=newdat2.3bw3$lower,Ciup=newdat2.3bw3$upper,Month=newdat2.3bw3$month,MonthBin="9"),
                       data.frame(Windur=newdat2.3cw3$windur, SoilP=newdat2.3cw3$soilP,Tsd=newdat2.3cw3$tsd, PredRes=newdat2.3cw3$fit,CIlow=newdat2.3cw3$lower,Ciup=newdat2.3cw3$upper,Month=newdat2.3cw3$month,MonthBin="15"))
gam_pred_dataw3

# Plotting the predicted resilience as a function of soil P
summary(gam_pred_dataw3)

#Figure all month bins
egplot_allw3 <- ggplot(data = gam_pred_dataw3, aes(x=SoilP, y=PredRes, group=MonthBin)) + 
  geom_point(aes(group=MonthBin,col=MonthBin))+geom_line(aes(group=MonthBin,col=MonthBin)) + labs(y="Predicted resilience", x="Standardized soil P")+
  theme(axis.title.x=element_text(vjust = 0.5,size=18),
        axis.text=element_text(vjust = 1,size=16),
        axis.title.y=element_text(vjust = 1,size=18))+ ggtitle("Table 3 model 2a - 3rd quartile wind duration")
egplot_allw3

# Combining all three plots

FinalFigS8 <- egplot_all+egplot_allw1+egplot_allw3+plot_layout(ncol=1)
FinalFigS8

FinalFigS8<-egplot_all/(egplot_allw1+egplot_allw3)+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 18,face="bold"),plot.tag.position =c(0.05,1))
FinalFigS8

##Final Figure S8####
ggsave(filename = "FinalFigS8.png",
       plot = FinalFigS8, width = 16, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

#Sample if trying to exclude random effects

head(
  cbind(
    gam1 = predict_gamm(gamm_2y_mixed_1, re_form = NA),
    gam2 = predict_gamm(gamm_2y_mixed_1,
                        exclude = c("s(Subject)", "s(Days,Subject)")
    )
  )
)

## END ###

