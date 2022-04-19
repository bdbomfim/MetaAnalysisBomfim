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

#Example

summary(tot_lit_amb_1to21_final$soilP)
summary(tot_lit_amb_1to21_final$Other_soil_P)
summary(tot_lit_amb_1to21_final$tsd)
summary(tot_lit_amb_1to21_final$TSD_months)
rec5 <- tot_lit_amb_1to21_final %>% filter(TSD_months == 5)
head(rec5)

# Generating some new data for which you'd like predictions:
#synthetic data
newdat <- data.frame(soilP = runif(100), tsd = runif(100))

#data based on existing soil P and tsd values
#untransformed data
newdat2 <- data.frame(soilP = c(130, 210, 275, 330), tsd = c(9,9,9,9))
newdat2.1 <- data.frame(soilP = c(130, 210, 275, 330), tsd = c(1,5,9,15))
newdat2.1a<- data.frame(soilP = c(130, 210, 275, 330), tsd = c(12,12,12,12))

#transformed data
newdat2.2 <- data.frame(soilP = c(-0.48852, -0.21887, 0.05414, 0.18560), tsd = c(-0.06559,-0.06559,-0.06559,-0.06559))
#At 12 months
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
predictions = predict(gamm_2y_mixed_1$gam, newdata=newdat, se.fit = TRUE)
predictions
predictions2 = predict(gamm_2y_mixed_1$gam, newdata=newdat2, se.fit = TRUE)
predictions2
predictions2.1 = predict(gamm_2y_mixed_1$gam, newdata=newdat2.1, se.fit = TRUE)
predictions2.1
predictions2.1a = predict(gamm_2y_mixed_1$gam, newdata=newdat2.1a, se.fit = TRUE)
predictions2.1a
predictions2.2 = predict(gamm_2y_mixed_1$gam, newdata=newdat2.2, se.fit = TRUE)
predictions2.2
predictions2.2a = predict(gamm_2y_mixed_1$gam, newdata=newdat2.2a, se.fit = TRUE)
predictions2.2a
newdat2.2a
predictions2.3 = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3, se.fit = TRUE)
predictions2.3
predictions2.3a = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3a, se.fit = TRUE)
predictions2.3a
predictions2.3b = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3b, se.fit = TRUE)
predictions2.3b
predictions2.3c = predict(gamm_2y_mixed_1$gam, newdata=newdat2.3c, se.fit = TRUE)
predictions2.3c

#predicions2 = predict_gam(gamm_2y_mixed_1$gam, values = list(tsd = c(-0.73489, -0.06559, 0.43638)),newdata=newdat2, se.fit = TRUE)
#predictions2

# Consolidating new data and predictions

newdat2.2a = cbind(newdat2.2a, predictions2.2a) 
newdat2.2a$month<- c(12,12,12,12)
newdat2.2a
#newdat2.2a = newdat2.2a %>% 
  #rename(soilP12 = soilP, tsd12 = tsd, fit12 = fit, se.fit12 = se.fit)
newdat2.3 = cbind(newdat2.3, predictions2.3)
newdat2.3 #1 month
newdat2.3$month<- c(1,1,1,1)
#newdat2.3 = newdat2.3 %>% 
  #rename(soilP1 = soilP, tsd1 = tsd)
newdat2.3a = cbind(newdat2.3a, predictions2.3a) # 5 months
newdat2.3a
newdat2.3a$month<- c(5,5,5,5)
#newdat2.3a = newdat2.3a %>% 
  #rename(soilP5 = soilP, tsd5 = tsd)
newdat2.3b = cbind(newdat2.3b, predictions2.3b) # 9 months
newdat2.3b$month<- c(9,9,9,9)
#newdat2.3b = newdat2.3b %>% 
  #rename(soilP9 = soilP, tsd9 = tsd)
newdat2.3b
newdat2.3c = cbind(newdat2.3c, predictions2.3c) # 15 months
newdat2.3c$month<- c(15,15,15,15)
#newdat2.3c = newdat2.3c %>% 
  #rename(soilP15 = soilP, tsd15 = tsd)
newdat2.3c

#Combining dataframes by row
gam_pred_data<-rbind(data.frame(SoilP=newdat2.2a$soilP,Tsd=newdat2.2a$tsd, PredRes=newdat2.2a$fit,CIlow=newdat2.2a$lower,Ciup=newdat2.2a$upper,Month=newdat2.2a$month,MonthBin="12"),
                    data.frame(SoilP=newdat2.3$soilP,Tsd=newdat2.3$tsd, PredRes=newdat2.3$fit,CIlow=newdat2.3$lower,Ciup=newdat2.3$upper,Month=newdat2.3$month,MonthBin="1"),
                    data.frame(SoilP=newdat2.3a$soilP,Tsd=newdat2.3a$tsd, PredRes=newdat2.3a$fit,CIlow=newdat2.3a$lower,Ciup=newdat2.3a$upper,Month=newdat2.3a$month,MonthBin="5"),
                    data.frame(SoilP=newdat2.3b$soilP,Tsd=newdat2.3b$tsd, PredRes=newdat2.3b$fit,CIlow=newdat2.3b$lower,Ciup=newdat2.3b$upper,Month=newdat2.3b$month,MonthBin="9"),
                    data.frame(SoilP=newdat2.3c$soilP,Tsd=newdat2.3c$tsd, PredRes=newdat2.3c$fit,CIlow=newdat2.3c$lower,Ciup=newdat2.3c$upper,Month=newdat2.3c$month,MonthBin="15"))
gam_pred_data

# Calculating CIs
newdat <- within(newdat, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2 <- within(newdat2, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.1 <- within(newdat2.1, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.1a <- within(newdat2.1a, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
newdat2.2 <- within(newdat2.2, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})
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

#Combining dataframes by row to prepare single figure ## CONTINUE HERE!

# Plot the predicted resilience as a function of soil P

summary(gam_pred_data)

#Figure all months
egplot_all <- ggplot(data = gam_pred_data, aes(x=SoilP, y=PredRes, group=MonthBin)) + 
   geom_point(aes(group=MonthBin,col=MonthBin))+geom_line(aes(group=MonthBin,col=MonthBin)) + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
    theme(axis.title.x=element_text(vjust = 0.5,size=20),
        axis.text=element_text(vjust = 1,size=18),
        axis.title.y=element_text(vjust = 1,size=20))
egplot_all

##New Figure9####
ggsave(filename = "Fig9_new.png",
       plot = Fig_9, width = 10, height = 10, units = 'cm',
       scale = 2, dpi = 1000)

#Figure at 12 months
newdat2.2a
egplot2.2a <- ggplot(data = newdat2.2a, aes(x=soilP, y=fit)) + 
  geom_smooth_ci() + geom_point() + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=22),
        axis.text =element_text(vjust = 1,size=20),
        axis.title.y =element_text(vjust = 1,size=22))
egplot2.2a
#Figure at 1 month
egplot2.3 <- ggplot(newdat2.3, aes(x=soilP, y=fit)) + 
  geom_smooth_ci() + geom_point() + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=22),
        axis.text =element_text(vjust = 1,size=20),
        axis.title.y =element_text(vjust = 1,size=22))
egplot2.3
#Figure at 5 month
egplot2.3a <- ggplot(newdat2.3a, aes(x=soilP, y=fit)) + 
  geom_smooth_ci() + geom_point() + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=22),
        axis.text =element_text(vjust = 1,size=20),
        axis.title.y =element_text(vjust = 1,size=22))
egplot2.3a
#Figure at 9 month
egplot2.3b <- ggplot(newdat2.3b, aes(x=soilP, y=fit)) + 
  geom_smooth_ci() + geom_point() + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=22),
        axis.text =element_text(vjust = 1,size=20),
        axis.title.y =element_text(vjust = 1,size=22))
egplot2.3b
#Figure at 15 month
egplot2.3c <- ggplot(newdat2.3c, aes(x=soilP, y=fit)) + 
  geom_smooth_ci() + geom_point() + labs(y="Predicted resilience", x="Soil P (ln mg P/kg)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=22),
        axis.text =element_text(vjust = 1,size=20),
        axis.title.y =element_text(vjust = 1,size=22))
egplot2.3c




#Calculating resilience values
((exp(0.047369)-1)*100)
0.1944441 - 0.1470751
0.047369/0.1944441

#Excluding random effects

head(
  cbind(
    gam1 = predict_gamm(gamm_2y_mixed_1, re_form = NA),
    gam2 = predict_gamm(gamm_2y_mixed_1,
                        exclude = c("s(Subject)", "s(Days,Subject)")
    )
  )
)
#To plot the smooths across a few values of a continuous predictor, 
#we can use the values argument in predict_gam().
#predict_gam(gamm_2y_mixed_1$gam, values = list(tsd = c(1, 5, 12))) %>%
  #ggplot(data=newdat2, aes(x=soilP, y=fit)) #+geom_smooth_ci(tsd)

