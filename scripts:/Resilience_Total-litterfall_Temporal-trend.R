##Post-cyclone temporal dynamics of total litterfall resilience##
#GAMMs##

##Packages
library(performance)
library(lme4)
library(lmerTest)
library(metafor)
library(forestmodel)
library(ggplot2)
library(ggcorrplot)
library(DescTools)
library(nlme)
library(gapminder)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(magrittr)
library(mgcv)
library(nlme)
library(dplyr)
library(gamm4) 
library(AICcmodavg)
library(lattice)
library(RColorBrewer)

##Upload data####
metadat<-read.csv(file.choose())#Litterfall_Mass data file
str(metadat)#2370 obs of 62 variables

##Data wrangling####
#subset including observations 1 month and on post-cyclone#
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
rec

#subset inculding 1 to 36 months post-cyclone and excluding duplicated studies and obs. in Bisley
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36####
data_esall <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_all, measure = "ROM")
str(data_esall)#1527 obs. including all litterfall mass fractions

##Total Litterfall mass flux resilience####
#subset of total litterfall 1 to 36 months
tot_lit_all<-data_esall %>% filter(Fraction=="TotLitfall")
str(tot_lit_all)#425 obs
summary(tot_lit_all$Treatment)

#Filtering data to include Ambient and CTE
tot_lit_amb_cte<-data_esall %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")
str(tot_lit_amb_cte)#285 obs Amb + CTE
tot_lit_amb<-data_esall %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient")#|Treatment=="TrimDeb")
str(tot_lit_amb)#249 obs Amb

##Analyzing resilience for 1 to 21 months, wherein most data is concentrated
#Data 1 to 21 months - Ambient + CTE
tot_lit_amb_1to21_final<-tot_lit_amb_cte %>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_amb_1to21_final)#213
summary(tot_lit_amb_1to21_final$Treatment)
#creating case study column
tot_lit_amb_1to21_final$Case_study= paste(tot_lit_amb_1to21_final$Site, tot_lit_amb_1to21_final$DisturbanceName,sep=" | ")
unique(levels(as.factor(tot_lit_amb_1to21_final$Case_study)))

#Transforming the numerical moderators####
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
tot_lit_amb_1to21_final$soilP<-z.trans(tot_lit_amb_1to21_final$Other_soil_P)
tot_lit_amb_1to21_final$hurrwind<-z.trans(tot_lit_amb_1to21_final$HURRECON_wind_ms)
tot_lit_amb_1to21_final$windur<-z.trans(tot_lit_amb_1to21_final$Gale_wind_duration_minutes)
tot_lit_amb_1to21_final$tsd<-z.trans(tot_lit_amb_1to21_final$TSD_months)


#Excluding CTE for sensitivity analysis
tot_lit_amb_1to21_amb<-tot_lit_amb_1to21_final %>% filter(Treatment!="TrimDeb")
summary(tot_lit_amb_1to21_amb$Treatment)

##Calculating Weights to later use in the GAMMs####

#Fitting multilevel meta-analysis models to obtain the weights included in the final GAMMs
#cross random effects for Site and Cyclone
tot_meta<- rma.mv(yi,vi,random = list(~1|Site,~1|DisturbanceName),
                          tdist = TRUE,
                           data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta)
#these are the sigma2 values used to calculate weight2
tot_meta$sigma2

#Calculate weights for GAMMs: adding the values obtained by running tot_meta_amb$sigma2
tot_lit_amb_1to21_final$weight2<-(1/(tot_lit_amb_1to21_final$vi+0.02803322+0.17168580))
tot_lit_amb_1to21_final$weight2

#other possibility to calculate weights
weight<-weights(tot_meta,type="matrix")
tot_lit_amb_1to21_final$weight<-colSums(weight)/sum(weight)

#cyclone nested within region as random effects
tot_meta_b<- rma.mv(yi,vi,random = ~1|Region/DisturbanceName,
                  tdist = TRUE,
                  data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_b)
tot_meta_b$sigma2
#weights
weight_b<-weights(tot_meta_b,type="matrix")
tot_lit_amb_1to21_final$weight_b<-colSums(weight_b)/sum(weight_b)
tot_lit_amb_1to21_final$weight_b
tot_lit_amb_1to21_final$weight_b2<-(1/(tot_lit_amb_1to21_final$vi+0.00000006+0.084))

#cyclone as random effect
tot_meta_c<- rma.mv(yi,vi,random = ~1|DisturbanceName,
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_c)
#weights
weight_c<-weights(tot_meta_c,type="matrix")
tot_lit_amb_1to21_final$weight_c<-colSums(weight_c)/sum(weight_c)
tot_lit_amb_1to21_final$weight_c

#region as random effect
tot_meta_d<- rma.mv(yi,vi,random = ~(1|Region),
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_d)
#weights
weight_d<-weights(tot_meta_d,type="matrix")
tot_lit_amb_1to21_final$weight_d<-colSums(weight_d)/sum(weight_d)
tot_lit_amb_1to21_final$weight_d

#site as random effect
tot_meta_e<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_e)
#weights
weight_e<-weights(tot_meta_e,type="matrix")
tot_lit_amb_1to21_final$weight_e<-colSums(weight_e)/sum(weight_e)
tot_lit_amb_1to21_final$weight_e


#Fitting GAMMs using weights generated from rma.mv function####

#Site as random effect with site-based weights
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP, tsd)                              ,weights=(1/tot_lit_amb_1to21_final$vi), random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1.2 <- gamm4(yi ~ s(soilP, tsd)                              ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )

#using AICc to check best approach to weights
mods_1 <- list(mixed_1 = gamm_2y_mixed_1$mer, mixed_1.2 = gamm_2y_mixed_1.2$mer)
atab_1 <- aictab(mods_1)
atab_1#MUCH BETTER TO USE THE SECOND APPROACH TO WEIGHTS

#Fitting several GAMMs 
gamm_2y_mixed_2 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=(1/tot_lit_amb_1to21_final$vi), random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2.2 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_b2, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2.1 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_b, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_3 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_5 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_e, random = ~(1|Site), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1a <- gamm4(yi ~ s(soilP) +s(tsd)                           ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2a <- gamm4(yi ~ s(Other_soil_P,by=TSD_months)              ,weights=tot_lit_amb_1to21_final$weight_b, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_3a <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4a <- gamm4(yi ~ s(Other_soil_P, TSD_months,k=20)            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_5a <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_e, random = ~(1|Site), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4b <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1c <- gamm4(yi ~ soilP + s(tsd)                                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1d <- gamm4(yi ~ soilP + tsd                                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mix_1e <- gamm4(yi ~ s(soilP, tsd,by=windur)                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )

#Comparing models using AICc
mods_1to21 <- list(mixed_1 = gamm_2y_mixed_1$mer, mixed_1a = gamm_2y_mixed_1a$mer,mixed_1b = gamm_2y_mixed_1b$mer,
                   mixed_1c = gamm_2y_mixed_1c$mer,mixed_1d = gamm_2y_mixed_1d$mer,mix_1e = gamm_2y_mix_1e$mer)
atab_1to21 <- aictab(mods_1to21)
atab_1to21

#Table 3 GAMMs####
summary(tot_lit_amb_1to21_final$Country)
#Model 1a - ambient + CTE####
gamm_2y_mixed_0 <- gamm4(yi ~ s(soilP)                             ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_0$gam)#R2=0.06
gamm_2y_mixed_0.1 <- gamm4(yi ~ soilP                              ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_0.1$gam)#R2=0.01
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP,tsd, k=20)                       ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1$gam)#R2=0.25
gamm_2y_mixed_1.2 <- gamm4(yi ~ s(soilP)+s(tsd)                     ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1.2$gam)#R2 = 0.34
summary(gamm_2y_mixed_1.2$mer)
gamm_2y_mixed_1.3 <- gamm4(yi ~ s(tsd)                              ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1.3$gam)#R2 = 0.34

#checking best Model 1a - ambient + CTE####
mods_1a <- list(mixed_0 = gamm_2y_mixed_0$mer, mixed_0.1 = gamm_2y_mixed_0.1$mer, 
                mixed_1 = gamm_2y_mixed_1$mer, mixed_1.2 = gamm_2y_mixed_1.2$mer,
                mixed_1.3 = gamm_2y_mixed_1.3$mer)
atab_1a <- aictab(mods_1a) #Mixed 1 is the best model
atab_1a

#Model 1a - ambient####
summary(tot_lit_amb_1to21_amb$Country)
gamm_2y_mixed_0_amb <- gamm4(yi ~ s(soilP)                             ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_0_amb$gam)#R2=0.02

gamm_2y_mixed_0.1_amb <- gamm4(yi ~ soilP                             ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_0.1_amb$gam)#R2=0.01

gamm_2y_mixed_1_amb <- gamm4(yi ~ s(soilP,by=tsd)                             ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1_amb$gam)#R2=0.16

gamm_2y_mixed_1.2_amb <- gamm4(yi ~ s(soilP)+s(tsd)                              ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1.2_amb$gam)#R2 = 0.29
summary(gamm_2y_mixed_1.2_amb$mer)

gamm_2y_mixed_1.3_amb <- gamm4(yi ~ s(soilP,tsd)                              ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1.3_amb$gam)#R2 = 0.43

gamm_2y_mixed_1.4_amb <- gamm4(yi_new ~ s(tsd)                              ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1.4_amb$gam)#R2 = 0.29

#checking best Model 1a - ambient + CTE####
mods_1a_amb <- list(mixed_0_amb = gamm_2y_mixed_0_amb$mer, mixed_0.1_amb = gamm_2y_mixed_0.1_amb$mer, 
                mixed_1_amb = gamm_2y_mixed_1$mer, mixed_1.2_amb = gamm_2y_mixed_1.2$mer,
                mixed_1.3_amb = gamm_2y_mixed_1.3$mer, mixed_1.4_amb = gamm_2y_mixed_1.4_amb$mer)
atab_1a_amb <- aictab(mods_1a_amb)
atab_1a_amb

#Table 3 Model 2a (ambient + CTE)####
gamm_2y_mixed_1e <- gamm4(yi ~ s(soilP, tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1e$gam)#R2 = 0.4
summary(gamm_2y_mixed_1e$mer)

#Other models
gamm_2y_mixed_1e_new.1 <- gamm4(yi ~ s(soilP)+s(tsd)+s(windur)                 ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1e_new.1$gam)#R2 = 0.59
summary(gamm_2y_mixed_1e_new.1$mer)

gamm_2y_mixed_1f_new <- gamm4(yi ~ s(tsd,by=windur, k=20)                 ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1f_new$gam)#R2 = 0.33
summary(gamm_2y_mixed_1f_new$mer)

gamm_2y_mixed_1g_new <- gamm4(yi ~ s(windur)                             ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1g_new$gam)#R2 = 0.15
summary(gamm_2y_mixed_1g_new$mer)

gamm_2y_mixed_1h_new <- gamm4(yi ~ s(soilP)+s(windur)                    ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F)
summary(gamm_2y_mixed_1h_new$gam)#R2 = 0.13
summary(gamm_2y_mixed_1h_new$mer)


#checking best Model 2a - ambient + CTE####
mods_2a <- list(mixed_1e = gamm_2y_mixed_1e$mer, mixed_1enew = gamm_2y_mixed_1e_new.1$mer,
                mixed_1f_new = gamm_2y_mixed_1f_new$mer, mixed_1g_new = gamm_2y_mixed_1g_new$mer,
                mixed_1h_new = gamm_2y_mixed_1h_new$mer)
atab_2a <- aictab(mods_2a)
atab_2a

#Model 2a - ambient####
summary(tot_lit_amb_1to21_amb$windur)
gamm_2y_mixed_1e_amb <- gamm4(yi ~ s(soilP, tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1e_amb$gam)#R2 = 0.33
summary(gamm_2y_mixed_1e_amb$mer)

#Best model for ambient 
gamm_2y_mixed_1e_new_amb <- gamm4(yi_new ~ s(soilP, tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1e_new_amb$gam)#R2 = 0.39
summary(gamm_2y_mixed_1e_new_amb$mer)

gamm_2y_mixed_1e_new_amb2 <- gamm4(yi_new ~ s(soilP) +s(tsd)+s(windur)                 ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1e_new_amb2$gam)#R2 = 0.39
summary(gamm_2y_mixed_1e_new_amb2$mer)

gamm_2y_mixed_1f_new_amb <- gamm4(yi_new ~ s(tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1f_new_amb$gam)#R2 = 0.32
summary(gamm_2y_mixed_1f_new_amb$mer)

gamm_2y_mixed_1g_new_amb <- gamm4(yi_new ~ s(windur)                             ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1g_new_amb$gam)#R2 = 0.32
summary(gamm_2y_mixed_1g_new_amb$mer)

gamm_2y_mixed_1h_new_amb <- gamm4(yi_new ~ s(soilP)+s(windur)                    ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1h_new_amb$gam)#R2 = 0.32
summary(gamm_2y_mixed_1h_new_amb$mer)

gamm_2y_mixed_1i_new_amb <- gamm4(yi_new ~ s(tsd)                    ,weights=tot_lit_amb_1to21_amb$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_amb, REML=F)
summary(gamm_2y_mixed_1i_new_amb$gam)#R2 = 0.32
summary(gamm_2y_mixed_1i_new_amb$mer)

#checking best Model 2a - ambient####
mods_2a <- list(mixed_1e_new_amb = gamm_2y_mixed_1e_new_amb$mer, mixed_1e_new_amb2 =gamm_2y_mixed_1e_new_amb2$mer,
                mixed_1f_new_amb = gamm_2y_mixed_1f_new_amb$mer, mixed_1g_new_amb = gamm_2y_mixed_1g_new_amb$mer,
                mixed_1h_new_amb = gamm_2y_mixed_1h_new$mer,mixed_1i_new_amb = gamm_2y_mixed_1i_new$mer)
atab_2a <- aictab(mods_2a)
atab_2a

#Checking correlation among predictors

names(tot_lit_amb_1to21_final)
#data
cor_1to21<-tot_lit_amb_1to21_final[,c(83,84,85,86)]
str(cor_1to21)## Final
#updating column names
names(cor_1to21)[names(cor_1to21) == "soilP"] <- "Soil P"
names(cor_1to21)[names(cor_1to21) == "hurrwind"] <- "HURRECON Wind speed"
names(cor_1to21)[names(cor_1to21) == "windsp"] <- "WMO wind speed"
names(cor_1to21)[names(cor_1to21) == "tsd"] <- "Time since cyclone"

#calculating pearson's r
corr_1to21 <- round(cor(cor_1to21,method="pearson"),1)
corr_1to21
p.mat_1to21 <- cor_pmat(cor_1to21)

#Plotting correlation plot
Fig_1to21<-ggcorrplot(corr_1to21, hc.order = TRUE, type = "lower",hc.method = "ward.D2",
                  outline.col = "white", p.mat = p.mat_1to21,method="square",ggtheme=ggplot2::theme_bw(),show.legend=TRUE, legend.title="Pearson's r", 
                  lab=TRUE, lab_size=8, tl.cex=20,colors = c("#003f5c", "white", "#ffa600",pch.cex=22,nbreaks = 8,legend.text.cex=22))+font("legend.text",size=16)+font("legend.title", size=16)
Fig_1to21
#saving the figure
ggsave(filename = "Fig_corr_1to21.png",
       plot = Fig_1to21, width = 12, height = 12, units = 'cm',
       scale = 2, dpi = 600)

##Predictions####
#The best model Table 3 Model 2a is gamm_2y_mixed_1e

#Predictions from gamm_2y_mixed_1e####
mypreds_1to21_final_1e<-predict(gamm_2y_mixed_1e$gam,newdata=tot_lit_amb_1to21_final,se.fit=T)
mypreds_1to21_final_1e
tot_lit_amb_1to21_final$Pred_1e<-mypreds_1to21_final_1e$fit
tot_lit_amb_1to21_final$Se_1e<-mypreds_1to21_final_1e$se.fit

#Final Plot Resilience of Total Litterfall by Soil P and Cyclone Intensity####

##Data wrangling - Time since disturbance as factor for plotting
tot_lit_amb_1to21_final$Months<-factor(tot_lit_amb_1to21_final$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                      "18","19","20","21"))
levels(tot_lit_amb_1to21_final$Months)

#Figure 9 Total Litterfall Mass Resilience by Time, soil P and wind duration####
Fig_res1to21 <- ggplot(tot_lit_amb_1to21_final, aes(y=Pred_1e, x=Months,group=Case_study))
Fig_res1to21<-Fig_res1to21+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.7,size=1.2)
Fig_res1to21<-Fig_res1to21+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(tot_lit_amb_1to21_final$vi)),alpha=0.7,shape=21,stroke=1.3)#size=1.5)#color="darkgray")
Fig_res1to21<-Fig_res1to21+geom_ribbon(aes(ymin=Pred_1e-(1.96*Se_1e),ymax=Pred_1e+(1.96*Se_1e),col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.02,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))+scale_fill_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))#+scale_color_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))+scale_fill_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))#+ scale_color_manual(values=c("#FFD700","#A07717","#F06E05","#268FE1","#265E89","#192DB5","#85EE85","#409E40","#0C620C","#D37EA5","#D14082","#EF0672"))+ scale_fill_manual(values=c("#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#EF0672"))
#Fig_res1to21<-Fig_res1to21+geom_line(aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
Fig_res1to21<-Fig_res1to21+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21))+scale_shape_discrete(solid=F)+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=22),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=20),legend.box="horizontal",legend.position="top",legend.justification="center")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  labs(color="Soil P \n(ln mg/kg)",size="Effect size variance",y="Resilience",x="Time since cyclone (Months)")+guides(fill=FALSE,size=FALSE)+ 
  annotate("text", x = 1, y = 2.2, label = "Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
Fig_res1to21

##Figure9 - saving in high res####
ggsave(filename = "Fig9_Resilience_Pred_Tot.png",
       plot = Fig_res1to21, width = 12, height = 10, units = 'cm',
       scale = 2, dpi = 1000)

##NEW FIGURE 9 WITH ONLY 3 LINES AND RIBBONS##
new_fig9<- cbind(data.frame(Case_study=tot_lit_amb_1to21_final$Case_study, Other_soil_P=tot_lit_amb_1to21_final$Other_soil_P),Pred_1e=tot_lit_amb_1to21_final$Pred_1e,Se_1e=tot_lit_amb_1to21_final$Se_1e,Months=tot_lit_amb_1to21_final$Months)
summary(new_fig9$Case_study)
summary(new_fig9$Case_study)

new_fig9_filtered<-new_fig9 %>% filter(Case_study == "Wooroonooran Basalt | Larry | Ambient"|Case_study=="Mt Spec | Charlie | Ambient"|Case_study=="Chamela-Cuixmala | Jova | Ambient")
summary(new_fig9_filtered$Case_study)
            
#new figure 9####
Fig_9 <- ggplot(tot_lit_amb_1to21_final, aes(y=Pred_1e, x=Months,group=Case_study))
Fig_9<-Fig_9+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.7,size=1.2)
Fig_9
Fig_9<-Fig_9+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(tot_lit_amb_1to21_final$vi)),alpha=0.7,shape=21,stroke=1.2)#size=1.5)#color="darkgray")
Fig_9
Fig_9<-Fig_9+geom_ribbon(data=new_fig9_filtered,aes(ymin=Pred_1e-(1.96*Se_1e),ymax=Pred_1e+(1.96*Se_1e),x=Months,col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.02,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))+scale_fill_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))
Fig_9
Fig_9<-Fig_9+geom_line(data=new_fig9_filtered,aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
Fig_9
Fig_9<-Fig_9+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9, 12, 15, 18, 21))+scale_shape_discrete(solid=F)+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=22),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=20),legend.box="horizontal",legend.position="top",legend.justification="center")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  labs(color="Soil P \n(ln mg/kg)",size="Effect size variance",y="Resilience",x="Time since cyclone (Months)")+guides(fill=FALSE,size=FALSE)+ 
  annotate("text", x = 1, y = 2.2, label = "Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
Fig_9


##New Figure9####
ggsave(filename = "Fig9_new.png",
       plot = Fig_9, width = 10, height = 10, units = 'cm',
       scale = 2, dpi = 1000)

#END###