##20210316##

##Post-cyclone temporal dynamics of total litterfall resilience##

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
library(colorRamps)
install.packages("remotes")
remotes::install_github("easystats/performance")
library(remotes)

##Upload data####

metadat<-read.csv(file.choose())#20210520_Litterfall_Mass
attach(metadat)
str(metadat)#2370 obs of 77 variables
names(metadat)
metadat$Case_study= paste(metadat$Site, metadat$DisturbanceName,metadat$Treatment, sep=" | ")

##Data wrangling####

#subset from 1 to 36 months and 1 to 24 months post-cyclone####
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
rec

#1 to 36 months
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)
res_amb<-res_all %>% filter(Treatment=="Ambient")
str(res_amb)#948

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36####
data_esall <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_all, measure = "ROM")
str(data_esall)#1530 obs. including all litterfall mass fractions

##Total Litterfall mass resilience including all Treatments####
#1 to 36 months
tot_lit_all<-data_esall %>% filter(Fraction=="TotLitfall")
str(tot_lit_all)#426 obs
levels(tot_lit_all$Treatment)

#1 to 36 months - Total litterfall Ambient only
tot_lit_amb<-data_esall %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")
str(tot_lit_amb)#250 obs
tot_lit_amb$TSD_months=as.integer(tot_lit_amb$TSD_months)

##Plotting the total litterfall data 1 to 36 months - All Treatments
dres <- tot_lit_all %>%group_by(Site,DisturbanceName,Treatment,TSD_months)  %>%
  dplyr::summarise(counts = dplyr::n())
dres 
p<- dres %>% ggplot(aes(x = TSD_months, fill=Treatment)) +geom_histogram(binwidth = 0.5)+theme_bw()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text = element_text(angle=0, hjust=0.5,vjust = 1,size=18),axis.title=element_text(size=22),
        legend.position="right",legend.text =  element_text(size=20,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  scale_fill_brewer(palette="Paired")+
  labs(x="Time since disturbance (months)",y="Number of observations")#+ annotate("text", x = 10, y = 25, label = "3-year resilience of total litterfall",size=6,colour="black")
FigSupTotTreat<-p + facet_wrap(~Site)+theme(strip.text.x = element_text(size=18))
FigSupTotTreat
ggsave(filename = "SupFig_Res_Tot_Obs.png",
       plot = FigSupTotTreat, width = 22, height = 19, units = 'cm',
       scale = 2, dpi = 600)

##Plotting the total litterfall data 1 to 36 months - Ambient only
dres_amb <- tot_lit_amb %>%group_by(Country, Site,DisturbanceName,Treatment,TSD_months)  %>%
  dplyr::summarise(counts = dplyr::n())
dres_amb
p_amb<- dres_amb %>% ggplot(aes(x = TSD_months, fill=Treatment)) +geom_histogram(binwidth = 0.5)+theme_bw()+guides(color = guide_legend(title = "Country"),legend.key.width=32)+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text= element_text(angle=0, hjust=0.5,vjust = 1,size=18),axis.title=element_text(size=22),
        legend.position="right",legend.text =  element_text(size=20,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+
  #scale_fill_brewer(palette="Paired")+
  labs(x="Time since disturbance (months)",y="Number of observations")#+ annotate("text", x = 10, y = 25, label = "3-year resilience of total litterfall",size=6,colour="black")
FigSupTotAmb<-p_amb + facet_wrap(~Site)+theme(strip.text.x = element_text(size=18))
FigSupTotAmb
ggsave(filename = "SupFig_TotAmbRes_Obs.png",
       plot = FigSupTotAmb, width = 22, height = 19, units = 'cm',
       scale = 2, dpi = 600)

##Analyzing Ambient Conditions only, 1 to 21 months, wherein most data is concentrated

#Data 1 to 21 months - Ambient only
tot_lit_amb_1to21_final<-tot_lit_amb %>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_amb_1to21_final)#218
tot_lit_amb_1to21_final$Case_study= paste(tot_lit_amb_1to21_final$Site, tot_lit_amb_1to21_final$DisturbanceName, sep="|")
unique(levels(as.factor(tot_lit_amb_1to21_final$Case_study)))

#Transforming the numerical predictors
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
tot_lit_amb_1to21_final$soilP<-z.trans(tot_lit_amb_1to21_final$Other_soil_P)
tot_lit_amb_1to21_final$windsp<-z.trans(tot_lit_amb_1to21_final$WMO_wind_kts)
tot_lit_amb_1to21_final$hurrwind<-z.trans(tot_lit_amb_1to21_final$HURRECON_wind_ms)
tot_lit_amb_1to21_final$windur<-z.trans(tot_lit_amb_1to21_final$Gale_wind_duration_minutes)
tot_lit_amb_1to21_final$d2track<-z.trans(tot_lit_amb_1to21_final$Distance_to_Disturb_km)
tot_lit_amb_1to21_final$tsd<-z.trans(tot_lit_amb_1to21_final$TSD_months)

summary(tot_lit_amb_1to21_final)
##########Plot Data ###################################
xyplot(yi~TSD_months, tot_lit_amb_1to21, groups=factor(Other_soil_P), auto.key = TRUE)
xyplot(yi~TSD_months, tot_lit_amb_1to21, groups=Par_Mat, auto.key = TRUE)

##Weights####

#checking case study levels
unique(levels(as.factor(tot_lit_amb_1to21_final$Case_study)))

tot_meta<- rma.mv(yi,vi,random = list(~1|Site,~1|DisturbanceName),
                          tdist = TRUE,
                          data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta)
tot_meta$sigma2
weight<-weights(tot_meta,type="matrix")
tot_lit_amb_1to21_final$weight<-colSums(weight)/sum(weight)
#Best method
tot_lit_amb_1to21_final$weight2<-(1/(tot_lit_amb_1to21_final$vi+0.02666+0.1402))
tot_lit_amb_1to21_final$weight2

tot_meta_b<- rma.mv(yi,vi,random = ~1|Region/DisturbanceName,
                  tdist = TRUE,
                  data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_b)
tot_meta_b$sigma2
weight_b<-weights(tot_meta_b,type="matrix")
tot_lit_amb_1to21_final$weight_b<-colSums(weight_b)/sum(weight_b)
tot_lit_amb_1to21_final$weight_b
tot_lit_amb_1to21_final$weight_b2<-(1/(tot_lit_amb_1to21_final$vi+0.00000006+0.084))

tot_meta_c<- rma.mv(yi,vi,random = ~1|DisturbanceName,
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_c)
weight_c<-weights(tot_meta_c,type="matrix")
tot_lit_amb_1to21_final$weight_c<-colSums(weight_c)/sum(weight_c)
tot_lit_amb_1to21_final$weight_c

tot_meta_d<- rma.mv(yi,vi,random = ~(1|Region),
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_d)
weight_d<-weights(tot_meta_d,type="matrix")
tot_lit_amb_1to21_final$weight_d<-colSums(weight_d)/sum(weight_d)
tot_lit_amb_1to21_final$weight_d

tot_meta_e<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = tot_lit_amb_1to21_final,struct = "HAR",method = "REML")
summary(tot_meta_e)
weight_e<-weights(tot_meta_e,type="matrix")
tot_lit_amb_1to21_final$weight_e<-colSums(weight_e)/sum(weight_e)
tot_lit_amb_1to21_final$weight_e

# Sample GAMMs to analyse trends across months
###Sample GAMMs from Walker et al. 2019####
# specify the number of knots in the smooth terms of the GAMM, this sets how wiggly the gamms can be, 4 is used here to avoid tracking strong inter-annual variation
k <- 4 

# full model with fixed co2 effect, smoothed fixed effect of year, and a fixed interaction effect between co2 and year, random effect on intercept of sub-sample nested within plot 
gamm_full <- gamm4(NPP ~ co2 + s(YEAR, k=k) + s(YEAR, by=co2, k=k), random = ~(1|plot/sub), data=df1, REML=F )
# additive model with fixed co2 effect and smoothed fixed effect of year, random effect on intercept of sub-sample nested within plot 
gamm_add  <- gamm4(NPP ~ co2 + s(YEAR, k=k)                       , random = ~(1|plot/sub), data=df1, REML=F )
# additive model with fixed co2 effect and linear fixed effect of year, random effect on intercept of sub-sample nested within plot 
gamm_lin  <- gamm4(NPP ~ co2 + YEAR                               , random = ~(1|plot/sub), data=df1, REML=F )


#Fitting GAMMs using weights generated from rma.mv function####
tot_lit_amb_1to21_final$weight
w<- (1/tot_lit_amb_1to21_final$vi)
w
#Site as random effect with site-based weights
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP, tsd)                              ,weights=(1/tot_lit_amb_1to21_final$vi), random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1.2 <- gamm4(yi ~ s(soilP, tsd)                              ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
mods_1 <- list(mixed_1 = gamm_2y_mixed_1$mer, mixed_1.2 = gamm_2y_mixed_1.2$mer)
atab_1 <- aictab(mods_1)
atab_1#MUCH BETTER TO USE THE SECOND APPROACH TO WEIGHTS

gamm_2y_mixed_2 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=(1/tot_lit_amb_1to21_final$vi), random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2.2 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_b2, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2.1 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_b, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
mods_2 <- list(mixed_2 = gamm_2y_mixed_2$mer, mixed_2.1 = gamm_2y_mixed_2.1$mer,mixed_2.2 = gamm_2y_mixed_2.2$mer)
atab_2 <- aictab(mods_2)
atab_2

gamm_2y_mixed_3 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_5 <- gamm4(yi ~ s(soilP, by=tsd)                           ,weights=tot_lit_amb_1to21_final$weight_e, random = ~(1|Site), data=tot_lit_amb_1to21_final, REML=F )

gamm_2y_mixed_1a <- gamm4(yi ~ s(soilP) +s(tsd)                           ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2a <- gamm4(yi ~ s(Other_soil_P,by=TSD_months)              ,weights=tot_lit_amb_1to21_final$weight_b, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_3a <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4a <- gamm4(yi ~ s(Other_soil_P, TSD_months,k=20)            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_5a <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_e, random = ~(1|Site), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4b <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )

gamm_2y_mixed_1.1a <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=Fujita_scale),weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_2.1a <- gamm4(yi ~ s(soilP, tsd,by=Fujita_scale)                               ,weights=tot_lit_amb_1to21_final$weight_b, random = ~(1|Region/DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_3.1a <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=Fujita_scale)                 ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_4.1a <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=Fujita_scale)            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_5.1a <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=Fujita_scale)                 ,weights=tot_lit_amb_1to21_final$weight_e, random = ~(1|Site), data=tot_lit_amb_1to21_final, REML=F )

gamm_2y_mixed_1c <- gamm4(yi ~ soilP + s(tsd)                                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1d <- gamm4(yi ~ soilP + tsd                                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )
gamm_2y_mixed_1e <- gamm4(yi ~ s(soilP, tsd,by=windur)                 ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=F )

summary(gamm_2y_mixed_1$gam)
summary(gamm_2y_mixed_1.2$gam)
summary(gamm_2y_mixed_2a$gam)
summary(gamm_2y_mixed_3$gam)
summary(gamm_2y_mixed_4$gam)
summary(gamm_2y_mixed_4b$gam)
summary(gamm_2y_mixed_5$gam)

#Comparing models
mods_1to21 <- list(mixed_1 = gamm_2y_mixed_1$mer, mixed_1a = gamm_2y_mixed_1a$mer,mixed_1b = gamm_2y_mixed_1b$mer,
                   mixed_1c = gamm_2y_mixed_1c$mer,mixed_1d = gamm_2y_mixed_1d$mer,mixed_1e = gamm_2y_mixed_1e$mer)
                   
                   #mixed_1a = gamm_2y_mixed_1a$mer,mixed_2 = gamm_2y_mixed_2$mer,mixed_2a = gamm_2y_mixed_2a$mer, 
                   #mixed_3 = gamm_2y_mixed_3$mer,mixed_3a = gamm_2y_mixed_3a$mer,mixed_4 = gamm_2y_mixed_4$mer,
                   #mixed_4a = gamm_2y_mixed_4a$mer,mixed_4b = gamm_2y_mixed_4b$mer,mixed_5 = gamm_2y_mixed_5$mer,mixed_5a = gamm_2y_mixed_5a$mer)
atab_1to21 <- aictab(mods_1to21)
atab_1to21

##Both are the best models
summary(gamm_2y_mixed_1a$gam)
summary(gamm_2y_mixed_2a$gam)
summary(gamm_2y_mixed_3a$gam)
summary(gamm_2y_mixed_4a$gam)
summary(gamm_2y_mixed_5a$gam)

unique(levels(as.factor(tot_lit_amb_1to21_final$Case_study)))

#New Table 5####
#Model 1a
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP,by=tsd)                             ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_1$gam)
gamm_2y_mixed_1.2 <- gamm4(yi ~ s(soilP, tsd,k=20)                              ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_1.2$gam)
summary(gamm_2y_mixed_1.2$mer)

#Model 2a
gamm_2y_mixed_1e <- gamm4(yi ~ s(soilP, tsd,by=windur,k=20)                 ,weights=tot_lit_amb_1to21_final$weight2, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_1e$gam)#R2 = 0.4
summary(gamm_2y_mixed_1e$mer)
vis.gam(gamm_2y_mixed_1e$gam,plot.type="contour",view=c("tsd","soilP"))


#
gamm_2y_mixed_1 <- gamm4(yi ~ s(soilP, tsd)                           ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_1$gam)
summary(gamm_2y_mixed_1$mer)

#
gamm_2y_mixed_1_A <- gamm4(yi ~ s(soilP, tsd,windur)                           ,weights=tot_lit_amb_1to21_final$weight, random = ~(1|Site)+(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_1_A$gam)
summary(gamm_2y_mixed_1_A$mer)

#
gamm_2y_mixed_3a <- gamm4(yi ~ s(Other_soil_P, TSD_months)                 ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|DisturbanceName), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_3a$gam)

#
gamm_2y_mixed_3b <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=factor(Fujita_scale))                 ,weights=tot_lit_amb_1to21_final$weight_c, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_3b$gam)

anova()

#mixed_1.1a = gamm_2y_mixed_1.1a$mer,mixed_2.1a = gamm_2y_mixed_2.1a$mer,mixed_3.1a = gamm_2y_mixed_3.1a$mer,
#mixed_4.1a = gamm_2y_mixed_4.1a$mer,mixed_5.1a = gamm_2y_mixed_5.1a$mer,mixed_1.2a = gamm_2y_mixed_1.2a$mer)

gamm_2y_mixed_4.1a <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=factor(Fujita_scale))            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE)
summary(gamm_2y_mixed_4.1a$gam)
gamm_2y_mixed_4a <- gamm4(yi ~ s(soilP, tsd,by=wind_cat)            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4a$gam)#R2 = 0.15
summary(gamm_2y_mixed_4a$mer)

gamm_2y_mixed_4a.1 <- gamm4(yi ~ s(soilP, tsd,by=factor(Fujita_scale))            ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4a.1$gam)#R2 = 0.22
summary(gamm_2y_mixed_4a.1$mer)

#Adding wind speed variable does not improve the model
gamm_2y_mixed_4b <- gamm4(yi ~ s(Other_soil_P, TSD_months,k=20)+HURRECON_wind_ms           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4b$gam)
summary(gamm_2y_mixed_4b$mer)

gamm_2y_mixed_4b.1 <- gamm4(yi ~ s(Other_soil_P, TSD_months,HURRECON_wind_ms)           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4b.1$gam)
summary(gamm_2y_mixed_4b$mer)

#Adding Fujita variable does not improve the model
gamm_2y_mixed_4c <- gamm4(yi ~ s(Other_soil_P, TSD_months)+Fujita_scale           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4c$gam)
summary(gamm_2y_mixed_4c$mer)

gamm_2y_mixed_4c.1 <- gamm4(yi ~ s(Other_soil_P, TSD_months,by=Fujita_scale)           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4c.1$gam)
summary(gamm_2y_mixed_4c$mer)

#Adding Wind duration variable does not improve the model
gamm_2y_mixed_4c <- gamm4(yi ~ s(soilP, tsd,k=20)+windur           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4c$gam)
summary(gamm_2y_mixed_4c$mer)

gamm_2y_mixed_4d <- gamm4(yi ~ s(soilP, tsd)+s(windur)           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4d$gam)
summary(gamm_2y_mixed_4d$mer)

gamm_2y_mixed_4e <- gamm4(yi ~ s(soilP, tsd,windur)           ,weights=tot_lit_amb_1to21_final$weight_d, random = ~(1|Region), data=tot_lit_amb_1to21_final, REML=TRUE )
summary(gamm_2y_mixed_4e$gam)#R2=0.24
summary(gamm_2y_mixed_4e$mer)


#There is significant interaction between soil P and Time since disturbance

#Now testing interaction with intensity

vis.gam(gamm_2y_mixed_4e$gam,plot.type="contour",view=c("soilP","windur"))
vis.gam(gamm_2y_mixed_4a.1$gam,plot.type="contour",view=c("soilP","tsd"))
#vis.gam(gamm_2y_mixed_4e, view=c("soilP","tsd","windur"),labcex = 1.1, cex.lab = 1.4, cex.axis =1.4, cex.main = 1.8, main = "Best GAM model: ln ANF ~ s(CN)*** + s(TP)*** +s(moisture)***", xlab = "soil CN (molar ratio)", ylab = "Soil Total P (log mol kg-1)", ticktype="detailed", color="terrain", theta=-35, phi=30, type="response", plot.type="persp", nCol=50, contour.col="black")
#plot.gam(gamm_2y_mixed_4a.1, select=1, cex.axis=1.4,cex.lab=1.4, xlab="Total N (log mol kg-1)",ylab="s(log TN,1.0):Eucalypt", se=TRUE, all.terms=T,col='#FF8000', shade=TRUE,shade.col='gray90')

#Checking correlation among predictors

dat_1to21<- tot_lit_amb_1to21_final %>%group_by(Case_study)  %>%
  dplyr::summarise(counts = dplyr::n())
dat_1to21
names(tot_lit_amb_1to21_final)
tot_lit_amb_1to21_final$Case

names(tot_lit_amb_1to21_final)
tot_lit_amb_1to21_final$fujita=as.numeric(tot_lit_amb_1to21_final$fujita)
cor_1to21<-tot_lit_amb_1to21_final[,c(40,81,83,84,86)]
str(cor_1to21)## Final!
names(cor_1to21)[names(cor_1to21) == "soilP"] <- "Soil P"
names(cor_1to21)[names(cor_1to21) == "hurrwind"] <- "HURRECON Wind speed"
names(cor_1to21)[names(cor_1to21) == "windsp"] <- "WMO wind speed"
names(cor_1to21)[names(cor_1to21) == "tsd"] <- "Time since cyclone"
names(cor_1to21)[names(cor_1to21) == "fujita"] <- "Fujita scale"

corr_1to21 <- round(cor(cor_1to21,method="pearson"),1)
corr_1to21
p.mat_1to21 <- cor_pmat(cor_1to21)

Fig_1to21<-ggcorrplot(corr_1to21, hc.order = TRUE, type = "lower",hc.method = "ward.D2",
                  outline.col = "white", p.mat = p.mat_1to21,method="square",ggtheme=ggplot2::theme_bw(),show.legend=TRUE, legend.title="Pearson's r", 
                  lab=TRUE, lab_size=8, tl.cex=20,colors = c("#003f5c", "white", "#ffa600",pch.cex=22,nbreaks = 8,legend.text.cex=22))+font("legend.text",size=16)+font("legend.title", size=16)
Fig_1to21

ggsave(filename = "Fig_corr_1to21.png",
       plot = Fig_1to21, width = 12, height = 12, units = 'cm',
       scale = 2, dpi = 600)

##Predictions####
# GAMM plotting best model gamm_2y_mixed_4e
summary(gamm_2y_mixed_4a.1$gam)

#Predictions from gamm_2y_mixed_1e####
mypreds_1to21_final_1e<-predict(gamm_2y_mixed_1e$gam,newdata=tot_lit_amb_1to21_final,se.fit=T)
mypreds_1to21_final_1e
tot_lit_amb_1to21_final$Pred_1e<-mypreds_1to21_final_1e$fit
tot_lit_amb_1to21_final$Se_1e<-mypreds_1to21_final_1e$se.fit

#Data wrangling for plotting
Est_res1to21 <- cbind(data.frame(Estimate=tot_lit_amb_1to21_final$Pred,Se=tot_lit_amb_1to21_final$Se,Months=factor(tot_lit_amb_1to21_final$TSD_months),
                                  SoilP=log(tot_lit_amb_1to21_final$Other_soil_P),fac_SoilP=factor(tot_lit_amb_1to21_final$Other_soil_P),Wind_speed=tot_lit_amb_1to21_final$HURRECON_wind_ms,
                                  Region=tot_lit_amb_1to21_final$Country,Case_study=tot_lit_amb_1to21_final$Case_study,Cyclone=tot_lit_amb_1to21_final$DisturbanceName,Fujita=factor(tot_lit_amb_1to21_final$fujita)))
str(Est_res1to21)
Est_res1to21$Months<-factor(Est_res1to21$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                  "18","19","20","21"))
levels(Est_res1to21$Months)

#Final Plot Resilience of Total Litterfall by Soil P and Cyclone Intensity####

##Data wrangling - Time since disturbance as factor for plotting
tot_lit_amb_1to21_final$Months<-factor(tot_lit_amb_1to21_final$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                      "18","19","20","21"))
levels(tot_lit_amb_1to21_final$Months)

#Figure 9 Total Litterfall Mass Resilience by Time, soil P and wind duration####
Fig_res1to21 <- ggplot(tot_lit_amb_1to21_final, aes(y=Pred_1e, x=Months,group=Case_study))
Fig_res1to21<-Fig_res1to21+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.9,size=1.2)
Fig_res1to21<-Fig_res1to21+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(tot_lit_amb_1to21_final$vi)),alpha=0.9,shape=21,stroke=1.3)#size=1.5)#color="darkgray")
Fig_res1to21<-Fig_res1to21+geom_ribbon(aes(ymin=Pred_1e-(1.96*Se_1e),ymax=Pred_1e+(1.96*Se_1e),col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.02,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))+scale_fill_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))#+scale_color_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))+scale_fill_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))#+ scale_color_manual(values=c("#FFD700","#A07717","#F06E05","#268FE1","#265E89","#192DB5","#85EE85","#409E40","#0C620C","#D37EA5","#D14082","#EF0672"))+ scale_fill_manual(values=c("#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#EF0672"))
Fig_res1to21
Fig_res1to21<-Fig_res1to21+geom_line(aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
Fig_res1to21<-Fig_res1to21+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21))+scale_shape_discrete(solid=F)+
  ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("")+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text=element_text(size=26),legend.text =  element_text(size=20),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=20),legend.box="horizontal",legend.position="top",legend.justification="center")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  labs(color="Soil P \n(ln mg/kg)",size="Effect size variance")+guides(fill=FALSE,size=FALSE)+ 
  annotate("text", x = 1, y = 2.2, label = "a Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
Fig_res1to21

Fig9<-Fig_res1to21+Fig_res1to21_l+plot_layout(ncol=1,heights=c(1,1))
Fig9

#Fig_res1to21_facetb<-Fig_res1to21+ facet_wrap(~Fujita_scale)+theme(strip.text.x = element_text(size=20))
#Fig_res1to21_facetb
##Figure9a####
ggsave(filename = "Fig9_Resilience_Pred_Tot_Leaf.png",
       plot = Fig9, width = 12, height = 14, units = 'cm',
       scale = 2, dpi = 1200)

##Soil P and TSD only

#Predictions from gamm_2y_mixed_1e####
mypreds_1to21_final_1.2<-predict(gamm_2y_mixed_1.2$gam,newdata=tot_lit_amb_1to21_final,se.fit=T)
mypreds_1to21_final_1.2
tot_lit_amb_1to21_final$Pred_1.2<-mypreds_1to21_final_1.2$fit
tot_lit_amb_1to21_final$Se_1.2<-mypreds_1to21_final_1.2$se.fit

##Data wrangling - Time since disturbance as factor for plotting
tot_lit_amb_1to21_final$Months<-factor(tot_lit_amb_1to21_final$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                                      "18","19","20","21"))
levels(tot_lit_amb_1to21_final$Months)

#Figure 9 Total Litterfall Mass Resilience by Time, soil P and wind duration####
Fig_res1to21a <- ggplot(tot_lit_amb_1to21_final, aes(y=Pred_1.2, x=Months,group=Case_study))
Fig_res1to21a<-Fig_res1to21a+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.9,size=1.2)
Fig_res1to21a<-Fig_res1to21a+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(tot_lit_amb_1to21_final$vi)),alpha=0.9,shape=21,stroke=1.3)#size=1.5)#color="darkgray")
Fig_res1to21a<-Fig_res1to21a+geom_ribbon(aes(ymin=Pred_1.2-(1.96*Se_1.2),ymax=Pred_1.2+(1.96*Se_1.2),col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.02,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))+scale_fill_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,6,7))#+scale_color_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))+scale_fill_manual(values=c("#4575b4","#abd9e9","#00876c","#3c986d","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#a50026"))#+ scale_color_manual(values=c("#FFD700","#A07717","#F06E05","#268FE1","#265E89","#192DB5","#85EE85","#409E40","#0C620C","#D37EA5","#D14082","#EF0672"))+ scale_fill_manual(values=c("#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#EF0672"))
Fig_res1to21a
Fig_res1to21a<-Fig_res1to21a+geom_line(aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
Fig_res1to21a<-Fig_res1to21a+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21))+scale_shape_discrete(solid=F)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("")+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y = element_text(size=22),legend.text =  element_text(size=20),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=20),legend.box="horizontal",legend.position="top",legend.justification="center")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  labs(color="Soil P \n(ln mg/kg)",size="Effect size variance")+guides(fill=FALSE,size=FALSE)+ 
  annotate("text", x = 1, y = 2.2, label = "a Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
Fig_res1to21a

##Figure9####
Fig9a<-Fig_res1to21a+Fig_res1to21_l+plot_layout(ncol=1,heights=c(1,1))
Fig9a
#Saving in high res
ggsave(filename = "Fig9_Resilience_Pred_Tot_Leaf_v2.png",
       plot = Fig9a, width = 12, height = 14, units = 'cm',
       scale = 2, dpi = 1200)


##END##