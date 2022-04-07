###Meta-Analysis of cyclone resistance of forest litterfall across the tropics###

library(ggplot2)
library(ggstatsplot)
library(nlme)
library(gapminder)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(magrittr)
library(dplyr)
library(hrbrthemes)
library(patchwork)
library(lme4)
library(tidyr)
library(plyr)
library(ggridges)
library(lmerTest)
library(metafor)
library(dmetar)
library(leaps)
library(glmulti)
library(metaviz)
library(Matrix)
library(ggExtra)
library(RColorBrewer)
library(AICcmodavg)

####STEP 0 uploading data####

#Litterfall mass flux data
metadat<-read.csv(file.choose())#Litterfall_Mass in data folder
str(metadat)#2367 obs of 62 variables

#transforming variables to numeric
metadat$HURRECON_wind_ms=as.numeric(metadat$HURRECON_wind_ms)
metadat$Gale_wind_duration_minutes=as.numeric(metadat$Gale_wind_duration_minutes)

#Create Case study column
metadat$Case_study= paste(metadat$Site, metadat$DisturbanceName,sep=" | ")

#Nutrient flux and concentration data
nutmeta<-read.csv(file.choose())#Litterfall_Nutrients in data folder
str(nutmeta) #2551 obs. of  78 variables

#create Case study column
nutmeta$Case_study= paste(nutmeta$Site, nutmeta$DisturbanceName, sep="| ")
unique(levels(as.factor(nutmeta$Case_study)))

####STEP 1 Data Wrangling ####

###Litterfall Mass flux####

##Total Litterfall####
#Annual-based, excluding CTE
data0a<-metadat %>% filter(Fraction=="TotLitfall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")#%>%filter(Treatment!="P+")#|Treatment!="N+"|Treatment!="NP+")
str(data0a)#48 observations
unique(levels(as.factor(data0a$Case_study)))
#Sub-annual
data0aS<- data0a%>% filter(Pre_Mean_MonthSpecific!="NA")#to exclude a factor level
str(data0aS)#23observations
data0aS$Effectsize_ID
#Annual with same observations as in sub-annual
data0aSS<-data0a %>% filter(Effectsize_ID == "34"|Effectsize_ID == "44"|Effectsize_ID == "47"|Effectsize_ID == "1735"|Effectsize_ID == "1845"|Effectsize_ID == "2099"|
                                    Effectsize_ID == "2103"|Effectsize_ID == "2107"|Effectsize_ID == "2111"|
                              Effectsize_ID == "2115"|Effectsize_ID == "2119"|Effectsize_ID=="2123"|Effectsize_ID == "2127"|Effectsize_ID == "2131"|Effectsize_ID == "2140"|
                                    Effectsize_ID == "2343"|Effectsize_ID == "2344"|
                                    Effectsize_ID == "2345"|Effectsize_ID == "2346"|Effectsize_ID == "2347"|Effectsize_ID == "2348"|Effectsize_ID == "2349")
str(data0aSS)#22 observations including the same effect sizes with and without sub-annual pre and post data

##Leaf litterfall###
#Annual-based, excluding CTE
data0alf<-metadat %>% filter(Fraction=="Leaf fall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0alf)#30 observations
#Sub-annual
data0alfS<- data0alf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0alfS)#16 observations
data0alfS$Effectsize_ID
#Annual with same observations as in sub-annual
data0alfSS<-data0alf %>% filter(Effectsize_ID == "40"|Effectsize_ID == "45"|Effectsize_ID == "1711"|Effectsize_ID == "1757"|Effectsize_ID == "1867"|Effectsize_ID == "2100"|
                                        Effectsize_ID == "2104"|Effectsize_ID == "2108"|Effectsize_ID == "2112"|Effectsize_ID=="2116"|Effectsize_ID=="2120"|Effectsize_ID=="2124"|
                                        Effectsize_ID == "2128"|Effectsize_ID == "2149"|Effectsize_ID == "2151"|Effectsize_ID == "2167")
str(data0alfSS)#16 observations

##Wood litterfall####
#Annual-based, excluding CTE
data0awf<-metadat %>% filter(Fraction=="Wood fall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0awf)#29 observations
#Sub-annual
data0awfS<- data0awf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0awfS)#14 observations
data0awfS$Effectsize_ID
#Annual with same observations as in sub-annual
data0awfSS<-data0awf %>% filter(Effectsize_ID == "41"|Effectsize_ID == "42"|Effectsize_ID == "1779"|Effectsize_ID == "1889"|
                                        Effectsize_ID == "2101"|Effectsize_ID == "2105"|Effectsize_ID == "2109"|Effectsize_ID == "2113"|
                                        Effectsize_ID == "2117"|Effectsize_ID == "2121"|Effectsize_ID == "2125"|Effectsize_ID == "2129"|Effectsize_ID == "2153"|Effectsize_ID == "2155")
str(data0awfSS)#14 observations

## FFS litterfall ####
#Annual-based, excluding CTE
data0aff<-metadat %>% filter(Fraction=="FFS fall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")%>% filter(Pre_Mean!="NA")
str(data0aff)#14 observations
#Sub-annual
data0affS<-data0aff %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0affS)#10 observations
data0affS$Effectsize_ID
#Annual with same observations as in sub-annual
data0affSS<-data0aff %>% filter(Effectsize_ID == "1801"|Effectsize_ID == "1911"|Effectsize_ID == "2102"|Effectsize_ID == "2106"|
                                        Effectsize_ID == "2110"|Effectsize_ID == "2114"|Effectsize_ID == "2118"|Effectsize_ID == "2122"|
                                        Effectsize_ID == "2126"|Effectsize_ID == "2130")
str(data0affSS)#10 observations

## Miscellaneous litterfall ####
#Annual-based, excluding CTE
data0amf<-metadat %>% filter(Fraction=="Misc fall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0amf)#9 observations
#Sub-annual
data0amfS<-data0amf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0amfS)#4 obs
data0amfS$Effectsize_ID
data0amfSS<-data0amf %>% filter(Effectsize_ID == "43"|Effectsize_ID == "46"|Effectsize_ID == "1823"|Effectsize_ID == "1933")
str(data0amfSS)#4 obs

####Litterfall Nutrients####

####P flux####
str(nutmeta)
unique(levels(as.factor(nutmeta$Treatment)))

#Total litterfall P flux
data0tpf<-nutmeta %>%filter(Fraction=="TotLitfall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0tpf)#10 observations
unique(levels(as.factor(data0tpf$Case_study)))
#Leaf litterfall P flux
data0lpf<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0lpf)#9
unique(levels(as.factor(data0lpf$Case_study)))
#Wood litterfall
data0wpf<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0wpf)#8
#FFS litterfall
data0ffspf<-nutmeta %>%filter(Fraction=="FFS fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0ffspf)#1
#Miscellaneous litterfall
data0mpf<-nutmeta %>%filter(Fraction=="Misc fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0mpf)#3

#### N flux####

#Total litterfall
data0tnf<-nutmeta %>%filter(Fraction=="TotLitfall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0tnf)#9 observations
unique(levels(as.factor(data0tnf$Case_study)))
#Leaf litterfall
data0lnf<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0lnf)#9 observations
#Wood litterfall
data0wnf<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0wnf)#7 observations
#FFS litterfall
data0ffsnf<-nutmeta %>%filter(Fraction=="FFS fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0ffsnf)#1 observation
#Miscellaneous litterfall
data0mnf<-nutmeta %>%filter(Fraction=="Misc fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0mnf)#3 observations

####P concentration####

#Leaf litterfall
data0lpc<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0lpc)#10
unique(levels(as.factor(data0lpc$Case_study)))
#Wood litterfall
data0wpc<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Response_variable=="P") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0wpc)#3

##N concentration ####

#Leaf litterfall
data0lnc<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0lnc)#10
unique(levels(as.factor(data0lnc$Case_study)))
#Wood litterfall
data0wnc<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Response_variable=="N") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
str(data0wnc)#3

####STEP 2 Individual Effect size calculation####

##Mass Flux####

#Total Litterfall
str(data0a)#48 observations

#Annual
data_es0ia <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = data0a, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).
str(data_es0ia)
data0a_all<-metadat %>% filter(Fraction=="TotLitfall")%>%filter(Cat_TSD_months=="0-0.5")
data_es0ia_all <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = data0a_all, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).
data_es0ia_amb<-data_es0ia%>%filter(Treatment!="P+")%>%filter(Treatment!="N+")%>%filter(Treatment!="NP+")
str(data_es0ia_amb)
#Sub-annual data
data_es0iaS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                      sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0aS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)
str(data_es0iaS)
#long-term pre mean for sites that also have subannual data
data_es0iaSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0aSS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)

#Leaf fall
data_es0ilf <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                      sd1i = Post_SD, sd2i = Pre_SD, data = data0alf, measure = "ROM")
data_es0ilfS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                       sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0alfS, measure = "ROM")
data_es0ilfSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                        sd1i = Post_SD, sd2i = Pre_SD, data = data0alfSS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)
#Wood fall
data_es0iwf <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                      sd1i = Post_SD, sd2i = Pre_SD, data = data0awf, measure = "ROM")
data_es0iwfS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                       sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0awfS, measure = "ROM")
data_es0iwfSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                        sd1i = Post_SD, sd2i = Pre_SD, data = data0awfSS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)
#FFS fall
data_es0iff <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                      sd1i = Post_SD, sd2i = Pre_SD, data = data0aff, measure = "ROM")
data_es0iffS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                       sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0affS, measure = "ROM")
data_es0iffSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                        sd1i = Post_SD, sd2i = Pre_SD, data = data0affSS, measure = "ROM")
#Misc. fall
data_es0imf <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                      sd1i = Post_SD, sd2i = Pre_SD, data = data0amf, measure = "ROM")
data_es0imfS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                       sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0amfS, measure = "ROM")
data_es0imfSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                        sd1i = Post_SD, sd2i = Pre_SD, data = data0amfSS, measure = "ROM")

#Total P flux and fractions####
data_es0itpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0tpf, measure = "ROM")
unique(levels(as.factor(data_es0itpf$Case_study)))
data_es0ilpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lpf, measure = "ROM")
data_es0iwpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wpf, measure = "ROM")
data_es0iffspf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = data0ffspf, measure = "ROM")
data_es0impf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0mpf, measure = "ROM")

##P concentration####
#Leaf litterfall
data_es0ilpc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lpc, measure = "ROM")
unique(levels(as.factor(data_es0ilpc$Case_study)))
#Wood litterfall
data_es0iwpc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wpc, measure = "ROM")

##N flux####
#Total and fractions
data_es0itnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0tnf, measure = "ROM")
unique(levels(as.factor(data_es0itnf$Case_study)))
data_es0ilnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lnf, measure = "ROM")
data_es0iwnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wnf, measure = "ROM")
data_es0iffsnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = data0ffsnf, measure = "ROM")
data_es0imnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0mnf, measure = "ROM")

##N concentration####
#Leaf fall
data_es0ilnc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lnc, measure = "ROM")
#Wood fall
data_es0iwnc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wnc, measure = "ROM")

#### STEP 3 Overall Effect Size calculations ####

## Total Litterfall Mass Flux##

##Final model for Pantropical response of total litterfall mass####
full.model3 <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                       tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                       data = data_es0ia,method = "REML")

summary(full.model3)
((exp(full.model3$b)-1)*100)
(exp(full.model3$se)-1)*100

#Testing if study 53 changes the mean pantropical resistance 
data_no53<-data_es0ia %>% filter(Study_ID!=53)
str(data_no53)
full.model_no_53 <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_no53,method = "REML")

#Testing if four case studies in Hawaii make a difference in the pantropical resistance
red_sites2 <- data_es0ia %>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")
full.model_red2 <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                           tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                           data = red_sites2,method = "REML")

summary(full.model_red2)
((exp(full.model_red2$b)-1)*100)
(exp(full.model_red2$se)-1)*100

#testing alternative random effects
full.model3_b<- rma.mv(yi, vi,random = ~1|Region/DisturbanceName,
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ia,method = "REML")

summary(full.model3_b)
((exp(full.model3_b$b)-1)*100)
(exp(full.model3_b$se)-1)*100

#testing with a random sample of the case studies
resistance_PR<- data_es0ia %>% filter(Country == "Puerto Rico")
summary(resistance_PR$Country)
random_sample_PR<-sample_n(resistance_PR, 10)
str(random_sample_PR)# 10 obs. of  68 variables
names2<-names(random_sample_PR)
names(random_sample_PR)<-c(names2)

resistance_rest<-data_es0ia %>% filter(Country!= "Puerto Rico")
str(resistance_rest)# 29 obs of 68 variables
names1<-names(resistance_rest)
names(resistance_rest)<-c(names1)

#Final data frame combining 10 random case studies from PR and the remainder of the extra-PR case studies
new_resistance_data<-rbind(data.frame(resistance_rest),data.frame(random_sample_PR))
summary(new_resistance_data$Country)

#new pantropical resistance with a random sample from PR including 10 case studies
full.model_randomPR <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = new_resistance_data,method = "REML")#N = 39

summary(full.model_randomPR)
((exp(full.model_randomPR$b)-1)*100)
(exp(full.model_randomPR$se)-1)*100

#soil P as moderator
full.model_randomPR_P <- rma.mv(-1*yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                              tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                              data = new_resistance_data, mods =~soilP, method = "REML")#N = 39

summary(full.model_randomPR)
(sum(full.model_randomPR$sigma2) - sum(full.model_randomPR_P$sigma2)) / sum(full.model_randomPR$sigma2)


#Pre-mean, sd, se
summary(data0a)
mean(data0a$Pre_Mean)
sd<-mean(data0a$Pre_SD)
n<-mean(data0a$S_size)
n
se<-sd/(sqrt(n))
se

#contrasting with other random effects options
full.model3a <- rma.mv(yi, vi,random = ~ 1 | DisturbanceName,
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ia,method = "REML")

summary(full.model3a)

#site as a random effect
full.model3b <- rma.mv(yi, vi,random = ~ 1 | Site,
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ia,method = "REML")

summary(full.model3b)

full.model3c <- rma.mv(yi, vi,random = list(~1|Country,~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ia,method = "REML")

summary(full.model3c)
((exp(full.model3c$b)-1)*100)
(exp(full.model3c$se)-1)*100

anova(full.model3,full.model3c)

## diagnostics
mlm.variance.distribution(x = full.model3)

#Checking profile likelihood plots of the variance components of the model

par(mfrow=c(2,1))
profile(full.model3, sigma2=1)
profile(full.model3, sigma2=2)

#New data frame removing the five sites where total soil P was estimated from available P
red_sites <- data_es0ia %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")
red_sites2 <- data_es0ia %>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")

#Running same model with reduced data set
full.model3_red<- rma.mv(yi, 
                         vi, 
                         random = list(~ 1 | Site,~1|DisturbanceName),
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = red_sites,
                         method = "REML")
summary(full.model3_red)

##Calculating Overall Effect Sizes for sub-annual data
full.model3S <- rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName),
                       tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                       data = data_es0iaS,method = "REML")
summary(full.model3S)
(exp(full.model3S$b)-1)*100
(exp(full.model3S$se)-1)*100

full.model3SS <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iaSS,method = "REML")
summary(full.model3SS)
(exp(full.model3SS$b)-1)*100
(exp(full.model3SS$se)-1)*100

## Leaf fall ####
#Annual data
full.model3lf <- rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0ilf,method = "REML")
summary(full.model3lf)
(exp(full.model3lf$b)-1)*100
(exp(full.model3lf$se)-1)*100

full.model3lf_b <- rma.mv(yi,vi,random = ~1|Region/DisturbanceName, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0ilf,method = "REML")
summary(full.model3lf_b)

#Sub-annual data
full.model3lfS <- rma.mv(yi,vi, 
                         random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilfS,method = "REML")
summary(full.model3lfS)

full.model3lfSS <- rma.mv(yi, 
                          vi, 
                          random =  list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0ilfSS,method = "REML")
summary(full.model3lfSS)

###Wood fall####

#Annual data
full.model3wf <- rma.mv(yi,vi,random =  list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iwf,method = "REML")
summary(full.model3wf)
(exp(full.model3wf$b)-1)*100
(exp(full.model3wf$se)-1)*100

#Sub-annual data
full.model3wfS <- rma.mv(yi,vi,random =list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwfS,method = "REML")
summary(full.model3wfS)

full.model3wfSS <- rma.mv(yi,vi,random =list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0iwfSS,method = "REML")
summary(full.model3wfSS)

### Miscellaneous fall####
data_es0imf$Case_study
#Annual
full.model3mf <- rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0imf,method = "REML")
summary(full.model3mf)
(exp(full.model3mf$b)-1)*100
(exp(full.model3mf$se)-1)*100

#Sub-annual
full.model3mfS <- rma.mv(yi,vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0imfS,method = "REML")
summary(full.model3mfS)

full.model3mfSS <- rma.mv(yi,vi, 
                          #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0imfSS,method = "REML")
summary(full.model3mfSS)

##FFS fall####

#Annual
full.model3ff <- rma(yi,vi, 
                        #random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iff,method = "REML")
summary(full.model3ff)

full.model3ffS <- rma.mv(yi,vi, 
                         #random = ~ 1 |Treatment, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iffS,method = "REML")
summary(full.model3ffS)

full.model3ffSS <- rma.mv(yi,vi, 
                          #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0iffSS,method = "REML")
summary(full.model3ffSS)

####Overall Effect Sizes Nutrients####

##P flux####

#Total litterfall P flux
full.model3tpf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0itpf,method = "REML")
summary(full.model3tpf)
(exp(full.model3tpf$b)-1)*100
(exp(full.model3tpf$se)-1)*100

#Leaf litterfall P flux
full.model3lpf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilpf,method = "REML")
summary(full.model3lpf)
(exp(full.model3lpf$b)-1)*100
(exp(full.model3lpf$se)-1)*100

#Wood litterfall P flux
full.model3wpf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwpf,method = "REML")
summary(full.model3wpf)
(exp(full.model3wpf$b)-1)*100
(exp(full.model3wpf$se)-1)*100

#Misc P flux
full.model3mpf <- rma.mv(yi,vi,#random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0impf,method = "REML")
summary(full.model3mpf)

##N flux####

#Total Litterfall N flux
full.model3tnf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0itnf,method = "REML")
summary(full.model3tnf)
(exp(full.model3tnf$b)-1)*100
(exp(full.model3tnf$se)-1)*100

#Leaf litterfall N flux
full.model3lnf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilnf,method = "REML")
summary(full.model3lnf)
(exp(full.model3lnf$b)-1)*100
(exp(full.model3lnf$se)-1)*100

#Wood litterfall N flux
full.model3wnf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwnf,method = "REML")
summary(full.model3wnf)
(exp(full.model3wnf$b)-1)*100
(exp(full.model3wnf$se)-1)*100

#Misc litterfall N flux
full.model3mnf <- rma.mv(yi,vi,#random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0imnf,method = "REML")
summary(full.model3mnf)

## P concentration####

#Leaf litterfall P concentration
full.model3lpc <- rma.mv(yi,vi,#random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilpc,method = "REML")
summary(full.model3lpc)
(exp(full.model3lpc$b)-1)*100
(exp(full.model3lpc$se)-1)*100

#Wood litterfall P concentration
full.model3wpc <- rma.mv(yi,vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwpc,method = "REML")
summary(full.model3wpc)
(exp(full.model3wpc$b)-1)*100
(exp(full.model3wpc$se)-1)*100

## N concentration####

#Leaf litterfall N concentration
full.model3lnc <- rma.mv(yi,vi,#random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilnc,method = "REML")
summary(full.model3lnc)
(exp(full.model3lnc$b)-1)*100
(exp(full.model3lnc$se)-1)*100

#Wood litterfall N concentration
full.model3wnc <- rma.mv(yi,vi,#random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwnc,method = "REML")
summary(full.model3wnc)
(exp(full.model3wnc$b)-1)*100
(exp(full.model3wnc$se)-1)*100

#### STEP 4 - Meta-regressions Total Litterfall Mass Flux ####

#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ia$soilP=z.trans(data_es0ia$Other_soil_P)
data_es0ia$hurrwind=z.trans(data_es0ia$HURRECON_wind_ms)
data_es0ia$tsls=z.trans(data_es0ia$YearsSinceLastStorm)

##Case studies with Hurrecon wind data
hurr_sites<-data_es0ia %>% filter(hurrwind!="NA")
str(hurr_sites)#n = 45
hurr_sites$hurrwind=z.trans(hurr_sites$HURRECON_wind_ms)
hurr_sites_no53<-hurr_sites %>% filter(Study_ID!=53)

##Mass and Nutrient Responses comparisons####
#Filtering the data
red_resp_mass<-data_es0ia %>% filter(Site=="Bisley"|Site=="Kokee"|Site=="Guanica"|Site=="Birthday Creek"|Site=="Chamela-Cuixmala"|Site=="Wooroonooran Basalt"|Site=="Wooroonooran Schist"|Site=="Lienhuachi")%>% 
  filter(Treatment=="Ambient")%>% filter(DisturbanceName=="Hugo"|DisturbanceName=="Georges"|DisturbanceName=="Charlie"|DisturbanceName=="Jova"|DisturbanceName=="Larry"|DisturbanceName=="Kalmaegi"|DisturbanceName=="Jangmi")
red_tpf<-data_es0itpf %>% filter(Study_ID!="17")
str(red_tpf)

#Calculating overall responses for both reduced datasets####
full.model_red_tpf <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = red_tpf,method = "REML")
summary(full.model_red_tpf)
tpf_res<-(exp(full.model_red_tpf$b)-1)*100
tpf_res
(exp(full.model_red_tpf$se)-1)*100

full.model_red_mass <- rma.mv(yi,vi,random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                             tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                             data = red_resp_mass,method = "REML")
summary(full.model_red_mass)
tm_res<-(exp(full.model_red_mass$b)-1)*100
tm_res
(exp(full.model_red_mass$se)-1)*100

tnf_res<-(exp(full.model3tnf$b)-1)*100
tnf_res

#Comparing resistance of total P flux to that of total mass flux
tpf_res/tm_res#P flux 1.31 times the 3mass flux

tnf_res/tm_res#N flux 1.1 times the mass flux

## Table 2 - Total Litterfall Response as a function of soil P####

#Table2 - Model1a####
#Multiplying by -1 to represent the resistance
model.mods3full<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = data_es0ia,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(model.mods3full)
tab_model(model.mods3full,p.val="kr")#Significant intercept and regression coefficient
#The intercept represents the mean effect of mean soil P
### Model R^2
(sum(full.model3$sigma2) - sum(model.mods3full$sigma2)) / sum(full.model3$sigma2)


#Testing if removing study 53 changes the effect of soil P on the resistance
model.mods_no53<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = data_no53,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(model.mods_no53)

#Testing if removing four case studies in Hawaii makes a difference
model.mods_red2<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = red_sites2,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(model.mods_red2)

#full model and calculation or R2
fullmodel_red2<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = red_sites2,
                        method = "REML")#using z transformed soil P
summary(fullmodel_red2)
(sum(fullmodel_red2$sigma2) - sum(model.mods_red2$sigma2)) / sum(fullmodel_red2$sigma2)

#Testing alternative random effects
model.mods3full_b<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                        tdist = TRUE,data = data_es0ia,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(model.mods3full_b)

#sites where hurrecon data is available
model.mods3_hurr<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = hurr_sites,
                        method = "REML",
                        mods = ~soilP)#using z transformed soil P
summary(model.mods3_hurr)#similar results

#testing other random effects structure - 
model.mods3full.1<-rma.mv(yi,vi,random = ~ 1 |Region/DisturbanceName, 
                        tdist = TRUE,data = data_es0ia,
                        method = "REML",
                        mods = ~soilP)#using z transformed soil P
summary(model.mods3full.1)

model.mods3full.2<-rma.mv(yi,vi, 
                          random = ~ 1 | DisturbanceName, 
                          tdist = TRUE,data = data_es0ia,
                          method = "REML",
                          mods = ~soilP)#using z transformed soil P
summary(model.mods3full.2)

#AICc to check best model with Site and Cyclone as crossed random effects is the best
BIC(model.mods3full,model.mods3full.1,model.mods3full.2)

#Testing the effect of soil P with a random sample from PR with fewer case studies
summary(new_resistance_data$HURRECON_wind_ms)
new_resistance_data$hurrwind
new_resistance_data<-new_resistance_data%>% filter(HURRECON_wind_ms!="NA")
#z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
#new_resistance_data$soilP=z.trans(new_resistance_data$Other_soil_P)
#new_resistance_data$hurrwind=z.trans(new_resistance_data$HURRECON_wind_ms)
#new_resistance_data$tsls=z.trans(new_resistance_data$YearsSinceLastStorm)

mixed_soilP_redPR<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = new_resistance_data,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(mixed_soilP_redPR)

#Sensitivity analysis - testing the effect of soil P when 5 sites are removed####
#Table S8 - Model 1a
model.mods3full_red<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = red_sites,
                        method = "REML",mods = ~soilP)#using z transformed soil P
summary(model.mods3full_red)
#No change in significance of soil P when sites are deleted

model.mods3full_red2<-rma.mv(yi,vi, 
                            random = ~ 1 |Region/DisturbanceName, 
                            tdist = TRUE,data = red_sites,
                            method = "REML",
                            mods = ~soilP)#using z transformed soil P
summary(model.mods3full_red2)

#Bringing disturbance intensity metrics along with soil P
model_a<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~hurrwind)#standardized HURRECON site-level wind speed
summary(model_a)#Hurrecon wind speed has a positive significant effect on litterfall response

#Table 2 - Model2a####
model_tab2<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*hurrwind)
summary(model_tab2)#Both positive predictors, no significant interaction


#Testing the influence of study 53
model_tab2_no53<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                   data = hurr_sites_no53,method = "REML",
                   mods = ~soilP*hurrwind)
summary(model_tab2_no53)
#checking time since last storm
model_tab2.1<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                   data = hurr_sites,method = "REML",
                   mods = ~soilP+tsls)
summary(model_tab2.1)

model_tab4S<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*hurrwind)
summary(model_tab4S)

#changing random effects structure
model_b2<-rma.mv(yi,vi,random = ~ 1 | Site, tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*windsp*d2track)
summary(model_b2)

model_a3<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                 data = hurr_sites,method = "REML",
                 mods = ~soilP+hurrwind)
summary(model_a3)

model_a.1<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*pred_site_wind)
summary(model_a.1)

model_a1<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP+windsp+d2track)
summary(model_a1)

model_a.0<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                                data = data_es0ia,method = "REML",
                                mods = ~soilP*timesincestorm)
summary(model_a.0)

model_hur<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                tdist = TRUE, 
                data = data_es0ia,
                method = "REML",
                mods = ~soilP*HURRECON_wind_ms)
summary(model_c)

#changing random effects to see if results change
model_c.1<-rma.mv(yi, 
                vi, 
                random = ~ 1|Site, 
                tdist = TRUE, 
                data = data_es0ia,
                method = "REML",
                mods = ~soilP*windsp+d2track)
summary(model_c.1)

AICc(model_c,model_c.1)#Better to include Site/DisturbanceName

#sensitivity analysis on best model - removing sites where available soil P was considered to be 10% of total soil P

levels(data_es0ia$Site)
red_sites <- data_es0ia %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")

red_sites_hurr <- hurr_sites %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")
str(red_sites_hurr)#includes only data points with measured total soil P and HURRECON wind speed

#Table S8 Model 2a####
model_a2_red<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), 
                tdist = TRUE,data = red_sites_hurr,
                method = "REML",
                mods = ~soilP*hurrwind)
summary(model_a2_red)#40 case studies - when the 5 sites are removed, there is a significant positive interaction betwen soil P and wind speed

model_a2_red_b<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                     tdist = TRUE,data = red_sites_hurr,
                     method = "REML",
                     mods = ~soilP*hurrwind)
summary(model_a2_red_b)

#### Leaf fall Response Meta-regression model ####

#standardizing variables 2x sd per Gelman's reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ilf$soilP=z.trans(data_es0ilf$Other_soil_P)
data_es0ilf$hurrwind=z.trans(data_es0ilf$HURRECON_wind_ms)
data_es0ilf$stormfreq=z.trans(data_es0ilf$StormFrequencyNorm)
data_es0ilf$timesincestorm=z.trans(data_es0ilf$YearsSinceLastStorm)

hurr_siteslf<-data_es0ilf %>% filter(hurrwind!="NA")
str(hurr_siteslf)#n = 30

#Leaf fall Table 2####
unique(levels(as.factor(data_es0ilf$Case_study)))

mixed_lf_tab2_1b<-rma.mv(yi,vi,random = list(~ 1|Site,~1|DisturbanceName), 
               tdist = TRUE,data = data_es0ilf,
               method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab2_1b)

mixed_lf_tab2_1b.1<-rma.mv(-1*yi,vi,random = ~ 1|Site/DisturbanceName, 
                         tdist = TRUE,data = data_es0ilf,
                         method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab2_1b.1)

mixed_lf_tab2_1b.2<-rma.mv(-1*yi,vi,random = ~ 1|Region/DisturbanceName, 
                           tdist = TRUE,data = data_es0ilf,
                           method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab2_1b.2)

#Leaf fall Alternative random effects#

mixed_lf_tab4a<-rma.mv(-1*yi,vi,random = ~ 1 | DisturbanceName, 
                  tdist = TRUE,data = data_es0ilf,
                  method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4a)

mixed_lf_tab4b<-rma.mv(-1*yi,vi,random = list(~ 1|Country, ~1|DisturbanceName), 
                 tdist = TRUE,data = data_es0ilf,
                 method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4b)#second best option

mixed_lf_tab4c<-rma.mv(-1*yi,vi,random = list(~ 1|Region, ~1|DisturbanceName), 
                       tdist = TRUE,data = data_es0ilf,
                       method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4c)

mixed_lf_tab4d<-rma.mv(-1*yi,vi,random = ~1|Region, 
                       tdist = TRUE,data = data_es0ilf,
                       method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4d)#BIC is highest

#Table 2 model 1b####
mixed_lf_tab4e<-rma.mv(-1*yi,vi,random = list(~1|Site, ~1|DisturbanceName),
                       tdist = TRUE,data = data_es0ilf,
                       method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4e)#BIC is highest

#comparing leaf fall models based on BIC
BIC(mixed_lf_tab2_1b,mixed_lf_tab2_1b.1,mixed_lf_tab2_1b.2,mixed_lf_tab4a,mixed_lf_tab4b,mixed_lf_tab4c,mixed_lf_tab4d,mixed_lf_tab4e)

#Sensitivity analysis - removing sites where total soil P was estimated from available soil P
red_siteslf <- data_es0ilf %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")
str(red_siteslf)#25
levels(as.factor(red_siteslf$Case_study))

hurr_red_siteslf <- hurr_siteslf %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")
str(red_siteslf)#25 same as red_siteslf 

#Table S8 model 1b####
mixed_lf_S8_1b<-rma.mv(-1*yi,vi,random = list(~1|Site, ~1|DisturbanceName), 
                     tdist = TRUE,data = red_siteslf,
                     method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_S8_1b)#N = 25

#Adding Cyclone metrics

##Table 2 Model 2b - Final Leaf Meta Regression####
mixed_lf_tab2_2b<-rma.mv(-1*yi,vi,random = ~ 1|Site, 
                      tdist = TRUE,data = data_es0ilf,
                      method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab2_2b) 

mixed_lf_tab2_2b.1<-rma.mv(-1*yi,vi,random = list(~ 1|Site, ~1|DisturbanceName), 
                         tdist = TRUE,data = data_es0ilf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab2_2b.1)

#Table S8 model 2b####
mixed_lf_tabS8_2b<-rma.mv(-1*yi,vi,random = list(~ 1|Site, ~1|DisturbanceName), 
                         tdist = TRUE,data = red_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tabS8_2b) 

mixed_lf_tabS8_2b.1<-rma.mv(yi,vi,random = ~ 1|DisturbanceName, 
                          tdist = TRUE,data = red_siteslf,
                          method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tabS8_2b.1) 
str(red_siteslf)

#### Multimodel Inference Meta-Regression ####

summary(model.mods3full_ext1log)

results<-multimodel.inference(TE = "yi", 
                              seTE = "vi",
                              data = data_es0ia,
                              predictors = c("soilP", "windsp", "d2track"),
                              interaction = FALSE,
                              method="REML",
                              test = "knha",eval.criterion="BIC")

summary(results)

plot(results)
results$predictor.importance.plot
names(data_es0ia)
chart.Correlation(data_es0ia[,79:81])

##Checking for multicollinearity

names(dataes0_new)
chart.Correlation(dataes0_new[,c(48,50,79,42)])

names(dataes0_new)
chart.Correlation(data_es0ia[,c(81,50,48)])


####RESPONSE PREDICTION Figures####

##Predictions based on best Meta-Regression model for total litterfall####
#The best model is model_tab4
preds_x2<-predict(model_tab2,levels=0, addx=TRUE)
preds_x2<-data.frame(preds_x2)
str(preds_x2)

#Data frame
metaregplot.1_x2<- cbind(data.frame(hurr_sites$Site, hurr_sites$yi, hurr_sites$vi,preds_x2$pred,preds_x2$se,hurr_sites$soilP,hurr_sites$hurrwind,fac_soilP=factor(hurr_sites$Other_soil_P),hurr_sites$Other_soil_P,hurr_sites$HURRECON_wind_ms))
str(metaregplot.1_x2)
levels(as.factor(hurr_sites$soilP))

##Figure 5a

#Changing color scheme to include unidirectional and log soil P
preg_tot2<-ggplot(metaregplot.1_x2, aes(x=hurr_sites.HURRECON_wind_ms, y=preds_x2.pred))+geom_point(shape=21,aes(col=log(hurr_sites.Other_soil_P),size=preds_x2.se),stroke=1.6)
preg_tot2<-preg_tot2+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_tot2<-preg_tot2+scale_color_gradient(low="#FFFD7D",high="#D8001F")#scale_color_continuous(type="viridis") #FCFF00
preg_tot2
preg_tot2<-preg_tot2+ labs(x="Wind speed (m/s)", y="Predicted resistance")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
preg_tot2<-preg_tot2+guides(size = guide_legend(override.aes = list(col = "black",shape=21)))#,color = guide_legend(override.aes = list(size = 8)))
preg_tot2
Fig5a_new2<- preg_tot2+theme_pubr()+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(legend.position = "top", legend.justification = "center", legend.title = element_text (size = 19), legend.text = element_text (size = 18),axis.title=element_text(size=26),
        axis.text=element_text(size=24))+
  labs(size="",color="Soil P \n(ln mg/kg)")+guides(size=FALSE)+
  annotate("text", x = 5, y = 0.4, label = "a Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#Fig5a_new2####
Fig5a_new2

ggsave(filename = "Fig5a.png",
       plot = Fig5a_new2, width = 8, height = 10, units = 'cm',
       scale = 2, dpi = 1200)

##Leaf fall response to cyclone####

#Predictions
preds_x_leaf<-predict(mixed_lf_tab2_2b.1,levels=0, addx=TRUE) # let's transfer LRR to RR
preds_x_leaf<-data.frame(preds_x_leaf)
str(preds_x_leaf)

#Data frame
metaregplot.1_x_leaf<- cbind(data.frame(data_es0ilf$Site, data_es0ilf$yi, preds_x_leaf$pred,preds_x_leaf$se,data_es0ilf$soilP,fac_soilP=factor(data_es0ilf$Other_soil_P),data_es0ilf$hurrwind,data_es0ilf$Other_soil_P,data_es0ilf$HURRECON_wind_ms))
names(metaregplot.1_x_leaf)

#Figure 5b
preg_lf2<-ggplot(metaregplot.1_x_leaf, aes(x=data_es0ilf.HURRECON_wind_ms, y=preds_x_leaf.pred))+geom_point(shape=21,aes(col=log(data_es0ilf.Other_soil_P),size=preds_x_leaf.se),stroke=1.6)
preg_lf2<-preg_lf2+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_lf2<-preg_lf2+scale_color_gradient(low="#FFFD7D",high="#D8001F")#scale_color_continuous(type="viridis")
preg_lf2
preg_lf2<-preg_lf2+ labs(x="Wind speed (m/s)", y="Predicted resistance")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=data_es0ilf.HURRECON_wind_ms,y=preds_x_leaf$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=data_es0ilf.HURRECON_wind_ms,y=preds_x_leaf$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
preg_lf2<-preg_lf2+guides(size = guide_legend(override.aes = list(col = "black",shape=21)))#,color = guide_legend(override.aes = list(size = 8)))
Fig5b_new2<-preg_lf2+theme_pubr()+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(legend.position = "top",legend.justification = "center",legend.title = element_text (size = 19), legend.text = element_text (size = 20),axis.title=element_text(size=26),axis.text=element_text(size=24))+
  labs(size="",color="Soil P \n(ln mg/kg)")+guides(size=FALSE)+ annotate("text", x = 5, y = 1, label = "b Leaf litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5b_new2####
Fig5b_new2

#Final Fig5a-b with new color gradient
Fig5_Response2<-Fig5a_new2+Fig5b_new2+plot_layout(ncol=2)#+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig5_Response2

#Saving in High Res
ggsave(filename = "Fig5ab_resistance.png",
       plot = Fig5_Response2, width = 18, height = 10, units = 'cm',
       scale = 2, dpi = 1200)


Fig5.2lf_b<-ggplot(data_es0ilf, aes(x=Other_soil_P, y=data_es0ilf$Pre_Mean))+geom_point(aes(color=Other_soil_P),size=6,alpha=0.9)+geom_smooth(method = 'gam', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=22))+ labs(y = bquote('Pre-cyclone leaf fall'~(g/m^2/day)),x = "Total soil P (ln mg/kg)")
rsq_label.alf_b <- paste('R^2 == 0.18')
rsq_label.alf_b
Fig5.2lf_b<-Fig5.2lf_b+scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  annotate("text", x = 3, y = 6, label = rsq_label.alf_b, size=8,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

FigS3b_2<-Fig5.2lf_b+ annotate("text", x = 3, y = 5.3, label = "p<0.02", size=8,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3b_2

lmp<-lmer(Pre_Mean~Other_soil_P+(1|Site),data=data0alf)
summary(lmp)
tab_model(lmp)

FigS5.2<-ggplot(data0alf, aes(x=Other_soil_P, y=Pre_Mean))+geom_point(aes(color=Other_soil_P),size=6,alpha=0.9)+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=22))+ labs(y = bquote('Pre-cyclone Leaf fall'~(g/m^2/day)),x = "Total soil P (mg/kg)")
rsq_label.lf <- paste('R^2 == 0.16')
rsq_label.lf
FigS5.2<-FigS5.2+scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,size="Ln wind speed",color="Soil P")+ annotate("text", x = 3, y = 6, label = rsq_label.lf, size=10,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

FigS5.2


Fig5b<-ggplot(metaregplot, aes(x=data_es0ia.logsoilP, y=fitted.lm_Rp2b.1.))+geom_point(size=6,alpha=0.9,col="black",fill="#2F2747")+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=22),
        axis.text=element_text(size=20))+ labs(x="Log soil P", y="Predicted cyclone response")
Fig5b<-Fig5b+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 22), legend.text = element_text (size = 20),legend.position = "right",legend.justification = c("right"))+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,size="Ln wind speed",color="Ln soil P")#+ annotate("text", x = 0.2, y = 6, label = rsq_label, size=8,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

Fig5b #labs(y = bquote('Density Litterfall mass'~(g/m^2/day))
ggsave(filename = "Fig5b.png",
       plot = Fig5b, width = 13, height = 11, units = 'cm',
       scale = 2, dpi = 600)

Fig5c<-ggplot(metaregplot, aes(x=data_es0ia.logwindsp, y=fitted.lm_Rp2b.1.))+geom_point(size=6,alpha=0.9,col="black",fill="#2F2747")+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=22),
        axis.text=element_text(size=20))+ labs(x="Log wind speed", y="Predicted cyclone response")
Fig5c<-Fig5c+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 22), legend.text = element_text (size = 20),legend.position = "right",legend.justification = c("right"))#+scale_size_continuous(range=c(1,20))+
  #labs(fill=NULL,size="Ln wind speed",color="Ln soil P")#+ annotate("text", x = 0.2, y = 6, label = rsq_label, size=8,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

Fig5c #labs(y = bquote('Density Litterfall mass'~(g/m^2/day))
ggsave(filename = "Fig5c.png",
       plot = Fig5b, width = 13, height = 11, units = 'cm',
       scale = 2, dpi = 600)

Fig5d<-ggplot(metaregplot, aes(x=data_es0ia.logdist, y=fitted.lm_Rp2b.1.))+geom_point(size=6,alpha=0.9,col="black",fill="#2F2747")+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=22),
        axis.text=element_text(size=20))+ labs(x="Log distance to track", y="Predicted cyclone response")
Fig5d<-Fig5d+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 22), legend.text = element_text (size = 20),legend.position = "right",legend.justification = c("right"))+ annotate("text", x = 1, y = 6, label = "d", size=9,hjust=0,colour="black",parse=TRUE)#+scale_size_continuous(range=c(1,20))+
#labs(fill=NULL,size="Ln wind speed",color="Ln soil P")+ annotate("text", x = 1.2, y = 6, label = "d", size=8,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

Fig5d #labs(y = bquote('Density Litterfall mass'~(g/m^2/day))
ggsave(filename = "Fig5d.png",
       plot = Fig5b, width = 13, height = 11, units = 'cm',
       scale = 2, dpi = 600)

FinalFig5 <- #NEW FIG PRE MEAN VS LOG SOIL P+Fig5b+Fig5c+Fig5d+plot_layout(ncol=2)
FinalFig5


library(mgcv)
library(tidymv)
gam_resp2<-gam(dataes0_new.yi~s(dataes0_new.logsoilP,dataes0_new.WMO_wind_kts),data=metaregplot)
summary(gam_resp2)
model_p <- predict_gam(gam_resp2)
model_p

prega<-ggplot(metaregplot, aes(x=dataes0_new.logsoilP, y=norm, color=dataes0_new.logsoilP))+geom_point(aes(size=dataes0_new.Distance_to_Disturb_km),alpha=0.8,fill="white")+geom_smooth(method = 'gam', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=22),
        axis.text=element_text(size=20))+ labs(x="Log soil P", y="Observed Response/Wind speed")
prega
Fig5a<-prega+scale_color_gradient(low="black", high="#EB0526")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 22), legend.text = element_text (size = 20),legend.position = "right",legend.justification = c("right"))+scale_size_continuous(range=c(2,20))+
  labs(fill=NULL,size="Distance to track (km)",color="Log soil P")#+ annotate("text", x = 2, y = 5, label = "R2adj = 0.5 p<0.001", size=8,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

Fig5a

####Leaf Fall MetaRegression and ggplot####

model.mods2lf<-rma.mv(yi, 
                      vi, 
                      random = ~ 1 |Site_ID/Effectsize_ID, 
                      tdist = TRUE, 
                      data = data_es0ilf,
                      method = "REML",
                      mods = ~ Other_soil_P)
summary(model.mods2lf)

lmlf<-lm(data_es0ilf$yi~data_es0ilf$Other_soil_P)
summary(lmlf)
shapiro.test(residuals(lmlf))
lmlf$fitted.values

metaregplotlf<- cbind(data.frame(data_es0ilf$yi, lmlf$fitted.values))
metaregplotlf

preg2lf<-ggplot(metaregplotlf, aes(x=lmlf.fitted.values, y=data_es0ilf.yi))+geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=24))+ labs(x="Predicted Effect Size", y="Observed Effect Size")
preg2lf
library("ggrepel")
preg2lf + geom_label_repel(aes(label=data_es0ilf$Region,
                               fill = factor(data_es0ilf$Region)), color = 'black',
                           size = 4.5) +
  theme(legend.position = "bottom")

####Wood fall MetaRegression####

model.mods2wf<-rma.mv(yi, 
                      vi, 
                      random = ~ 1 |Site_ID, 
                      tdist = TRUE, 
                      data = data_es0iwf,
                      method = "REML",
                      mods = ~ Max_sust_wspeed_mph)
summary(model.mods2wf)

lmwf<-lm(data_es0iwf$yi~data_es0iwf$Max_sust_wspeed_mph+data_es0iwf$Other_soil_P)
summary(lmwf)
shapiro.test(residuals(lmwf))
lmwf$fitted.values

metaregplotwf<- cbind(data.frame(data_es0iwf$yi, lmwf$fitted.values))
metaregplotwf

preg2wf<-ggplot(metaregplotwf, aes(x=lmwf.fitted.values, y=data_es0iwf.yi))+geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=24))+ labs(x="Predicted Effect Size", y="Observed Effect Size")
preg2wf
library("ggrepel")
preg2wf + geom_label_repel(aes(label=data_es0iwf$Region,
                               fill = factor(data_es0iwf$Region)), color = 'black',
                           size = 4.5) +
  theme(legend.position = "bottom")

####Figure 3####

##Fig3a

#### Fig3a Response Case Studies by Region ####
theme_set(theme_bw(base_size=10))

#Data
summmary(data_es0ia)
data_es0ia$Case_study2= paste(data_es0ia$Site, data_es0ia$DisturbanceName,data_es0ia$Treatment,sep="| ")
data_es0ia$Case_study2
data_es0ia$Case_study
data_es0ia_all$Case_study<-c("Cubuy| Georges","Bisley| Georges","East Peak| Georges","East Peak| Georges","East Peak| Georges","East Peak| Georges","Bisley| Hugo","El Verde| Hugo",            
"East Peak| Hugo","Kokee| Iniki","Kokee P+| Iniki","Kokee NP+| Iniki","Kokee N+| Iniki","El Verde| CTE","Utuado| Georges","Guanica| Georges","Milolii| Iniki","Makaha 1| Iniki","Kumuwela| Iniki","Halemanu| Iniki",           
"Wooroonooran Basalt| Larry", "Wooroonooran Schist| Larry", "Mt Spec| Charlie","Mt Spec Disturbed| Charlie","Birthday Creek| Charlie","Birthday Creek| Aivu","Birthday Creek| Ivor","Bisley| Irma",              
"Guanica| Irma","Guayama| Irma","Rio Abajo| Irma","Bisley| Maria","Guanica| Maria","Guayama| Maria","Rio Abajo| Maria","Chamela-Cuixmala| Jova","Chamela-Cuixmala| Patricia","Grande-Terre| Hugo","Lienhuachi| Kalmaegi","Lienhuachi| Fungwong",     
 "Lienhuachi| Sinlaku", "Lienhuachi| Jangmi", "Kengting III| Mindulle", "Kengting IV| Mindulle", "Kengting III| Haima", "Kengting IV| Haima", "Kengting III| Nanmadol", "Kengting IV| Nanmadol", "Gadgarra| Keith")
data_es0ia_all$yi
forrest_data_C<-rbind(data.frame(ES=(-1*data_es0ia_all$yi),SE=sqrt(data_es0ia_all$vi),Type="Case",Site=data_es0ia_all$Site, Case_study=data_es0ia_all$Case_study, Cyclone=data_es0ia_all$DisturbanceName, SoilP=data_es0ia_all$Other_soil_P, Region=data_es0ia_all$Country))
forrest_data_C
forrest_data_C$SoilP=as.numeric(forrest_data_C$SoilP)
forrest_data_C$Case_study=as.factor(forrest_data_C$Case_study)
forrest_data_C$Case2<-factor(forrest_data_C$Case_study, levels=rev(levels(forrest_data_C$Case_study)))
str(forrest_data_C)#aes(x=reorder(Site,-yi), y=yi))

#calculating overall mean and high and low limits of confidence interval
ESm=full.model3$b
ESm
SEm=sqrt(full.model3$se)
ESm+(1.96*SEm)
ESm-(1.96*SEm)

## Colors used for the Regions, respectively: Australia, Guadeloupe, Hawaii, Mexico, Puerto Rico, Taiwan
###+labs(y = bquote('Density Litterfall mass'~(g/m^2/day)),x = "Months since disturbance")+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "b", face="bold",size=8,colour="black")

#setting pantropical response shading
summary(full.model3)
b1 <- list(geom_hline(yintercept = -2.30, color = '#003f5c',alpha=0.4),
           annotate("rect", ymin = -1.39, ymax = -3.21,xmin = -Inf, xmax = Inf,alpha=0.1,linetype=2))
b1

## Fig3a #### Response of Total Litterfall Mass####
plot1D<-ggplot(data=forrest_data_C,aes(x=reorder(Case2,ES),y=ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE),color=Region))+geom_pointrange(alpha=0.9,size=1.2)+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
plot2D<-plot1D+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1.2, alpha = 0.6)#+geom_hline(aes(yintercept=3.64), lty=2, colour = "#003f5c", cex=0.9, alpha = 0.4)
plot3D<-plot2D+xlab("Case study")+theme_bw()+ylab("Resistance")#ylab(expression(Response~ln~(litterfall[ti]/litterfall[t0]))) #ylab(expression(Response~(ln~Total~litterfall~t[i]~t[0]^-1)))
finalFig3a<-plot3D+theme(axis.title=element_text(size=28),axis.text.x=element_text(size=22,angle=0, hjust=0.5),axis.text.y=element_text(size=17,angle=0, hjust=1))+theme(legend.title = element_blank(), legend.text = element_text (size = 24),legend.position=c(0.3,0.92),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))                                                                                                                                                                        legend.box.background = element_rect(colour = "grey"))+b1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))
finalFig3a<-finalFig3a+ annotate("text", y = -7, x = 46, fontface="bold",label = "a", size=8,colour="black")
finalFig3a

#Saving in high res
ggsave(filename = "Fig3a_Response_New.png",
       plot = finalFig3a, width = 15, height = 18, units = 'cm',
       scale = 2, dpi = 1000)

####Fig3b-d####

##Fig3b Data
data_frac4 <- rbind(data.frame(group="FFS", variable="Annual", estimate=-1*full.model3ff$b, var="Mass flux",
                               ci_low=(-1*full.model3ff$b-(1.96*full.model3ff$se)),ci_up=(-1*full.model3ff$b+(1.96*full.model3ff$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Wood", variable="Annual",estimate=-1*full.model3wf$b, var="Mass flux",ci_low=(-1*full.model3wf$b-(1.96*full.model3wf$se)),ci_up=(-1*full.model3wf$b+(1.96*full.model3wf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Leaf", variable="Annual",estimate=-1*full.model3lf$b, var="Mass flux",ci_low=(-1*full.model3lf$b-(1.96*full.model3lf$se)),ci_up=(-1*full.model3lf$b+(1.96*full.model3lf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Misc.", variable="Annual",estimate=-1*full.model3mf$b,var="Mass flux",ci_low=(-1*full.model3mf$b-(1.96*full.model3mf$se)),ci_up=(-1*full.model3mf$b+(1.96*full.model3mf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE))
data_frac4

## Fig3c Data
data_impact_frac2 <- rbind(data.frame(group="Total", var="P flux",estimate2=-1*full.model3tpf$b, ci_low2=(-1*full.model3tpf$b-(1.96*full.model3tpf$se)),ci_up2=(-1*full.model3tpf$b+(1.96*full.model3tpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Wood", var="P flux", estimate2=-1*full.model3wpf$b,
                                      ci_low2=(-1*full.model3wpf$b-(1.96*full.model3wpf$se)), ci_up2=(-1*full.model3wpf$b+(1.96*full.model3wpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Leaf", var="P flux", estimate2=-1*full.model3lpf$b,ci_low2=(-1*full.model3lpf$b-(1.96*full.model3lpf$se)), ci_up2=(-1*full.model3lpf$b+(1.96*full.model3lpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Misc.", var="P flux", estimate2=-1*full.model3mpf$b, ci_low2=(-1*full.model3mpf$b-(1.96*full.model3mpf$se)), ci_up2=(-1*full.model3mpf$b+(1.96*full.model3mpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Total", var="N flux",estimate2=-1*full.model3tnf$b,ci_low2=(-1*full.model3tnf$b-(1.96*full.model3tnf$se)), ci_up2=(-1*full.model3tnf$b+(1.96*full.model3tnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Wood", var="N flux", estimate2=-1*full.model3wnf$b, ci_low2=(-1*full.model3wnf$b-(1.96*full.model3wnf$se)), ci_up2=(-1*full.model3wnf$b+(1.96*full.model3wnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Leaf", var="N flux", estimate2=-1*full.model3lnf$b,ci_low2=(-1*full.model3lnf$b-(1.96*full.model3lnf$se)), ci_up2=(-1*full.model3lnf$b+(1.96*full.model3lnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Misc.", var="N flux", estimate2=-1*full.model3mnf$b,ci_low2=(-1*full.model3mnf$b-(1.96*full.model3mnf$se)), ci_up2=(-1*full.model3mnf$b+(1.96*full.model3mnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE))
data_impact_frac2

##Fig3d Data
data_impact_fracPN <- rbind(data.frame(group="Wood", var="P concentration", estimate2=-1*full.model3wpc$b,ci_low2=(-1*full.model3wpc$b-(1.96*full.model3wpc$se)),ci_up2=(-1*full.model3wpc$b+(1.96*full.model3wpc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Leaf", var="P concentration", estimate2=-1*full.model3lpc$b,ci_low2=(-1*full.model3lpc$b-(1.96*full.model3lpc$se)),ci_up2=(-1*full.model3lpc$b+(1.96*full.model3lpc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Wood", var="N concentration", estimate2=-1*full.model3wnc$b,
                                       ci_low2=(-1*full.model3wnc$b-(1.96*full.model3wnc$se)),ci_up2=(-1*full.model3wnc$b+(1.96*full.model3wnc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Leaf", var="N concentration", estimate2=-1*full.model3lnc$b,
                                       ci_low2=(-1*full.model3lnc$b-(1.96*full.model3lnc$se)),ci_up2=(-1*full.model3lnc$b+(1.96*full.model3lnc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE))

data_impact_fracPN
summary(full.model3lnc)

##Fig3b
pfrac4<-ggplot(data_frac4, aes(x=group,y=estimate,ymax=ci_up,ymin=ci_low,shape=var))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6))+
  geom_pointrange(size=1.2, alpha=0.8,position=position_dodge(width=c(0.7, 1.2)))+#coord_flip()+
  geom_hline(aes(yintercept=0), lty=2,size=1,col="magenta",alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="", x="") +#scale_shape_discrete(solid=F)+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),axis.title=element_text(size=20),
        axis.text=element_text(size=24),legend.text =  element_text(size=24),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.box="horizontal",legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.position = c(0.84,0.9),legend.box.background = element_blank())
pfrac4#+ scale_color_grey(start=0.65, end=0.25)
Fig3b<-pfrac4+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ annotate("text", y = 1.4, x = 0.5, fontface="bold",label = "b", size=8,colour="black")+ annotate("text", y = 0.4, x = 1.16, fontface="bold",label = "(11)", size=6,colour="black")+ annotate("text", y = 0.4, x = 2, fontface="bold",label = "(29)", size=6,colour="black")+ annotate("text", y = 0.4, x = 3, fontface="bold",label = "(30)", size=6,colour="black")+ annotate("text", y = 0.4, x =4.13, fontface="bold",label = "(9)", size=6,colour="black")
Fig3b

### Fig3c
data_impact_frac2
Fig3c<-ggplot(data_impact_frac2, aes(x=group,y=estimate2,ymax=ci_up2,ymin=ci_low2, shape = var,col=var))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6))+
  geom_pointrange(mapping=aes(shape=var),size=1.2, stroke=1.2,position=position_dodge(width=c(0.4, 0.8)))+
  geom_hline(aes(yintercept=0), lty=2,size=1.2,col="magenta", alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="Pantropical resistance", x="") +scale_shape_discrete(solid=F)+ scale_color_manual(values=c("#167923","#1C39A8"))+
  theme_bw() +# ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = 0.5),axis.text.x =element_text(vjust = 0.3),
        axis.title.y =element_text(vjust = 0.5,size=26),
        axis.text=element_text(size=24),legend.box="horizontal",legend.text =  element_text(size=24),legend.title = element_blank(),legend.position = c(.86,.2),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Fig3c<-Fig3c+ annotate("text", y = 1.2, x = 0.5, fontface="bold",label = "c", size=8,colour="black")+ annotate("text", y = -0.4, x = 0.8, fontface="bold",label = "(10)", size=6,colour="black")+ annotate("text", y = -0.4, x = 1.16, fontface="bold",label = "(9)", size=6,colour="black")+
  annotate("text", y = -0.4, x = 1.8, fontface="bold",label = "(8)", size=6,colour="black")+annotate("text", y = -0.4, x = 2.16, fontface="bold",label = "(7)", size=6,colour="black")+ annotate("text", y = -0.4, x = 2.8, fontface="bold",label = "(9)", size=6,colour="black")+annotate("text", y = -0.4, x = 3.16, fontface="bold",label = "(9)", size=6,colour="black")+ annotate("text", y = -0.4, x =3.8, fontface="bold",label = "(3)", size=6,colour="black")+ annotate("text", y = -0.4, x =4.16, fontface="bold",label = "(3)", size=6,colour="black")
Fig3c

##Fig3d
Fig3d<-ggplot(data_impact_fracPN, aes(x=group,y=estimate2,ymax=ci_up2,ymin=ci_low2, shape = var,col=var))+scale_y_continuous(breaks=c(0,-0.4,-0.8,-1.2))+
  geom_pointrange(mapping=aes(shape=var),size=1.2, position=position_dodge(width=c(0.3, 0.6)))+#coord_flip()+
  scale_shape_discrete(solid=F)+
  geom_hline(aes(yintercept=0), lty=2,size=1.2,col="magenta", alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="", x="Litterfall fraction") +scale_color_manual(values=c("#167923","#1C39A8"))+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),
        axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),
        axis.title=element_text(size=26),
        axis.text=element_text(size=24),legend.box="horizontal",legend.text =  element_text(size=24),legend.title = element_blank(),legend.position = c(.74,0.2),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Fig3d<-Fig3d+ annotate("text", y = 0.2, x = 0.45, fontface="bold",label = "d", size=8,colour="black")+ annotate("text", y = 0.06, x = 0.9, fontface="bold",label = "(3)", size=6,colour="black")+ annotate("text", y = 0.06, x = 1.09, fontface="bold",label = "(3)", size=6,colour="black")+
  annotate("text", y = 0.06, x = 1.85, fontface="bold",label = "(10)", size=6,colour="black")+annotate("text", y = 0.06, x = 2.09, fontface="bold",label = "(10)", size=6,colour="black")
Fig3d

#Figure 3b-d
FinalFig3 <- Fig3b+Fig3c+Fig3d+plot_layout(ncol=1)
FinalFig3

#Final Figure 3a-d####
Fig3<-finalFig3a+FinalFig3+plot_layout(ncol=2)
Fig3

#Saving in high res the Final Figure 3a-d
ggsave(filename = "Fig3ad_Final_Resistance2.png",plot = Fig3, width = 22, height = 22, units = 'cm',scale = 2, dpi = 1000)

####Mass and Nutrient Flux Fractions by % change####

data_impact_frac <- rbind(data.frame(group="Total", var="Mass flux",estimate=(exp(full.model3$b)-1)*100,estimate2=full.model3$b,
                                     ci_low=(exp(full.model3$ci.lb)-1)*100, ci_low2=full.model3$ci.lb,ci_up=(exp(full.model3$ci.ub)-1)*100,ci_up2=full.model3$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE), 
                          data.frame(group="FFS", var="Mass flux", estimate=-((1-exp(full.model3ff$b))*100),estimate2=full.model3ff$b,
                                     ci_low=-((1-exp(full.model3ff$ci.lb))*100), ci_low2=full.model3ff$ci.lb,ci_up=(exp(full.model3ff$ci.ub)-1)*100,ci_up2=full.model3ff$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Misc.", var="Mass flux", estimate=(exp(full.model3mf$b)-1)*100,estimate2=full.model3mf$b,
                                     ci_low=(exp(full.model3mf$ci.lb)-1)*100, ci_low2= full.model3mf$ci.lb,ci_up=(exp(full.model3mf$ci.ub)-1)*100,ci_up2=full.model3mf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Leaf", var="Mass flux", estimate=(exp(full.model3lf$b)-1)*100,estimate2=full.model3lf$b,
                                     ci_low=(exp(full.model3lf$ci.lb)-1)*100, ci_low2=full.model3lf$ci.lb,ci_up=(exp(full.model3lf$ci.ub)-1)*100,ci_up2=full.model3lf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Wood", var="Mass flux", estimate=(exp(full.model3wf$b)-1)*100, estimate2=full.model3wf$b,
                                     ci_low=(exp(full.model3wf$ci.lb)-1)*100, ci_low2=full.model3wf$ci.lb, ci_up=(exp(full.model3wf$ci.ub)-1)*100,ci_up2=full.model3wf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Total", var="P flux",estimate=(exp(full.model3tpf$b)-1)*100, estimate2=full.model3tpf$b,
                                     ci_low=(exp(full.model3tpf$ci.lb)-1)*100, ci_low2=full.model3tpf$ci.lb,ci_up=(exp(full.model3tpf$ci.ub)-1)*100,ci_up2=full.model3tpf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Misc.", var="P flux", estimate=-((1-exp(full.model3mpf$b))*100), estimate2=full.model3mpf$b,
                                     ci_low=-((1-exp(full.model3mpf$ci.lb))*100), ci_low2=full.model3mpf$ci.lb, ci_up=(exp(full.model3mpf$ci.ub)-1)*100,ci_up2=full.model3mpf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Leaf", var="P flux", estimate=(exp(full.model3lpf$b)-1)*100,estimate2=full.model3lpf$b, 
                                     ci_low=-((1-exp(full.model3lpf$ci.lb))*100), ci_low2=full.model3lpf$ci.lb, ci_up=(exp(full.model3lpf$ci.ub)-1)*100,ci_up2=full.model3lpf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Wood", var="P flux", estimate=(exp(full.model3wpf$b)-1)*100, estimate2=full.model3wpf$b,
                                     ci_low=(exp(full.model3wpf$ci.lb)-1)*100, ci_low2=full.model3wpf$ci.lb, ci_up=(exp(full.model3wpf$ci.ub)-1)*100,ci_up2=full.model3wpf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Total", var="N flux",estimate=(exp(full.model3tnf$b)-1)*100, estimate2=full.model3tnf$b,
                                     ci_low=(exp(full.model3tnf$ci.lb)-1)*100, ci_low2=full.model3tnf$ci.lb, ci_up=(exp(full.model3tnf$ci.ub)-1)*100,ci_up2=full.model3tnf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Misc.", var="N flux", estimate=-((1-exp(full.model3mnf$b))*100), estimate2=full.model3mnf$b,
                                     ci_low=-((1-exp(full.model3mnf$b))*100), ci_low2=full.model3mnf$ci.lb, ci_up=(exp(full.model3mnf$ci.ub)-1)*100,ci_up2=full.model3mnf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Leaf", var="N flux", estimate=(exp(full.model3lnf$b)-1)*100, estimate2=full.model3lnf$b,
                                     ci_low=(exp(full.model3lnf$ci.lb)-1)*100, ci_low2=full.model3lnf$ci.lb, ci_up=(exp(full.model3lnf$ci.ub)-1)*100,ci_up2=full.model3lnf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE),
                          data.frame(group="Wood", var="N flux", estimate=(exp(full.model3wnf$b)-1)*100, estimate2=full.model3wnf$b,
                                     ci_low=(exp(full.model3wnf$ci.lb)-1)*100, ci_low2 = full.model3wnf$ci.lb, ci_up=(exp(full.model3wnf$ci.ub)-1)*100, ci_up2=full.model3wnf$ci.ub,
                                     row.names=FALSE, stringsAsFactors=TRUE))
data_impact_frac

ggplot(model_dat_Amb_frac, aes(x=group,y=estimate,ymax=ci_up,ymin=ci_low)) +
  geom_pointrange(mapping=aes(color=group),size=1.2, alpha=0.8) +
  scale_color_manual(values=c("black","black","black", "black", "black"))+ylim(-100,300)+
  #scale_fill_manual(values=c("brown","black","green", "purple", "gray"))+
  #scale_colour_discrete(name="Method of enrichment")+
  geom_hline(aes(yintercept=0), lty=2,size=1.5,col="magenta") + # this adds a dotted line for effect size of 0
  labs(y="Summary effect size (%)", x="Litterfall fraction") +
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),
        axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),legend.text =  element_blank(),legend.title = element_blank())

####RESPONSE GGPLOT FIGURE N AND P CONCENTRATIONS####

names(nutmeta)

data_es0ilpc$yi

forrest_data_D<-rbind(data.frame(ES=full.model3$yi,SE=sqrt(full.model3$vi),Type="Case",Element="P",Fraction="Leaf",Site=data_es0ilpc$Site, Case_study=data_es0ilpc$Effectsize_ID, Cyclone=data_es0ilpc$Disturbance, Region=data_es0ilpc$Area),
                      data.frame(ES=full.model3lnc$yi,SE=sqrt(full.model3lnc$vi),Type="Case",Element="N",Fraction="Leaf",Site=data_es0ilnc$Site, Case_study=data_es0ilnc$Effectsize_ID, Cyclone=data_es0ilnc$Disturbance, Region=data_es0ilnc$Area),
                      data.frame(ES=full.model3wpc$yi,SE=sqrt(full.model3wpc$vi),Type="Case",Element="P",Fraction="Wood",Site=data_es0iwpc$Site, Case_study=data_es0iwpc$Effectsize_ID, Cyclone=data_es0iwpc$Disturbance, Region=data_es0iwpc$Area),
                      data.frame(ES=full.model3wnc$yi,SE=sqrt(full.model3wnc$vi),Type="Case",Element="N",Fraction="Wood",Site=data_es0iwnc$Site, Case_study=data_es0iwnc$Effectsize_ID, Cyclone=data_es0iwnc$Disturbance, Region=data_es0iwnc$Area))
forrest_data_D


forrest_data_C$SoilP=as.numeric(forrest_data_C$SoilP)
forrest_data_C$Case_study=as.factor(forrest_data_C$Case_study)
forrest_data_C$Case2<-factor(forrest_data_C$Case_study, levels=rev(levels(forrest_data_C$Case_study)))
str(forrest_data_C)
levels(forrest_data_C$Site2)#aes(x=reorder(Site,-yi), y=yi))

plot1C<-ggplot(data=forrest_data_C,aes(x=reorder(Case2,-ES),y=ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE)))+geom_pointrange(aes(color=Region),alpha=0.8)+guides(title="Region")

plot2C<-plot1C+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = 0.6)+geom_hline(aes(yintercept=3.88), lty=2, color = "black", cex=0.9, alpha = 0.6)

plot3C<-plot2C+xlab("Case study")+ylab("Response [ln(litterfall ti/t0)]")+theme_bw()                                                                                      #color="black")

plot3C+theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17,angle=0, hjust=1))+theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.85,0.9),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+b1+scale_color_discrete("Dark2")


#example to included shaded area on background of figure

b1 <- list(geom_vline(xintercept = 3.85, color = 'gray',alpha=0.1),
           annotate("rect", alpha = .1,ymin = 2.58, ymax = 5.12,xmin = -Inf, xmax = Inf))

modelplot(mod, background = b)

b <- list(geom_vline(xintercept = 3.85, color = 'black'),
          annotate("rect", alpha = .1,
                   xmin = -2.58, xmax = 5.12, 
                   ymin = -Inf, ymax = Inf),
          geom_point(aes(y = term, x = estimate), alpha = .3, 
                     size = 10, color = 'red', shape = 'square'))

modelplot(mod, background = b)

##END##
