###Meta-Analysis of cyclone response of forest litterfall across the tropics###

library(styler)
library(lintr)
library(magrittr)
library(ggcorrplot)
library(tidymv)
library(LMERConvenienceFunctions)
library(relaimpo)
library(lavaan)
library(eqs2lavaan)
library(ggplot2)
library(ggstatsplot)
library(nlme)
library(gapminder)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(magrittr)
library(mgcv)
library(nlme)
library(dplyr)
library(hrbrthemes)
library(PerformanceAnalytics)
library(patchwork)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(lattice)
library(lmerTest)
library(aod)
library(corrplot)
library(modelsummary)
library(tidyr)
library(plyr)
library(fs)
library(readr)
library(dygraphs)
library(xts)
library(plot3D)
library(forcats)
library(ggridges)
library("forcats")
library("broom")
library(lmerTest)
library(metafor)
library(forestmodel)
library(metaforest)
library(caret)
library(ranger)
library(dmetar)
library(MuMIn)
library(rJava)
library(leaps)
library(glmulti)
library(metaviz)
library(Matrix)
library(ggExtra)
library(RColorBrewer)
library(scales)
library(tibble)
library(AICcmodavg)
install.packages("devtools")

####STEP 0 uploading data####
#litterfall mass data
metadat<-read.csv(file.choose())#20210520_Litterfall_Mass
attach(metadat)
str(metadat)#2370 obs of 77 variables

#transforming as numeric
metadat$HURRECON_wind_ms=as.numeric(metadat$HURRECON_wind_ms)
metadat$Gale_wind_duration_minutes=as.numeric(metadat$Gale_wind_duration_minutes)

#create Case study column
metadat$Case_study= paste(metadat$Site, metadat$DisturbanceName,sep="| ")
unique(levels(as.factor(metadat$Case_study)))

#Nutrient flux and concentration data
nutmeta<-read.csv(file.choose())#20210520_Litterfall_Nutrients
attach(nutmeta)
str(nutmeta)

#create Case study column
nutmeta$Case_study= paste(nutmeta$Site, nutmeta$Disturbance, sep="| ")
levels(as.factor(nutmeta$Case_study))

####STEP 1 Data Wrangling ####

###Litterfall Mass flux##

##Total Litterfall mass####
#Annual
data0a<-metadat %>% filter(Fraction=="TotLitfall")%>%filter(Cat_TSD_months=="0-0.5")
str(data0a)#49 observations
levels(as.factor(data0a$Case_study))
#Sub-annual
data0aS<- data0a%>% filter(Pre_Mean_MonthSpecific!="NA")#to exclude a factor level
str(data0aS)#23observations
data0aS$Effectsize_ID
#Annual with same observations as in sub-annual
data0aSS<-data0a %>% filter(Effectsize_ID == "34"|Effectsize_ID == "44"|Effectsize_ID == "47"|Effectsize_ID == "1015"|
                              Effectsize_ID == "1738"|Effectsize_ID == "1848"|Effectsize_ID == "2102"|Effectsize_ID == "2106"|Effectsize_ID == "2110"|Effectsize_ID == "2114"|
                              Effectsize_ID == "2118"|Effectsize_ID == "2122"|Effectsize_ID=="2126"|Effectsize_ID == "2130"|Effectsize_ID == "2134"|Effectsize_ID == "2143"|Effectsize_ID == "2346"|Effectsize_ID == "2347"|Effectsize_ID == "2348"|Effectsize_ID == "2349"|Effectsize_ID == "2350"|Effectsize_ID == "2351"|Effectsize_ID == "2352")
str(data0aSS)#23 observations including the same effect sizes with and without subannual pre and post data

##Leaf fall####

#Annual
data0alf<-metadat %>% filter(Fraction=="Leaf fall")%>%filter(Cat_TSD_months=="0-0.5")
str(data0alf)#31 observations
#Sub-annual
data0alfS<- data0alf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0alfS)#17 observations
data0alfS$Effectsize_ID
#Annual with same observations as in sub-annual
data0alfSS<-data0alf %>% filter(Effectsize_ID == "40"|Effectsize_ID == "45"|Effectsize_ID == "1076"|Effectsize_ID == "1711"|Effectsize_ID == "1760"|Effectsize_ID == "1870"|Effectsize_ID == "2103"|Effectsize_ID == "2107"|Effectsize_ID == "2111"|Effectsize_ID == "2115"|Effectsize_ID=="2119"|Effectsize_ID=="2123"|Effectsize_ID=="2127"|Effectsize_ID == "2131"|Effectsize_ID == "2152"|Effectsize_ID == "2154"|Effectsize_ID == "2170")
str(data0alfSS)#17 obs

##Wood fall####

#Annual
data0awf<-metadat %>% filter(Fraction=="Wood fall")%>%filter(Cat_TSD_months=="0-0.5")
str(data0awf)#30 observations
#Sub-annual
data0awfS<- data0awf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0awfS)#15 observations
data0awfS$Effectsize_ID
#Annual with same observations as in sub-annua
data0awfSS<-data0awf %>% filter(Effectsize_ID == "41"|Effectsize_ID == "42"|Effectsize_ID == "1197"|Effectsize_ID == "1782"|Effectsize_ID == "1892"|Effectsize_ID == "2104"|Effectsize_ID == "2108"|Effectsize_ID == "2112"|Effectsize_ID == "2116"|Effectsize_ID == "2120"|Effectsize_ID == "2124"|Effectsize_ID == "2128"|Effectsize_ID == "2132"|Effectsize_ID == "2156"|Effectsize_ID == "2158")
str(data0awfSS)#15 obs

## FFS fall ####
data0aff<-metadat %>% filter(Fraction=="FFS fall")%>%filter(Cat_TSD_months=="0-0.5")
str(data0aff)#14 observations

data0affS<-data0aff %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0affS)#10 observations
data0affS$Effectsize_ID

data0affSS<-data0aff %>% filter(Effectsize_ID == "1804"|Effectsize_ID == "1914"|Effectsize_ID == "2105"|Effectsize_ID == "2109"|Effectsize_ID == "2113"|Effectsize_ID == "2117"|Effectsize_ID == "2117"|Effectsize_ID == "2121"|Effectsize_ID == "2125"|Effectsize_ID == "2129"|Effectsize_ID == "2133")
str(data0affSS)#10 obs

## Miscellaneous fall ####
data0amf<-metadat %>% filter(Fraction=="Misc fall")%>%filter(Cat_TSD_months=="0-0.5")
str(data0amf)#9 observations

data0amfSS<-data0amf %>% filter(Effectsize_ID == "43"|Effectsize_ID == "46"|Effectsize_ID == "1826"|Effectsize_ID == "1936")
str(data0amfSS)#4 obs

data0amfS<-data0amf %>% filter(Pre_Mean_MonthSpecific!="NA")
str(data0amfS)#4 obs

####Nutrient fluxes####

####P flux####
str(nutmeta)

#Total P flux
data0tpf<-nutmeta %>%filter(Fraction=="TotLitfall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0tpf)#11
#Leaf fall
data0lpf<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0lpf)#7
#Wood fall
data0wpf<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0wpf)#7
#FFS fall
data0ffspf<-nutmeta %>%filter(Fraction=="FFS fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0ffspf)#1
#Miscellaneous fall
data0mpf<-nutmeta %>%filter(Fraction=="Misc fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0mpf)#2

#### N flux####

#Total litterfall
data0tnf<-nutmeta %>%filter(Fraction=="TotLitfall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0tnf)#10
#Leaf fall
data0lnf<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0lnf)#7
#Wood fall
data0wnf<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0wnf)#6
#FFS fall
data0ffsnf<-nutmeta %>%filter(Fraction=="FFS fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0ffsnf)#1
#Miscellaneous fall
data0mnf<-nutmeta %>%filter(Fraction=="Misc fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/m2/day")%>%filter(Cat_TSD_months=="0-0.5")
str(data0mnf)#2

####P concentration####

#Leaf fall
data0lpc<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")
str(data0lpc)#10
#Wood fall
data0wpc<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Variable=="P") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")
str(data0wpc)#3

##N concentration ####

#Leaf fall
data0lnc<-nutmeta %>%filter(Fraction=="Leaf fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")
str(data0lnc)#10
#Wood fall
data0wnc<-nutmeta %>%filter(Fraction=="Wood fall")%>%filter(Variable=="N") %>%filter(Raw_Unit=="mg/g")%>%filter(Cat_TSD_months=="0-0.5")
str(data0wnc)#3

####STEP 2 Individual Effect size calculation####

##Total Litterfall##
str(data0a)#34 observations

#Annual
data_es0ia <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = data0a, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).
#Sub-annual data
data_es0iaS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                      sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = data0aS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)
str(data_es0iaS)
#long-term pre mean for sites that also have subannual data
data_es0iaSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0aSS, measure = "ROM") # ROM =  log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011).str(data_es0iaS)

#### Litterfall fractions ####
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
data_es0iff$yi
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

##P flux####

#Total and fractions
data_es0itpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0tpf, measure = "ROM")
data_es0ilpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lpf, measure = "ROM")
data_es0iwpf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wpf, measure = "ROM")
data_es0iffspf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = data0ffspf, measure = "ROM")
data_es0impf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0mpf, measure = "ROM")

##P concentration####
#Leaf fall
data_es0ilpc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0lpc, measure = "ROM")
#Wood fall
data_es0iwpc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0wpc, measure = "ROM")

##N flux####
#Total and fractions
data_es0itnf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0tnf, measure = "ROM")
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
data_es0imnc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0mnc, measure = "ROM")
data_es0itcf <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = data0tcf, measure = "ROM")

#### STEP 3 Overall Effect Size calculations ####

## Total Litterfall ##
#Data viz
df_TSD_tot_resp %>% ggplot(aes(x = yi, fill=DisturbanceName)) +geom_histogram(binwidth = 0.5)+theme_pubr()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18),axis.text=element_text(size=16),axis.title=element_text(size=20),
        legend.position="right",legend.text =  element_text(size=16,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.box.background = element_rect(colour = "black"))+
  #scale_fill_brewer(palette="Paired")+
  labs(x="Response effect sizes",y="Number of observations")+ annotate("text", x = 4, y = 8, label = "Response of total litterfall",size=6,colour="black")

##Final model for Pantropical response of total litterfall mass####
full.model3 <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                       tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                       data = data_es0ia,method = "REML")

summary(full.model3)
((exp(full.model3$b)-1)*100)
(exp(full.model3$se)-1)*100

##Sensitivity analysis - remove CTE
no_cte_totlitfall<-data_es0ia %>% filter(DisturbanceName!="CTE")

full.model3_red <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = no_cte_totlitfall,method = "REML")

summary(full.model3_red)
((exp(full.model3_red$b)-1)*100)
(exp(full.model3_red$se)-1)*100

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

red_sites <- data_es0ia %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")

full.model3_red<- rma.mv(yi, 
                         vi, 
                         random = list(~ 1 | Site,~1|DisturbanceName),
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = red_sites,
                         method = "REML")
summary(full.model3_red)

dataes0_new2

##Calculating Overall Effect Sizes for sub-annual data

full.model3S <- rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName),
                       tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                       data = data_es0iaS,
                       method = "REML")
summary(full.model3S)
(exp(full.model3S$b)-1)*100
(exp(full.model3S$se)-1)*100

levels(as.factor(data_es0iaSS$Case_study))
full.model3SS <- rma.mv(yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iaSS,
                        method = "REML")
summary(full.model3SS)
(exp(full.model3SS$b)-1)*100
(exp(full.model3SS$se)-1)*100

2837.5-2835.206

## Leaf fall ####
data_es0ilf$Case_study
full.model3lf <- rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0ilf,
                        method = "REML")
summary(full.model3lf)
(exp(full.model3lf$b)-1)*100
(exp(full.model3lf$se)-1)*100

levels(as.factor(data_es0ilfS$Case_study))

full.model3lfS <- rma.mv(yi, 
                         vi, 
                         random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilfS,
                         method = "REML")
summary(full.model3lfS)

full.model3lfSS <- rma.mv(yi, 
                          vi, 
                          random =  list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0ilfSS,
                          method = "REML")
summary(full.model3lfSS)

###Wood fall####
levels(as.factor(data_es0iwf$Case_study))
full.model3wf <- rma.mv(yi, 
                        vi, 
                        random =  list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iwf,
                        method = "REML")
summary(full.model3wf)
(exp(full.model3wf$b)-1)*100
(exp(full.model3wf$se)-1)*100

full.model3wfS <- rma.mv(yi, 
                         vi, 
                         random =list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwfS,
                         method = "REML")
summary(full.model3wfS)

full.model3wfSS <- rma.mv(yi, 
                          vi, 
                          random =list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0iwfSS,
                          method = "REML")
summary(full.model3wfSS)

### Miscellaneous fall
data_es0imf$Case_study
full.model3mf <- rma.mv(yi, 
                        vi, 
                        random = list(~ 1 | Site,~1|DisturbanceName), #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0imf,
                        method = "REML")
summary(full.model3mf)
(exp(full.model3mf$b)-1)*100
(exp(full.model3mf$se)-1)*100

data_es0imfS$Case_study
full.model3mfS <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0imfS,
                         method = "REML")
summary(full.model3mfS)

full.model3mfSS <- rma.mv(yi, 
                          vi, 
                          #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0imfSS,
                          method = "REML")
summary(full.model3mfSS)

##FFS fall
data_es0iff$Case_study
full.model3ff <- rma(yi, 
                        vi, 
                        #random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                        tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                        data = data_es0iff,
                        method = "REML")
summary(full.model3ff)

full.model3ffS <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 |Treatment, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iffS,
                         method = "REML")
summary(full.model3ffS)

full.model3ffSS <- rma.mv(yi, 
                          vi, 
                          #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                          tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                          data = data_es0iffSS,
                          method = "REML")
summary(full.model3ffSS)

####Overall Effect Sizes Nutrients####

##P flux
full.model3tpf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0itpf,
                         method = "REML")
summary(full.model3tpf)
(exp(full.model3tpf$b)-1)*100
(exp(full.model3tpf$se)-1)*100


full.model3lpf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilpf,
                         method = "REML")
summary(full.model3lpf)
(exp(full.model3lpf$b)-1)*100
(exp(full.model3lpf$se)-1)*100

full.model3wpf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwpf,
                         method = "REML")
summary(full.model3wpf)
(exp(full.model3wpf$b)-1)*100
(exp(full.model3wpf$se)-1)*100

full.model3mpf <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0impf,
                         method = "REML")
summary(full.model3mpf)

## P concentration
full.model3lpc <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilpc,
                         method = "REML")
summary(full.model3lpc)
tab_model(full.model3lpc)
(exp(full.model3lpc$b)-1)*100
(exp(full.model3lpc$se)-1)*100

full.model3wpc <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwpc,
                         method = "REML")
summary(full.model3wpc)
(exp(full.model3wpc$b)-1)*100
(exp(full.model3wpc$se)-1)*100

## N flux 
full.model3tnf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0itnf,
                         method = "REML")
summary(full.model3tnf)
(exp(full.model3tnf$b)-1)*100
(exp(full.model3tnf$se)-1)*100

full.model3lnf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilnf,
                         method = "REML")
summary(full.model3lnf)
(exp(full.model3lnf$b)-1)*100
(exp(full.model3lnf$se)-1)*100

full.model3wnf <- rma.mv(yi, 
                         vi, 
                         random = ~ 1 | Site, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwnf,
                         method = "REML")
summary(full.model3wnf)
(exp(full.model3wnf$b)-1)*100
(exp(full.model3wnf$se)-1)*100

full.model3mnf <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0imnf,
                         method = "REML")
summary(full.model3mnf)

## N concentration
full.model3lnc <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0ilnc,
                         method = "REML")
summary(full.model3lnc)
(exp(full.model3lnc$b)-1)*100
(exp(full.model3lnc$se)-1)*100

full.model3wnc <- rma.mv(yi, 
                         vi, 
                         #random = ~ 1 | Site_ID, #individual effect size nested within basin (Basin is level 3, effect size is level 2)
                         tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                         data = data_es0iwnc,
                         method = "REML")
summary(full.model3wnc)
(exp(full.model3wnc$b)-1)*100
(exp(full.model3wnc$se)-1)*100

#### STEP 4 Meta-regressions Total Litterfall ####

#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ia$stormfreq=z.trans(data_es0ia$StormFrequencyNorm)
data_es0ia$timesincestorm=z.trans(data_es0ia$YearsSinceLastStorm)
data_es0ia$soilP=z.trans(data_es0ia$Other_soil_P)
data_es0ia$windsp=z.trans(data_es0ia$WMO_wind_ms)
data_es0ia$distrain=z.trans(data_es0ia$Disturb_Rainfall_mm)

##Case studies with Hurrecon wind data
hurr_sites<-data_es0ia %>% filter(hurrwind!="NA")
str(hurr_sites)#n = 46
hurr_sites$hurrwind=z.trans(hurr_sites$HURRECON_wind_ms)

hurr_sites$pred_site_wind<- 9.40127+(0.64772*hurr_sites$WMO_wind_ms)-(0.03751*hurr_sites$Distance_to_Disturb_km)
hurr_sites$pred_site_wind

hurr_sites_amb<-data_es0ia %>% filter(hurrwind!="NA")%>% filter(Treatment=="Ambient")
str(hurr_sites_amb)#n = 42
hurr_sites_amb$hurrwind=z.trans(hurr_sites_amb$HURRECON_wind_ms)

#transforming the predicted site wind speed
hurr_sites$predwind=z.trans(hurr_sites$pred_site_wind)

## Total Litterfall Response as a function of soil P####

#Table4 - Model1a
model.mods3full<-rma.mv(yi,vi, 
                        random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = data_es0ia,
                        method = "REML",
                        mods = ~soilP)#using z transformed soil P
summary(model.mods3full)
tab_model(model.mods3full,p.val="kr")
#Significant intercept and regression coefficient
#The intercept represents the mean effect of mean soil P

#under ambient conditions
model.mods3_amb<-rma.mv(yi,vi, 
                        random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = hurr_sites,
                        method = "REML",
                        mods = ~soilP)#using z transformed soil P
summary(model.mods3_amb)

#testing other random effects structure - 
model.mods3full.1<-rma.mv(yi,vi, 
                        random = ~ 1 |Region/DisturbanceName, 
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

#sensitivity analysis - testing the effect of soil P when 5 sites are removed

model.mods3full_red<-rma.mv(yi,vi, 
                        random = list(~ 1 | Site, ~1|DisturbanceName), 
                        tdist = TRUE,data = red_sites,
                        method = "REML",
                        mods = ~soilP)#using z transformed soil P
summary(model.mods3full_red)
#No change in significance of soil P when sites are deleted

model.mods3full_red2<-rma.mv(yi,vi, 
                            random = ~ 1 |Region/DisturbanceName, 
                            tdist = TRUE,data = red_sites,
                            method = "REML",
                            mods = ~soilP)#using z transformed soil P
summary(model.mods3full_red2)

#no_cte_totlitfall

model.mods3full_nocte<-rma.mv(yi,vi, 
                            random = list(~ 1 | Site, ~1|DisturbanceName), 
                            tdist = TRUE,data = no_cte_totlitfall,
                            method = "REML",
                            mods = ~soilP)#using z transformed soil P
summary(model.mods3full_nocte)#no difference

#Bringing disturbane intensity metrics along with soil P
model_a<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~hurrwind)#standardized HURRECON site-level wind speed
summary(model_a)#Hurrecon wind speed has a positive significant effect on litterfall response

model_b<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~predwind)#not great
summary(model_b)#not so great

#Table 2 - Model2a####
model_tab4<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*hurrwind)
summary(model_tab4)#Both positive predictors, no significant interaction

model_tab4S<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*hurrwind)
summary(model_tab4S)#strange

#changing random effects structure
model_b2<-rma.mv(yi,vi,random = ~ 1 | Site, tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*windsp*d2track)
summary(model_b2)

model_a3<-rma.mv(yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                 data = hurr_sites,method = "REML",
                 mods = ~soilP+hurrwind)
summary(model_a3)

hurr_sites$pred_site_wind<- 8.24618+0.67714*hurr_sites$WMO_wind_ms-0.04986*hurr_sites$Distance_to_Disturb_km
hurr_sites$pred_site_wind
xyplot(hurr_sites$pred_site_wind~hurr_sites$HURRECON_wind_ms)

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
red_sites_hurr_amb <- hurr_sites_amb %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")

#Table S3####
model_a2_red<-rma.mv(yi, 
                vi, 
                random = list(~ 1 | Site, ~1|DisturbanceName), 
                tdist = TRUE, 
                data = red_sites_hurr,# n = 41
                method = "REML",
                mods = ~soilP*hurrwind)
summary(model_a2_red)#41 case studies - when the 5 sites are removed, there is a significant positive interaction betwen soil P and wind speed

model_a2_red_amb<-rma.mv(yi, 
                     vi, 
                     random = list(~ 1 | Site, ~1|DisturbanceName), 
                     tdist = TRUE, 
                     data = red_sites_hurr_amb,# n = 41
                     method = "REML",
                     mods = ~soilP*hurrwind)
summary(model_a2_red_amb)


#Cook's distance for best model
hatvalues(model_c,type="diagonal")
x <- cooks.distance(model_c)
x
plot(x, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance")

#Testing correlations among moderator variables in best model (model_c)
names(data_es0ia)
chart.Correlation(data_es0ia[,c(85,86,87)])#No correlation among moderators

#Funnel plot - not so useful because does not resolve the non-independence that the random effects resolves
funnel(model_c,xlab = "Litterfall Mass Response",studlab = TRUE,cex=2,cex.lab=1.5,cex.axis=1.4)

#### Leaf fall Response Meta-regression model ####
#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ilf$soilP=z.trans(data_es0ilf$Other_soil_P)
data_es0ilf$hurrwind=z.trans(data_es0ilf$HURRECON_wind_ms)
data_es0ilf$stormfreq=z.trans(data_es0ilf$StormFrequencyNorm)
data_es0ilf$timesincestorm=z.trans(data_es0ilf$YearsSinceLastStorm)

hurr_siteslf<-data_es0ilf %>% filter(hurrwind!="NA")
str(hurr_siteslf)#n = 31

#Leaf fall Table 2####
mixed_lf_tab4<-rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), 
               tdist = TRUE,data = data_es0ilf,
               method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4)

#Leaf fall Table S8#### alt random effect
mixed_lf_tab4a<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                  tdist = TRUE,data = data_es0ilf,
                  method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4a)

mixed_lf_tab4b<-rma.mv(yi,vi,random = ~ 1 | Country/DisturbanceName, 
                 tdist = TRUE,data = data_es0ilf,
                 method = "ML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4b) 

#comparing leaf fall models based on BIC
BIC(mixed_lf_tab4,mixed_lf_tab4a,mixed_lf_tab4b)

#sensitivity analysis - removing sites where total soil P was estimated from available soil P
red_siteslf <- data_es0ilf %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")
str(red_siteslf)#26
levels(as.factor(red_siteslf$Case_study))

hurr_red_siteslf <- hurr_siteslf %>% filter(Site!="Grande-Terre")%>% filter(Site!="Kumuwela")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")#%>% filter(Site!="Grande-Terre")
str(red_siteslf)#26

mixed_lf_red<-rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), 
                 tdist = TRUE,data = red_siteslf,
                 method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_red) 

#Same for observations with HURRECON data

#Table 4 Best Option for Leaf fall
mixed_lf_tab4.0<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                        tdist = TRUE,data = data_es0ilf,
                        method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.0) 
mixed_lf_tab4.1<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                 tdist = TRUE,data = hurr_siteslf,
                 method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.1) 

#Adding Cyclone metrics

##Predictions
##Table 4 Model 2b - Final Leaf Meta Regression####
mixed_lf_tab4.2<-rma.mv(yi,vi,random = list(~ 1 |Site,~1|DisturbanceName), 
                      tdist = TRUE,data = hurr_siteslf,
                      method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.2) 

mixed_lf_tab4.2a<-rma.mv(yi,vi,random = ~ 1 | Region/DisturbanceName, 
                   tdist = TRUE,data = hurr_siteslf,
                   method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.2a)#this works

mixed_lf_tab4.3a<-rma.mv(yi,vi,random = ~ 1 | Country/DisturbanceName, 
                         tdist = TRUE,data = hurr_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.3a)#this works

mixed_lf_tab4.2b<-rma.mv(yi,vi,random = list(~ 1 | Region, ~1|DisturbanceName), 
                         tdist = TRUE,data = hurr_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.2b)

mixed_lf_tab4.3b<-rma.mv(yi,vi,random = list(~ 1 | Country, ~1|DisturbanceName), 
                         tdist = TRUE,data = hurr_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.3b)

mixed_lf_tab4.2c<-rma.mv(yi,vi,random = ~1|DisturbanceName, 
                         tdist = TRUE,data = hurr_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.2c)

mixed_lf_tab4.2d<-rma.mv(yi,vi,random = ~1|Site, 
                         tdist = TRUE,data = hurr_siteslf,
                         method = "REML",mods = ~soilP*hurrwind)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf_tab4.2d)

#with reduced sites
mixed_lf.0_red<-rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), 
                   tdist = TRUE,data = red_siteslf,
                   method = "REML",mods = ~soilP+windsp+d2track)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf.0_red) 

mixed_lf.1<-rma.mv(yi,vi,random = list(~ 1 | Site,~1|DisturbanceName), 
                 tdist = TRUE,data = data_es0ilf,
                 method = "REML",mods = ~soilP+stormfreq)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf.1) 

mixed_lf.2<-rma.mv(yi,vi,random = ~ 1 | DisturbanceName, 
                 tdist = TRUE,data = data_es0ilf,
                 method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf.2)

mixed_lf.3<-rma.mv(yi,vi,random = list(~1|Site,~ 1 | DisturbanceName), 
                   tdist = TRUE,data = data_es0ilf,
                   method = "REML",mods = ~soilP)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf.3)

# compare models using BIC
BIC(mixed_lf,mixed_lf.1,mixed_lf.2,mixed_lf.3) 


mixed_lf2<-rma.mv(yi,vi,random = ~ 1 | Site, 
                 tdist = TRUE,data = data_es0ilf,
                 method = "REML",mods = ~soilP*windsp+soilP*d2track)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf2) 



mixed_lf3<-rma.mv(yi,vi,random = ~ 1 | Site/DisturbanceName, 
                  tdist = TRUE,data = data_es0ilf,
                  method = "REML",mods = ~soilP+windsp+d2track+soilP*d2track+soilP*windsp)#+WMO_wind_kts+Distance_to_Disturb_km)
summary(mixed_lf3)

res_results_leaf<-multimodel.inference(TE = "yi", 
                                      seTE = "vi",
                                      data = data_es0ilf,
                                      predictors = c("soilP", "windsp", "d2track"),
                                      interaction = TRUE,method="REML",test = "knha",eval.criterion='BIC')

summary(res_results_leaf)


chart.Correlation(data_es0ilf[,c(48,50,81)])

(sum(full.model3_red2$sigma2) - sum(model.mods3fullb2$sigma2)) / sum(full.model3_red2$sigma2)#30% reduction in deviance

## BIAS ####

funnel(model.mods3fullb2,xlab = "Litterfall Mass Response",studlab = TRUE,cex=2,cex.lab=1.5,cex.axis=1.4)
model.mods3full_e
funnel(model.mods3full_e,xlab = "Litterfall Mass Response",studlab = dataes0_new2_b$Effectsize_ID,cex=2,cex.lab=1.5,cex.axis=1.4)

ggsave(filename = "FUNNEL.png",
       plot = FUNNEL, width = 14, height = 11, units = 'cm',
       scale = 2, dpi = 600)
funnel(model.mods4full, xlab="Response", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))
eggers.test(x = model.mods4full)
trimfill(model.mods4full)
#+
legend(2, 1.5, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
       fill=c("darkblue","blue","lightblue"))

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

###RESPONSE predition data ####

pred.bm.resp<-data.frame(predict.rma(model.mods3fullb2))
pred.bm.resp


names(data_es0ia$StormFrequencyNorm)


####RESPONSE PREDICTIONS linear model####

pred.bm.resp<-data.frame(predict.rma(model.mods3fullb2))
pred.bm.resp
summary(model.mods3fullb2)

pred.bm.resp_final<-data.frame(predict.rma(model.mods3full_e2))
pred.bm.resp_final

names(dataes0_new)
lm_Rp2<-glm(yi~logsoilP,data=dataes0_new2)
summary(lm_Rp2)
nagelkerke(lm_Rp2)
lm_Rp2$fitted.values
with(summary(lm_Rp2), 1 - deviance/null.deviance)

lm_Rp2.1<-glm(yi~logsoilP+WMO_wind_kts,data=dataes0_new2)
summary(lm_Rp2.1)
nagelkerke(lm_Rp2.1)
lm_Rp2$fitted.values
with(summary(lm_Rp2.1), 1 - deviance/null.deviance)
anova(lm_Rp2.1,null,test="Chisq")

lm_Rp2.2<-glm(yi~logsoilP+WMO_wind_kts+Distance_to_Disturb_km,data=dataes0_new2)
summary(lm_Rp2.2)
nagelkerke(lm_Rp2.2)

with(summary(lm_Rp2.2), 1 - deviance/null.deviance)
anova(lm_Rp2.2,null,test="Chisq")

lm_Rp2.3<-glm(yi~WMO_wind_kts+Distance_to_Disturb_km,data=dataes0_new2)
summary(lm_Rp2.3)
nagelkerke(lm_Rp2.3)

null=glm(yi~1,data=dataes0_new2)
anova(lm_Rp2.3,lm_Rp2.2,lm_Rp2,null,test="Chisq")

lme_Rp<-lme(yi~logsoilP,random=~1|Site, data=dataes0_new2)
summary(lme_Rp)
library(rcompanion)
nagelkerke(lme_Rp)

library(lmtest)
lrtest(lme_Rp)
sjPlot::tab_model(lme_Rp)


lme_Rp2<-lme(yi~logsoilP+WMO_wind_kts,random=~1|Site, data=dataes0_new2)
summary(lme_Rp2)
sjPlot::tab_model(lme_Rp2)


#### FINAL linear mixed model to predict the Response ####

lm_Rp2b.1<-lmer(yi~(1|Site)+soilP+windsp+soilP*d2track, data=data_es0ia)#+(1|Site)
summary(lm_Rp2b.1)
sjPlot::tab_model(lm_Rp2b.1,show.df = TRUE,digits.p = 2,p.threshold = 0.05)

lm_Rp2b.2<-lmer(yi~soilP+(1|Site), data=data_es0ia)
summary(lm_Rp2b.2)
r.squaredGLMM(lm_Rp2b.1)
predict(lm_Rp2b.1)
fitted(lm_Rp2b)

lm_Rp2b.1lf<-lmer(yi~(1|Site)+soilP+windsp+soilP*d2track, data=data_es0ilf)#+(1|Site)
summary(lm_Rp2b.1lf)
sjPlot::tab_model(lm_Rp2b.1lf,show.df = TRUE,digits.p = 2,p.threshold = 0.05)
sjPlot::tab_model(lm_Rp2b.1,lm_Rp2b.1lf,show.df = TRUE,digits.p = 2,p.threshold = 0.05)
fixef(lm_Rp2b.1lf)
anova(lm_Rp2b.1lf)

r.squaredGLMM(lm_Rp2b.1lf)
modelplot(lm_Rp2b.1lf)

lm_Rp2b.1lf_b<-lme(yi~logsoilP+logwindsp+logdist, random = ~ 1|Site, data=data_es0ilf)
summary(lm_Rp2b.1lf_b)
calc.relimp(lm_Rp2b.1lf_b, type = c("lmg"), rela = TRUE)


lm_Rp2c<-lmer(yi~logsoilP+(1|Site)+WMO_wind_kts+YearsSinceLastStorm, data=dataes0_new2)#+(1|Site)
summary(lm_Rp2c)
sjPlot::tab_model(lm_Rp2c)

lm_Rp2c<-lmer(yi~logsoilP+(1|Site), data=dataes0_new2)#+(1|Site)
summary(lm_Rp2c)

sjPlot::tab_model(lm_Rp2b)
sjPlot::tab_model(model.mods6b,p.val = "kr", show.df = TRUE)

gam_Rp2a<-gamm(yi~s(logsoilP)+s(WMO_wind_kts)+s(Distance_to_Disturb_km),random=list(Site=~1), data=dataes0_new2)
summary(gam_Rp2a$gam)

shapiro.test(residuals(lm_Rp2))
modelplot(lm_Rp2)
calc.relimp(lm_Rp2a, type = c("lmg"), rela = TRUE)
predictors2 <- c('Log soil P','Wind speed','Log soil P')
relimp2 <- c(10.9,83.9, 16.1)
relimp.data2<- data.frame(predictors2, relimp2)
relimp.data2 %>%
  ggplot(aes(x=reorder(predictors2,relimp2),y=relimp2)) +
  geom_col(fill="steelblue") + 
  coord_flip()+theme_minimal()+theme(axis.title=element_text(size=20),
                                     axis.text=element_text(size=18))+ labs(x="Predictors", y="Relative importance (%)")+
  geom_text(aes(label=relimp2), position = position_stack(vjust= 0.5),
            colour = "white", size = 10)

lm_Rp3<-lm(yi~WMO_wind_kts+logsoilP,data=data_es0ia)#best model
summary(lm_Rp3)
shapiro.test(residuals(lm_Rp3))
modelplot(lm_Rp3)

mm_Rp3<-lme(yi~WMO_wind_kts+logsoilP+Distance_to_Disturb_km,random=~1|Site_ID, data=dataes0_new)
summary(mm_Rp3)
sjPlot::tab_model(mm_Rp3)
sjPlot::tab_model(model.mods6b,p.val = "kr", show.df = TRUE)
modelplot(mm_Rp3)
calc.relimp(mm_Rp3, type = c("lmg"), rela = TRUE)
predictors2 <- c('Wind speed','Log soil P')
tab_model(mm_Rp3)
sem.model.fits(mm_Rp3)

#calculating relative importance for each predictor
calc.relimp(lm_Rp3, type = c("lmg"), rela = TRUE)

predictors <- c('Wind speed','Log soil P','Storm frequency','Log soil P * Storm frequency')
relimp <- c(53.4, 11, 26.5,9)
relimp.data<- data.frame(predictors, relimp)
relimp.data %>%
  ggplot(aes(x=reorder(predictors,relimp),y=relimp)) +
  geom_col(fill="steelblue") + 
  coord_flip()+theme_minimal()+theme(axis.title=element_text(size=20),
                                     axis.text=element_text(size=18))+ labs(x="Predictors", y="Relative importance (%)")+
  geom_text(aes(label=relimp), position = position_stack(vjust= 0.5),
            colour = "white", size = 10)

lm_Rp3a<-lm(yi~WMO_wind_kts+logsoilP+StormsPerYearSince1955,data=data_es0ia)
summary(lm_Rp3a)
shapiro.test(residuals(lm_Rp3a))

datared<- data_es0ia %>% filter(Site!="Halemanu") %>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")%>% filter(Site!="Kumuwela")

lm_Rp4<-lm(yi~WMO_wind_kts+logsoilP*StormsPerYearSince1955,data=datared)
summary(lm_Rp4)
shapiro.test(residuals(lm_Rp4))

AIC(lm_Rp2,lm_Rp3,lm_Rp4)

####RESPONSE PREDICTION Figures####

#example
preds_x<-predict(model_c,levels=0, addx=TRUE) # let's transfer LRR to RR
preds_x

pred.bm.resp<-data.frame(predict.rma(model.mods4full))
pred.bm.resp

pred.bm.resp_final<-data.frame(predict.rma(model.mods3full_e2))
pred.bm.resp_final

metaregplot<- cbind(data.frame(data_es0ia$Site, data_es0ia$yi, fitted(lm_Rp2b.1),data_es0ia$soilP,data_es0ia$windsp,data_es0ia$d2track))
metaregplot

preg2<-ggplot(metaregplot, aes(x=fitted.lm_Rp2b.1., y=data_es0ia.yi))+geom_point(size=6,alpha=0.9,col="black",fill="#2F2747")+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=20))+ labs(x="Predicted response", y="Observed response")
rsq_label <- paste('R^2 == 0.73')
Fig5.1<-preg2+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "right",legend.justification = c("right"))+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,size="Ln wind speed",color="Ln soil P")+ annotate("text", x = 0.2, y = 6, label = rsq_label, size=10,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Old Fig5a
Fig5.1 #labs(y = bquote('Density Litterfall mass'~(g/m^2/day))
ggsave(filename = "Fig5a_new.png",
       plot = Fig5.1, width = 8, height = 8, units = 'cm',
       scale = 2, dpi = 600)

model.mods3full_ext1_zt.1$pred
mypreds_resp_tot<-predict(model.mods3full_ext1_zt.1,addx=TRUE)
data_es0ia$pred<-mypreds_resp_tot$pred
data_es0ia$se<-mypreds_resp_tot$se 

data_es0ia$d2track

##Predictions based on best Meta-Regression model for total litterfall####
#The best model is model_tab4
preds_x2<-predict(model_tab4,levels=0, addx=TRUE) # let's transfer LRR to RR
preds_x2<-data.frame(preds_x2)
str(preds_x2)

#Data frame
metaregplot.1_x2<- cbind(data.frame(hurr_sites$Site, hurr_sites$yi, hurr_sites$vi,preds_x2$pred,preds_x2$se,hurr_sites$soilP,hurr_sites$hurrwind,fac_soilP=factor(hurr_sites$Other_soil_P),hurr_sites$d2track,hurr_sites$Other_soil_P,hurr_sites$HURRECON_wind_ms,hurr_sites$Distance_to_Disturb_km))
str(metaregplot.1_x2)
levels(as.factor(hurr_sites$soilP))
metaregplot.1_x2$hurr_sites.Other_soil_P

#Figure
preg_tot<-ggplot(metaregplot.1_x2, aes(x=hurr_sites.HURRECON_wind_ms, y=preds_x2.pred))+geom_point(shape=21,aes(col=fac_soilP,size=preds_x2.se),stroke=1.4)
preg_tot<-preg_tot+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_tot<-preg_tot+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#878787","#313695","#8073ac","#4575b4","#74add1","#abd9e9","#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027","#a50026"))
preg_tot<-preg_tot+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
#geom_smooth(method = 'glm', col="#9f8a89",alpha=0.1,se=FALSE)
#preg_tot
#preg_tot<-preg_tot#+labs(x="Predicted response",y="Observed response")
#rsq_label <- paste('R^2 == 0.73')
preg_tot<-preg_tot+guides(size = guide_legend(override.aes = list(col = "black",shape=21)),color = guide_legend(override.aes = list(size = 8)))
#preg_tot
Fig5a_new<-preg_tot+theme_pubr()+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+
  theme(legend.position = "right", legend.justification = "right",legend.title = element_text (size = 21), legend.text = element_text (size = 20),axis.title=element_text(size=28),axis.text=element_text(size=26))+
  labs(size="",color="Total soil P (mg/kg)")+guides(size=FALSE)+ annotate("text", x = 18, y = 6.5, label = "a Total litterfall", size=10,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5a_new####
Fig5a_new

ggsave(filename = "Fig5a_Final.png",
       plot = Fig5a_new, width = 12, height = 12, units = 'cm',
       scale = 2, dpi = 600)
metaregplot.1_x2$hurr_sites.Other_soil_P
#Changing color scheme to include unidirectional and log soil P
preg_tot2<-ggplot(metaregplot.1_x2, aes(x=hurr_sites.HURRECON_wind_ms, y=preds_x2.pred))+geom_point(shape=21,aes(col=log(hurr_sites.Other_soil_P),size=preds_x2.se),stroke=1.6)
preg_tot2<-preg_tot2+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_tot2<-preg_tot2+scale_color_gradient(low="#FFFD7D",high="#D8001F")#scale_color_continuous(type="viridis") #FCFF00
preg_tot2
preg_tot2<-preg_tot2+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
preg_tot2<-preg_tot2+guides(size = guide_legend(override.aes = list(col = "black",shape=21)))#,color = guide_legend(override.aes = list(size = 8)))
preg_tot2
Fig5a_new2<- preg_tot2+theme_pubr()+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(legend.position = "top", legend.justification = "center", legend.title = element_text (size = 19), legend.text = element_text (size = 18),axis.title=element_text(size=26),
        axis.text=element_text(size=24))+
  labs(size="",color="Soil P \n(ln mg/kg)")+guides(size=FALSE)+
  annotate("text", x = 5, y = 6.6, label = "a Total litterfall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#Fig5a_new2####
Fig5a_new2

##Leaf fall response to cyclone####

#Predictions
preds_x_leaf<-predict(mixed_lf_tab4.2,levels=0, addx=TRUE) # let's transfer LRR to RR
preds_x_leaf<-data.frame(preds_x_leaf)
str(preds_x_leaf)

#Data frame
metaregplot.1_x_leaf<- cbind(data.frame(hurr_siteslf$Site, hurr_siteslf$yi, preds_x_leaf$pred,preds_x_leaf$se,hurr_siteslf$soilP,fac_soilP=factor(hurr_siteslf$Other_soil_P),hurr_siteslf$hurrwind,hurr_siteslf$Other_soil_P,hurr_siteslf$HURRECON_wind_ms))
names(metaregplot.1_x_leaf)
hurr_siteslf$fac_P=as.factor(hurr_siteslf$Other_soil_P)
levels(hurr_siteslf$fac_P)
hurr_siteslf$fac_P

preg_tot<-preg_tot+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#878787","#313695","#8073ac","#4575b4","#74add1","#abd9e9","#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027","#a50026"))

#Figure
preg_lf<-ggplot(metaregplot.1_x_leaf, aes(x=hurr_siteslf.HURRECON_wind_ms, y=preds_x_leaf.pred))+geom_point(shape=21,aes(col=fac_soilP,size=preds_x_leaf.se),stroke=1.4)
preg_lf<-preg_lf+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_lf<-preg_lf+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#8073ac","#4575b4","#abd9e9","#00876c","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027"))
preg_lf
preg_lf<-preg_lf+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
#geom_smooth(method = 'glm', col="#9f8a89",alpha=0.1,se=FALSE)
#preg_tot
#preg_tot<-preg_tot#+labs(x="Predicted response",y="Observed response")
#rsq_label <- paste('R^2 == 0.73')
preg_lf<-preg_lf+guides(size = guide_legend(override.aes = list(col = "black",shape=21)),color = guide_legend(override.aes = list(size = 8)))
preg_lf
Fig5b_new<-preg_lf+theme_pubr()+
  theme(legend.position = "top", legend.justification = "center",legend.title = element_text (size = 21), legend.text = element_text (size = 20),axis.title=element_text(size=28),axis.text=element_text(size=26))+
  labs(size="",color="Total soil P \n(mg/kg)")+guides(size=FALSE,color=FALSE)+ annotate("text", x = 16, y = 4.6, label = "b Leaf fall", size=10,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5b_new####
Fig5b_new

Fig5_Response<-Fig5a_new+Fig5b_new+plot_layout(ncol=2)+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig5_Response

ggsave(filename = "Fig5ab_Final2.png",
       plot = Fig5_Response, width = 20, height = 12, units = 'cm',
       scale = 2, dpi = 600)

#Changing color scheme to include unidirectional and log soil P
preg_lf2<-ggplot(metaregplot.1_x_leaf, aes(x=hurr_siteslf.HURRECON_wind_ms, y=preds_x_leaf.pred))+geom_point(shape=21,aes(col=log(hurr_siteslf.Other_soil_P),size=preds_x_leaf.se),stroke=1.6)
preg_lf2<-preg_lf2+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_lf2<-preg_lf2+scale_color_gradient(low="#FFFD7D",high="#D8001F")#scale_color_continuous(type="viridis")
preg_lf2
preg_lf2<-preg_lf2+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
preg_lf2<-preg_lf2+guides(size = guide_legend(override.aes = list(col = "black",shape=21)))#,color = guide_legend(override.aes = list(size = 8)))
Fig5b_new2<-preg_lf2+theme_pubr()+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(legend.position = "top",legend.justification = "center",legend.title = element_text (size = 19), legend.text = element_text (size = 20),axis.title=element_text(size=26),axis.text=element_text(size=24))+
  labs(size="",color="Soil P \n(ln mg/kg)")+guides(size=FALSE)+ annotate("text", x = 5, y = 5.2, label = "b Leaf fall", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5b_new2####
Fig5b_new2

#Final Fig5a-b with new color gradient
Fig5_Response2<-Fig5a_new2+Fig5b_new2+plot_layout(ncol=2)#+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig5_Response2

#Saving in High Res
ggsave(filename = "Fig5ab_Final3.png",
       plot = Fig5_Response2, width = 16, height = 10, units = 'cm',
       scale = 2, dpi = 1200)

#Supplements Pre Total Literfall vs Soil P####
data_resp_tot<-data_es0ia %>% filter(Treatment=="Ambient") %>% filter(Site!="Gadgarra")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")
str(data_resp_tot) #42 Obs Ambient only
unique(levels(data_resp_tot$Site))#26 sites

pre_tot_gam<-gam(Pre_Mean~s(log(Other_soil_P),k=6), data=data_resp_tot)
summary(pre_tot_gam)
tab_model(pre_tot_gam)
plot.gam(pre_tot_gam)

pre_tot<-lm(Pre_Mean~log(Other_soil_P), data=data_resp_tot)
summary(pre_tot)
tab_model(pre_tot)

pre_tot_lmer<-lmer(Pre_Mean~log(Other_soil_P)+(1|Country), data=data_resp_tot)
summary(pre_tot_lmer)
tab_model(pre_tot_lmer)

Fig5.2<-ggplot(data_resp_tot, aes(x=log(Other_soil_P), y=Pre_Mean))+geom_point(aes(color=log(Other_soil_P)),size=6,alpha=0.9,shape=21,stroke=2)+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=20))+ labs(y = bquote('Total litterfall'~(g/m^2/day)),x = "")
rsq_label.a <- paste('R^2 == 0.16')
Fig5.2<-Fig5.2+scale_color_gradient(low="#FCFF00", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,color="Soil P")+ annotate("text", x = 3, y = 6.5, label = rsq_label.a, size=6,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3a<-Fig5.2+ annotate("text", x = 3, y = 6, label = "p=0.009 N=42", size=6,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS3a<-FigS3a+ annotate("text", x = 3, y = 5.6, label = "y = -1.216 + 0.686*log soil P", size=5,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3a<-FigS3a+ annotate("text", x = 3, y = 7, label = "a", size=8,hjust=0,fontface="bold",colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3a

##Leaf fall Sup Figure Pre vs Soil P####
data_resp_leaf<-data_es0ilf %>% filter(Treatment=="Ambient")%>% filter(Site!="Gadgarra")%>% filter(Site!="Halemanu")%>% filter(Site!="Milolii")%>% filter(Site!="Makaha 1")
str(data_resp_leaf) #27 Obs Ambient only
unique(levels(data_resp_leaf$Site))#26 sites

prelf<-lm(Pre_Mean~log(Other_soil_P), data=data_resp_leaf)
summary(prelf)

#Leaf fall
Fig5.2lf<-ggplot(data_resp_leaf, aes(x=log(Other_soil_P), y=Pre_Mean))+geom_point(aes(color=log(Other_soil_P)),size=6,alpha=0.9,shape=21,stroke=2)+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=20))+ labs(y = bquote('Leaf fall'~(g/m^2/day)),x = "")
rsq_label.b <- paste('R^2 == 0.15')
Fig5.2lf<-Fig5.2lf+scale_color_gradient(low="#FCFF00", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,color="Soil P")+ annotate("text", x = 3, y = 4, label = rsq_label.b, size=6,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3b<-Fig5.2lf+ annotate("text", x = 3, y = 3.5, label = "p=0.04 N=27", size=6,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS3b<-FigS3b+ annotate("text", x = 3, y = 5.6, label = "y = -1.411 + 0.571*log soil P", size=5,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3b<-FigS3b+ annotate("text", x = 3, y = 4.5, label = "b", size=8,hjust=0,fontface="bold",colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS3b

FinalFigS3<-FigS3a+FigS3b+ plot_layout(ncol=2)
FinalFigS3

#FigS3####
ggsave(filename = "Final_Fig_S3.png",
       plot = FinalFigS3, width = 8, height = 15, units = 'cm',
       scale = 2, dpi = 600)

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
data_es0ia$Case_study<-c("Cubuy| Georges","Bisley| Georges","East Peak| Georges","East Peak| Georges","East Peak| Georges","East Peak| Georges","Bisley| Hugo","El Verde| Hugo",            
"East Peak| Hugo","Kokee| Iniki","Kokee P+| Iniki","Kokee NP+| Iniki","Kokee N+| Iniki","El Verde| CTE","Utuado| Georges","Guanica| Georges","Milolii| Iniki","Makaha 1| Iniki","Kumuwela| Iniki","Halemanu| Iniki",           
"Wooroonooran Basalt| Larry", "Wooroonooran Schist| Larry", "Mt Spec| Charlie","Mt Spec Disturbed| Charlie","Birthday Creek| Charlie","Birthday Creek| Aivu","Birthday Creek| Ivor","Bisley| Irma",              
"Guanica| Irma","Guayama| Irma","Rio Abajo| Irma","Bisley| Maria","Guanica| Maria","Guayama| Maria","Rio Abajo| Maria","Chamela-Cuixmala| Jova","Chamela-Cuixmala| Patricia","Grande-Terre| Hugo","Lienhuachi| Kalmaegi","Lienhuachi| Fungwong",     
 "Lienhuachi| Sinlaku", "Lienhuachi| Jangmi", "Kengting III| Mindulle", "Kengting IV| Mindulle", "Kengting III| Haima", "Kengting IV| Haima", "Kengting III| Nanmadol", "Kengting IV| Nanmadol", "Gadgarra | Keith")

forrest_data_C<-rbind(data.frame(ES=full.model3$yi,SE=sqrt(full.model3$vi),Type="Case",Site=data_es0ia$Site, Case_study=data_es0ia$Case_study, Cyclone=data_es0ia$DisturbanceName, SoilP=data_es0ia$Other_soil_P, Region=data_es0ia$Country))
forrest_data_C
forrest_data_C$SoilP=as.numeric(forrest_data_C$SoilP)
forrest_data_C$Case_study=as.factor(forrest_data_C$Case_study)
forrest_data_C$Case2<-factor(forrest_data_C$Case_study, levels=rev(levels(forrest_data_C$Case_study)))
str(forrest_data_C)#aes(x=reorder(Site,-yi), y=yi))

#calculating overall mean and high and low limits of confidence interval
ESm=full.model3$b
ESm
SEm=sqrt(full.model3$se)
SEm
ESm+(1.96*SEm)
ESm-(1.96*SEm)

## Colors used for the Regions, respectively: Australia, Guadeloupe, Hawaii, Mexico, Puerto Rico, Taiwan

###+labs(y = bquote('Density Litterfall mass'~(g/m^2/day)),x = "Months since disturbance")+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "b", face="bold",size=8,colour="black")

#setting pantropical response shading
full.model3$ci
b1 <- list(geom_hline(yintercept = 2.53, color = '#003f5c',alpha=0.4),
           annotate("rect", ymin = 1.16, ymax = 3.91,xmin = -Inf, xmax = Inf,alpha=0.1,linetype=2))
b1

## Fig3a #### Response of Total Litterfall Mass####
plot1D<-ggplot(data=forrest_data_C,aes(x=reorder(Case2,-ES),y=ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE),color=Region))+geom_pointrange(alpha=0.9,size=1.2)+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
plot2D<-plot1D+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1.2, alpha = 0.6)#+geom_hline(aes(yintercept=3.64), lty=2, colour = "#003f5c", cex=0.9, alpha = 0.4)
plot3D<-plot2D+xlab("Case study")+theme_bw()+ylab(expression(Response~(ln~litterfall[ti]/litterfall[t0]))) #ylab(expression(Response~(ln~Total~litterfall~t[i]~t[0]^-1)))
finalFig3a<-plot3D+theme(axis.title=element_text(size=26),axis.text.x=element_text(size=22,angle=0, hjust=0.5),axis.text.y=element_text(size=17,angle=0, hjust=1))+theme(legend.title = element_blank(), legend.text = element_text (size = 24),legend.position=c(0.8,0.9),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),                                                                                                                                                                        legend.box.background = element_rect(colour = "grey"))+b1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))
finalFig3a<-finalFig3a+ annotate("text", y = -1, x = 46, fontface="bold",label = "a", size=8,colour="black")
finalFig3a

#Saving in high res
ggsave(filename = "Fig3a_Response_New.png",
       plot = finalFig3a, width = 15, height = 18, units = 'cm',
       scale = 2, dpi = 600)
####Fig3b-d####

##Fig3b Data
data_frac4 <- rbind(data.frame(group="FFS", variable="Annual", estimate=full.model3ff$b, var="Mass flux",
                               ci_low=(full.model3ff$b-(1.96*full.model3ff$se)),ci_up=(full.model3ff$b+(1.96*full.model3ff$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Wood", variable="Annual",estimate=full.model3wf$b, var="Mass flux",ci_low=(full.model3wf$b-(1.96*full.model3wf$se)),ci_up=(full.model3wf$b+(1.96*full.model3wf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Leaf", variable="Annual",estimate=full.model3lf$b, var="Mass flux",ci_low=(full.model3lf$b-(1.96*full.model3lf$se)),ci_up=(full.model3lf$b+(1.96*full.model3lf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Misc.", variable="Annual",estimate=full.model3mf$b,var="Mass flux",ci_low=(full.model3mf$b-(1.96*full.model3mf$se)),ci_up=(full.model3mf$b+(1.96*full.model3mf$se)),
                               row.names=FALSE, stringsAsFactors=TRUE))
data_frac4

## Fig3c Data
data_impact_frac2 <- rbind(data.frame(group="Total", var="P flux",estimate2=full.model3tpf$b, ci_low2=(full.model3tpf$b-(1.96*full.model3tpf$se)),ci_up2=(full.model3tpf$b+(1.96*full.model3tpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Wood", var="P flux", estimate2=full.model3wpf$b,
                                      ci_low2=(full.model3wpf$b-(1.96*full.model3wpf$se)), ci_up2=(full.model3wpf$b+(1.96*full.model3wpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Leaf", var="P flux", estimate2=full.model3lpf$b,ci_low2=(full.model3lpf$b-(1.96*full.model3lpf$se)), ci_up2=(full.model3lpf$b+(1.96*full.model3lpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Misc.", var="P flux", estimate2=full.model3mpf$b, ci_low2=(full.model3mpf$b-(1.96*full.model3mpf$se)), ci_up2=(full.model3mpf$b+(1.96*full.model3mpf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Total", var="N flux",estimate2=full.model3tnf$b,ci_low2=(full.model3tnf$b-(1.96*full.model3tnf$se)), ci_up2=(full.model3tnf$b+(1.96*full.model3tnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Wood", var="N flux", estimate2=full.model3wnf$b, ci_low2=(full.model3wnf$b-(1.96*full.model3wnf$se)), ci_up2=(full.model3wnf$b+(1.96*full.model3wnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Leaf", var="N flux", estimate2=full.model3lnf$b,ci_low2=(full.model3lnf$b-(1.96*full.model3lnf$se)), ci_up2=(full.model3lnf$b+(1.96*full.model3lnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE),
                           data.frame(group="Misc.", var="N flux", estimate2=full.model3mnf$b,ci_low2=(full.model3mnf$b-(1.96*full.model3mnf$se)), ci_up2=(full.model3mnf$b+(1.96*full.model3mnf$se)),
                                      row.names=FALSE, stringsAsFactors=TRUE))
data_impact_frac2

##Fig3d Data
data_impact_fracPN <- rbind(data.frame(group="Wood", var="P conc.", estimate2=full.model3wpc$b,ci_low2=(full.model3wpc$b-(1.96*full.model3wpc$se)),ci_up2=(full.model3wpc$b+(1.96*full.model3wpc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Leaf", var="P conc.", estimate2=full.model3lpc$b,ci_low2=(full.model3lpc$b-(1.96*full.model3lpc$se)),ci_up2=(full.model3lpc$b+(1.96*full.model3lpc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Wood", var="N conc.", estimate2=full.model3wnc$b,
                                       ci_low2=(full.model3wnc$b-(1.96*full.model3wnc$se)),ci_up2=(full.model3wnc$b+(1.96*full.model3wnc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE),
                            data.frame(group="Leaf", var="N conc.", estimate2=full.model3lnc$b,
                                       ci_low2=(full.model3lnc$b-(1.96*full.model3lnc$se)),ci_up2=(full.model3lnc$b+(1.96*full.model3lnc$se)),
                                       row.names=FALSE, stringsAsFactors=TRUE))

data_impact_fracPN
full.model3lnc$b
((exp(full.model3lpc$se))-1)*100

##Fig3b
pfrac4<-ggplot(data_frac4, aes(x=group,y=estimate,ymax=ci_up,ymin=ci_low,shape=var))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6))+
  geom_pointrange(size=1.2, alpha=0.8,position=position_dodge(width=c(0.7, 1.2)))+#coord_flip()+
  geom_hline(aes(yintercept=0), lty=2,size=1,col="magenta",alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="", x="") +#scale_shape_discrete(solid=F)+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),
        axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),
        axis.title=element_text(size=20),
        axis.text=element_text(size=24),legend.text =  element_text(size=24),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.box="horizontal",legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.position = c(0.84,0.96),legend.box.background = element_blank())
pfrac4#+ scale_color_grey(start=0.65, end=0.25)
Fig3b<-pfrac4+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ annotate("text", y = 5.2, x = 0.5, fontface="bold",label = "b", size=8,colour="black")+ annotate("text", y = 0.4, x = 1.14, fontface="bold",label = "(11)", size=6,colour="black")+ annotate("text", y = 0.4, x = 2, fontface="bold",label = "(30)", size=6,colour="black")+ annotate("text", y = 0.4, x = 3, fontface="bold",label = "(31)", size=6,colour="black")+ annotate("text", y = 0.4, x =4.1, fontface="bold",label = "(9)", size=6,colour="black")
Fig3b

### Fig3c
data_impact_frac2
Fig3c<-ggplot(data_impact_frac2, aes(x=group,y=estimate2,ymax=ci_up2,ymin=ci_low2, shape = var,col=var))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6))+
  geom_pointrange(mapping=aes(shape=var),size=1.2, stroke=1.2,position=position_dodge(width=c(0.4, 0.8)))+
  geom_hline(aes(yintercept=0), lty=2,size=1.2,col="magenta", alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="Pantropical response", x="") +scale_shape_discrete(solid=F)+ scale_color_manual(values=c("#167923","#1C39A8"))+
  theme_bw() +# ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = 0.5),axis.text.x =element_text(vjust = 0.3),
        axis.title.y =element_text(vjust = 0.5,size=26),
        axis.text=element_text(size=24),legend.box="horizontal",legend.text =  element_text(size=24),legend.title = element_blank(),legend.position = c(.86,.88),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Fig3c<-Fig3c+ annotate("text", y = 6.2, x = 0.5, fontface="bold",label = "c", size=8,colour="black")+ annotate("text", y = 0.4, x = 0.8, fontface="bold",label = "(11)", size=6,colour="black")+ annotate("text", y = 0.4, x = 1.16, fontface="bold",label = "(10)", size=6,colour="black")+
  annotate("text", y = 0.4, x = 1.8, fontface="bold",label = "(9)", size=6,colour="black")+annotate("text", y = 0.4, x = 2.16, fontface="bold",label = "(8)", size=6,colour="black")+ annotate("text", y = 0.4, x = 2.8, fontface="bold",label = "(10)", size=6,colour="black")+annotate("text", y = 0.4, x = 3.16, fontface="bold",label = "(10)", size=6,colour="black")+ annotate("text", y = 0.4, x =3.8, fontface="bold",label = "(3)", size=6,colour="black")+ annotate("text", y = 0.4, x =4.16, fontface="bold",label = "(3)", size=6,colour="black")
Fig3c

##Fig3d
Fig3d<-ggplot(data_impact_fracPN, aes(x=group,y=estimate2,ymax=ci_up2,ymin=ci_low2, shape = var,col=var))+scale_y_continuous(breaks=c(-2,-1,0,1,2,4))+
  geom_pointrange(mapping=aes(shape=var),size=1.2, position=position_dodge(width=c(0.3, 0.6)))+#coord_flip()+
  scale_shape_discrete(solid=F)+
  geom_hline(aes(yintercept=0), lty=2,size=1.2,col="magenta", alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="", x="Litterfall fraction") +scale_color_manual(values=c("#167923","#1C39A8"))+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),
        axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),
        axis.title=element_text(size=26),
        axis.text=element_text(size=24),legend.box="horizontal",legend.text =  element_text(size=24),legend.title = element_blank(),legend.position = c(.86,.88),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Fig3d<-Fig3d+ annotate("text", y = 1.2, x = 0.45, fontface="bold",label = "d", size=8,colour="black")+ annotate("text", y = 0.06, x = 0.9, fontface="bold",label = "(3)", size=6,colour="black")+ annotate("text", y = 0.06, x = 1.09, fontface="bold",label = "(3)", size=6,colour="black")+
  annotate("text", y = 0.06, x = 1.85, fontface="bold",label = "(10)", size=6,colour="black")+annotate("text", y = 0.06, x = 2.09, fontface="bold",label = "(10)", size=6,colour="black")
Fig3d

#Fig3b-d
FinalFig3 <- Fig3b+Fig3c+Fig3d+plot_layout(ncol=1)
FinalFig3

#Fig3a-d
Fig3<-finalFig3a+FinalFig3+plot_layout(ncol=2)
Fig3

#Saving in high res the Final Fig3a-d
ggsave(filename = "Fig3ad_Final_v2.png",plot = Fig3, width = 22, height = 18, units = 'cm',scale = 2, dpi = 1000)

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

ESm=full.model3$b
ESm
SEm=sqrt(full.model3$se)
SEm
ESm+(1.96*SEm)
ESm-(1.96*SEm)
ci_up=(full.model3wpf$b+(1.96*full.model3wpf$se))
ci_up

####RESPONSE GGPLOT FIGURE N AND P CONCENTRATIONS####

names(nutmeta)

data_es0ilpc$yi

forrest_data_D<-rbind(data.frame(ES=full.model3lpc$yi,SE=sqrt(full.model3lpc$vi),Type="Case",Element="P",Fraction="Leaf",Site=data_es0ilpc$Site, Case_study=data_es0ilpc$Effectsize_ID, Cyclone=data_es0ilpc$Disturbance, Region=data_es0ilpc$Area),
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


####SupFigS5 Litterfall mass fractions with seasonality ####

full.model3ffS$k
data_frac3 <- rbind(data.frame(group="Total", variable="Annual (all)",estimate=full.model3$b, 
                               ci_low=full.model3$ci.lb, ci_up=full.model3$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Wood", variable="Annual (all)",estimate=full.model3wf$b, 
                               ci_low=full.model3wf$ci.lb, ci_up=full.model3wf$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Leaf", variable="Annual (all)",estimate=full.model3lf$b, 
                               ci_low=full.model3lf$ci.lb, ci_up=full.model3lf$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="FFS", variable="Annual (all)", estimate=full.model3ff$b, 
                               ci_low=full.model3ff$ci.lb, ci_up=full.model3ff$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Misc.", variable="Annual (all)",estimate=full.model3mf$b, 
                               ci_low=full.model3mf$ci.lb, ci_up=full.model3mf$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Total", variable="Subannual (red.)",estimate=full.model3S$b, 
                               ci_low=full.model3S$ci.lb, ci_up=full.model3S$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Wood", variable="Subannual (red.)",estimate=full.model3wfS$b,
                               ci_low=full.model3wfS$ci.lb, ci_up=full.model3wfS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Leaf", variable="Subannual (red.)",estimate=full.model3lfS$b, 
                               ci_low=full.model3lfS$ci.lb, ci_up=full.model3lfS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="FFS", variable="Subannual (red.)",estimate=full.model3ffS$b,
                               ci_low=full.model3ffS$ci.lb, ci_up=full.model3ffS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Misc.", variable="Subannual (red.)",estimate=full.model3mfS$b, 
                               ci_low=full.model3mfS$ci.lb, ci_up=full.model3mfS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Total", variable="Annual (red.)",estimate=full.model3SS$b,
                               ci_low=full.model3SS$ci.lb, ci_up=full.model3SS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Wood", variable="Annual (red.)",estimate=full.model3wfSS$b,
                               ci_low=full.model3wfSS$ci.lb, ci_up=full.model3wfSS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Leaf", variable="Annual (red.)",estimate=full.model3lfSS$b, 
                               ci_low=full.model3lfSS$ci.lb, ci_up=full.model3lfSS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="FFS", variable="Annual (red.)",estimate=full.model3ffSS$b,
                               ci_low=full.model3ffSS$ci.lb, ci_up=full.model3ffSS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(group="Misc.", variable="Annual (red.)",estimate=full.model3mfSS$b, 
                               ci_low=full.model3mfSS$ci.lb, ci_up=full.model3mfSS$ci.ub,
                               row.names=FALSE, stringsAsFactors=TRUE))
data_frac3

####FigS5 ggplot####

data_frac3

pfrac3<-ggplot(data_frac3, aes(x=group,y=estimate,ymax=ci_up,ymin=ci_low, shape = variable))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8))+
  geom_pointrange(mapping=aes(color=variable),size=1, position=position_dodge(width=c(0.8, 1.2)))+coord_flip()+
  geom_hline(aes(yintercept=0), lty=2,size=1,col="magenta",alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="Pantropical response with 95% CI", x="Litterfall fraction") +scale_shape_discrete(solid=F)+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),
        axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),
        axis.title=element_text(size=28),
        axis.text=element_text(size=26),legend.text =  element_text(size=20),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.position = c(0.84,0.73),legend.box.background = element_rect(colour = "black"))
FigS5<-pfrac3+ scale_color_manual(values=c("#141212","#5E3FBA","#A81C38"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#scale_color_grey(start=0.65, end=0.25)
FigS5

ggsave(filename = "FigS5_Fractions_Mass_Resp.png",
       plot = FigS5, width = 13, height = 10, units = 'cm',
       scale = 2, dpi = 600)

#### Supplementary Figures ####

## Litterfall fractions mass flux

forrest_data_lf<-rbind(data.frame(ES=full.model3lf$yi,SE=sqrt(full.model3lf$vi),Type="Case",Site=data_es0ilf$Site, Case_study=data_es0ilf$Effectsize_ID, Cyclone=data_es0ilf$DisturbanceName, SoilP=data_es0ia$Other_soil_P, Region=data_es0ia$Country))
forrest_data_lf
forrest_data_C$SoilP=as.numeric(forrest_data_C$SoilP)
forrest_data_C$Case_study=as.factor(forrest_data_C$Case_study)
forrest_data_C$Case2<-factor(forrest_data_C$Case_study, levels=rev(levels(forrest_data_C$Case_study)))
str(forrest_data_C)
levels(forrest_data_C$Site2)#aes(x=reorder(Site,-yi), y=yi))

plot1C<-ggplot(data=forrest_data_C,aes(x=reorder(Case2,-ES),y=ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE)))+geom_pointrange(aes(color=Region),alpha=0.8)+guides(title="Region")

plot2C<-plot1C+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = 0.6)+geom_hline(aes(yintercept=3.88), lty=2, color = "black", cex=0.9, alpha = 0.6)

plot3C<-plot2C+xlab("Case study")+ylab("Response [ln(litterfall ti/t0)]")+theme_bw()                                                                                      #color="black")

plot3C+theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17,angle=0, hjust=1))+theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.85,0.9),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+b1+scale_color_discrete("Dark2")

summary(full.model3)
#calculating overall mean and high and low limits of confidence interval
ESm=full.model3$b
ESm
SEm=sqrt(full.model3$se)
SEm
ESm+(1.96*SEm)
ESm-(1.96*SEm)

## Colors used for the Regions, respectively: Australia, Guadeloupe, Hawaii, Mexico, Puerto Rico, Taiwan

###+labs(y = bquote('Density Litterfall mass'~(g/m^2/day)),x = "Months since disturbance")+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "b", face="bold",size=8,colour="black")


b1 <- list(geom_hline(yintercept = 3.64, color = '#003f5c',alpha=0.2),
           annotate("rect", alpha = .1,ymin = 2.37, ymax = 4.91,xmin = -Inf, xmax = Inf,alpha=0.2))
b1

##Fig3a

plotlf<-ggplot(data=forrest_data_C,aes(x=reorder(Case2,-ES),y=ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE),color=Region))+geom_pointrange(alpha=0.9,size=1.2)+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

plotlf

plotlf2<-plotlf+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = 0.6)+geom_hline(aes(yintercept=3.64), lty=2, colour = "#003f5c", cex=0.9, alpha = 0.4)
plotlf2

plotlf3<-plotlf2+xlab("Case study")+ylab("Response [ln(litterfall ti/t0)]")+theme_bw()                                                                                    

finalFig3a<-plot3D+theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17,angle=0, hjust=1))+theme(legend.title = element_blank(), legend.text = element_text (size = 22),legend.position=c(0.85,0.9),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
                                                                                                                                                                         legend.box.background = element_rect(colour = "grey"))+b1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))
finalFig3a
ggsave(filename = "Fig3a-Response.png",
       plot = finalFig3a, width = 12, height = 15, units = 'cm',
       scale = 2, dpi = 600)
