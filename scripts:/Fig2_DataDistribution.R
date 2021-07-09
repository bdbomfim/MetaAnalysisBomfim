#01-07-2020
#020520
#032720
#033120
#042220
#051920 updated litterfall for sites with disturbances occurring within 2 months, but still need to re-do disturbance metrics
#060120 updated Bisley input info
#071320 updated removed Tanguro and added dummy variables to fit multilevel model
#072820 adding treatment effect and calculating total soil P for hawaiian sites

getRversion()
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggcorrplot")
install.packages('ggcorrplot')
library(ggcorrplot)
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
library(viridis)
library(viridisLite)
library(lsmeans)
update.packages("insight")
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
citation()
#"black"   "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "gray62" 
show_col(palette("R4"))

library(hrbrthemes)
ggsave(filename, plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE, ...)

install.packages("devtools")
devtools::install_github("daniel1noble/metaDigitise")
library(metaDigitise)

devtools::install_github("MathiasHarrer/dmetar")

####STEP 0 Upload data####

metadat<-read.csv(file.choose())#20210209_Litterfall_Mass_Flux
attach(metadat)
str(metadat)#2370 obs of 76 variables
names(metadat)
summary(metadat)
levels(metadat$Site)
2370-414

#Colors for the countries +scale_colour_manual(values=c("#33FF99","#0099CC","#FF0099","#009966","#C00000","#FFFF66"))


metadat$Total_soil_P..Bicarb_NaOH_Pt.=as.numeric(metadat$Total_soil_P..Bicarb_NaOH_Pt.)
metadat$Other_soil_P=as.numeric(metadat$Other_soil_P)
metadat$MAP_mm=as.numeric(metadat$MAP_mm)
metadat$Pre_SD_MonthSpecific=as.numeric(metadat$Pre_SD_MonthSpecific)
metadat$Pre_Mean_MonthSpecific=as.numeric(metadat$Pre_Mean_MonthSpecific)
metadat$Pre_Mean_MonthSpecific=as.numeric(metadat$Pre_Mean_MonthSpecific)

nutmeta<-read.csv(file.choose())#20210204_Nutrient_Meta_updated
attach(nutmeta)
str(nutmeta)#2551 of 78 variables
summary(nutmeta)

nutmeta$Total_Soil_P..Bicarb....NaOH.Pt.=as.numeric(nutmeta$Total_Soil_P..Bicarb....NaOH.Pt.)
nutmeta$MAP_mm=as.numeric(nutmeta$MAP_mm)
nutmeta$Elevation_m=as.numeric(nutmeta$Elevation_m)

######## STEP 1 #### Wrangling Data for Calculations #######

#Data for Resistance index - impact 0-0.5

rt <- metadat %>% filter(Cat_TSD_months == "0-0.5") #including data on cyclone-induced litterfall
rt

rec <- metadat %>% filter(Cat_TSD_months == "Rec")#including data for recovery and resilience calculations
rec

#This is for Resilience calculations for litterfall mass flux, excluding the duplicated study in Bisley

res <- rec %>% filter(Fraction=="TotLitfall") %>% filter(Case_ID!="25.2") #%>%filter(DisturbanceName!="Luis")%>%filter(DisturbanceName!="Ivor") #%>% mutate (Res = (Resilience)*100, ResS = (ResilienceS)*100, deltaRes = Res - ResS)
str(res)

#This is for Recovery calculations for litterfall mass flux

allres <- recmeta %>% filter(Fraction=="TotLitfall") %>% filter(Case_ID!="25.2")%>%filter(DisturbanceName!="Luis")%>%filter(DisturbanceName!="Ivor") #%>% mutate (Res = (Resilience)*100, ResS = (ResilienceS)*100, deltaRes = Res - ResS)
str(allres)

recmeta <- rec %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 37)%>% filter(Site!="San Felipe") %>%filter(DisturbanceName!="Luis") #%>%filter(DisturbanceName!="Ivor") %>% mutate (Res = (Resilience)*100, ResS = (ResilienceS)*100, deltaRes = Res - ResS)
recmeta

recov<-rec %>% filter(Fraction=="TotLitfall")%>% filter (TSD_months < 37) %>% filter(Cyclone_Mean!="NA") %>% filter(Case_ID!="25.2")



figRR1 <- ggline(metadat, x = "Forest_type",y = "Pre_Mean",linetype = "Fraction_Component", shape = "Fraction_Component",
               palette ="rainbow",add="mean",
               width=0.2,ylab="Stock (kg/ha)",xlab="Months since disturbance",short.panel.labs = FALSE,font.label = list(size = 14, face = "bold"),font.main=c(16,"plain"),font.x=c(16,"bold"))
ggpar(figRR1,font.x=c(10),font.y=c(18),legend="top",font.xtickslab = 10,font.ytickslab = 18,font.legend = c(16))



p_1<- ggboxplot(data1, x = "Forest_type", y = "Pre_Mean",
                    palette ="npg",color="Site",alpha=0.5,shape = "Forest_type",size=0.6,
                    add = "jitter", width=0.8,xlab="Forest",short.panel.labs = FALSE,font.label = list(size = 14, face = "bold"),font.main=c(16,"plain"),font.x=c(16,"bold"))
p_1
p_2<-facet(p_1 + theme_light(),facet.by="Fraction_Component",
                short.panel.labs = FALSE,linetype="dashed",size=3,
                panel.labs.background = list(fill = "white", color = "white"),panel.labs.font = list(color = "black", size = 22),panel.labs=list(Fraction_Component=c("Leaf fall","Misc fall","Wood Fall","TotLitfall")))

ggpar(p_2,font.x=c(20,"bold"),font.y=c(20,"bold"),legend="top",font.xtickslab = 18,font.ytickslab = 18,xlab="Treatment",ylab="Mass flux (g/m2/day)",font.legend = c(20,"bold"))


p_4<- ggscatter(metadat, aes(x = "Forest_type", y = "Pre_Mean"),
                palette ="npg",color="Treatment",alpha=0.5,
                xlab="Forest",short.panel.labs = FALSE,font.label = list(size = 14, face = "bold"),font.main=c(16,"plain"),font.x=c(16,"bold"))
p_4
p_4a<-facet(p_4 + theme_light(),facet.by="Fraction_Component",
           short.panel.labs = FALSE,linetype="dashed",size=3,
           panel.labs.background = list(fill = "white", color = "white"),panel.labs.font = list(color = "black", size = 22),panel.labs=list(Fraction_Component=c("Leaf fall","Misc fall","Wood Fall","TotLitfall")))

ggpar(p_4a,font.x=c(20,"bold"),font.y=c(20,"bold"),legend="top",font.xtickslab = 18,font.ytickslab = 18,xlab="Treatment",ylab="Mass flux (g/m2/day)",font.legend = c(20,"bold"))



figM <- ggscatter(metadat, aes(x = "Forest_type", y = "dif.dif"),alpha=0.5,size=3,conf.int = TRUE,
                  ylab="Total Litterfall Production (g/m2/day)",xlab="Forest type",short.panel.labs = FALSE,font.label = list(size = 14, face = "bold"),font.main=c(16,"plain"),font.x=c(16,"bold"))
figM
ggpar(figM,font.x=c(20,"bold"),font.y=c(20,"bold"),legend="top",font.xtickslab = 18,font.ytickslab = 18,xlab="Forest type",ylab="Total Litterfall Production (g/m2/day)",font.legend = c(20,"bold"))

figM <- ggline(metadat, x = "Forest_type",y = "Pre_Mean",linetype = "Fraction_Component", shape = "Fraction_Component",
               palette ="rainbow",add="mean",
               width=0.2,ylab="Stock (kg/ha)",xlab="Months since disturbance",short.panel.labs = FALSE,font.label = list(size = 14, face = "bold"),font.main=c(16,"plain"),font.x=c(16,"bold"))
ggpar(figM,font.x=c(10),font.y=c(18),legend="top",font.xtickslab = 10,font.ytickslab = 18,font.legend = c(16))


#Calculate effect size

names(metadat)

data_es <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                  sd1i = sd1i, sd2i = sd2i, data = data, measure = "ROM")

#######Data Viz Meeting Lara 03/30 and updated to 04/13#####

s1<-metadat %>%
  group_by(Region, Site, Holdridge_life_zone) %>%
  count() %>%
  ungroup() %>%
  mutate(prop = n / sum(n))
s1
str(s1)

names(metadat)
levels(metadat$Holdridge_life_zone)

df_1 <- metadat%>%
  group_by(Region, Holdridge_life_zone) %>%
  summarise(counts = n())
df_1

dfex <- df_1 %>%
  arrange(desc(Region)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfex, 4)

s2<-metadat %>%
  group_by(Region) %>%
  count()
s2
880+7+54+152
880/1093
p <- s2 %>%
  ggplot(aes(n)) +geom_histogram(alpha = 0.4)+
  theme_bw()
p

ggplot(metadat, aes(Region)) +
  geom_bar(fill = "#0073C2FF") +
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20))

#Litterfall Mass Flux Data distribution by Region

df <- metadat %>%
  group_by(Region, Holdridge_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df

plot_f<-ggplot(df, aes(x=reorder(Region, -counts),y=counts, fill = Holdridge_life_zone)) 
plot_f+geom_bar(stat="identity") +
  geom_text(aes(label = counts), vjust = -0.5, size=7, color="black")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20))

plot_f2<-ggplot(df, aes(x=reorder(Region, -counts),y=counts, fill = Holdridge_life_zone)) 
plot_f2+geom_bar(stat="identity") +
  geom_text(aes(label = sprintf("%0.2f", round(percent, digits = 1))), vjust = 0.3, size=7, color="black", position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20))

plot_f3<-ggplot(df, aes(x=reorder(Region, -counts),y=counts, fill = Holdridge_life_zone)) 
plot_f3+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+labs(y = "Number of observations", x = "Region")+theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

#Same for subset data0a including Ambient, Total litterfall, 0-0.5

df0a <- data0a %>%
  group_by(Region, Holdridge_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df0a

plot_0a<-ggplot(df0a, aes(x=reorder(Region, -counts),y=counts, fill = Holdridge_life_zone)) 
plot_0a+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+labs(x="Region", y="Number of observations")+
  theme(legend.title = element_text (size = 30), legend.text = element_text (size = 28), legend.position = c(.90,.95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

#Same for subset data0a including Ambient, Leaf fall, 0-0.5

df0b <- data0b %>%
  group_by(Region, Holdridge_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df0b

plot_0b<-ggplot(df0b, aes(x=Region,y=counts, fill = Holdridge_life_zone)) 
plot_0b+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

pr <- metadat %>% filter (Region == "Puerto Rico")

plotpr<-ggplot(pr, aes(x=Site,y=MAP_mm)) 
plotpr
plotpr+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

  group_by(Region, Holdridge_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df0b

plot_0b<-ggplot(df0b, aes(x=Region,y=counts, fill = Holdridge_life_zone)) 
plot_0b+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


#Litterfall Nutrient Flux Data distribution by Region

names(nutmeta)

df_nutflux <- nutmeta %>% filter(Raw_Unit == "mg/m2/day")%>%
  group_by(Region, Hold_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df_nutflux

plot_nutflux<-ggplot(df_nutflux, aes(x=Region,y=counts, fill = Hold_life_zone)) 
plot_nutflux+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

#Site Characteristics Multi-panel

names(metadat)
str(metadat)
metadat$Elevation_m=as.numeric(metadat$Elevation_m)
metadat$MAP_mm

metadatlong <- gather(metadat, key="measure", value="value", c("Elevation_m", "TempPrec_ratio_x10.2_Cmmyr", "Other_soil_P"))

pr <- metadat %>% filter(Region == "Puerto Rico") %>% group_by(Site,Holdridge_life_zone,Disturbance) %>% summarize (Temp_Prec = mean(TempPrec_ratio_x10.2_Cmmyr))
pr
pr %>% ggplot() + geom_col(mapping=aes(x=reorder(Site, -Temp_Prec), y=Temp_Prec, fill=Holdridge_life_zone))+theme_pubclean()+theme(axis.title=element_text(size=26),
                                                                                                                                   axis.text=element_text(size=26))+guides(fill = guide_legend(title = "Holdridge life zone"))+
  labs(x="Site", y="MAT/MAP x 100")+theme(legend.title = element_text (size = 30), axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text (size = 26), legend.position = c(.90, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


pr %>% ggplot() + geom_col(mapping=aes(x=reorder(Site, -Temp_Prec), y=Temp_Prec, fill=Disturbance))+theme_pubclean()+theme(axis.title=element_text(size=26),
                                                                                                                                   axis.text=element_text(size=26))+guides(fill = guide_legend(title = "Wind disturbance"))+
  labs(x="Site", y="MAT/MAP x 100")+theme(legend.title = element_text (size = 30), axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text (size = 26), legend.position = c(.90, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

aus <- metadat %>% filter(Region == "Australia") %>% group_by(Site,Holdridge_life_zone,Disturbance) %>% summarize (Temp_Prec = mean(TempPrec_ratio_x10.2_Cmmyr))
aus

aus %>% ggplot() + geom_col(mapping=aes(x=Site, y=Temp_Prec,fill=Disturbance))+theme_pubclean()+theme(axis.title=element_text(size=26),
                                                                                                                                   axis.text=element_text(size=26))+guides(fill = guide_legend(title = "Wind disturbance"))+
  labs(x="Site", y="MAT/MAP x 100")+theme(legend.title = element_text (size = 30), axis.text.x = element_text(angle=25, hjust=1),legend.text = element_text (size = 26), legend.position = c(.80, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


haw <- metadat %>% filter(Region == "Hawaii") %>% group_by(Site,Holdridge_life_zone,Disturbance) %>% summarize (Temp_Prec = mean(TempPrec_ratio_x10.2_Cmmyr))
haw

haw %>% ggplot() + geom_col(mapping=aes(x=reorder(Site,-Temp_Prec), y=Temp_Prec,fill=reorder(Holdridge_life_zone,-Temp_Prec)))+theme_pubclean()+theme(axis.title=element_text(size=26),
                                                                                                      axis.text=element_text(size=26))+guides(fill = guide_legend(title = "Holdridge life zone"))+
  labs(x="Site", y="MAT/MAP x 100")+theme(legend.title = element_text (size = 30), axis.text.x = element_text(angle=25, hjust=1),legend.text = element_text (size = 26), legend.position = c(.80, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


haw %>% ggplot() + geom_col(mapping=aes(x=reorder(Site, -Temp_Prec), y=Temp_Prec, fill=Disturbance))+theme_pubclean()+theme(axis.title=element_text(size=26),
                                                                                                                           axis.text=element_text(size=26))+guides(fill = guide_legend(title = "Wind disturbance"))+
  labs(x="Site", y="MAT/MAP x 100")+theme(legend.title = element_text (size = 30), axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text (size = 26), legend.position = c(.90, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))




names(metadatlong)
head(metadatlong)

variable_names <- list(
  "Elevation_m" = "Elevation (m)" ,
  "TempPrec_ratio_x10.2_Cmmyr" = "MAT/MAP x 100",
  "Other_soil_P" = "Soil P (mg/kg)"
)

variable_labeller <- function(variable,value){
  return(variable_names[value])
}

ggplot(pr, aes(x=Site, y=value, fill = Region))+
  geom_bar(stat='identity', fill="royalblue")+
  facet_wrap(~measure, ncol=1,strip.position = "left", scales = "free")+theme_bw()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Region"))+
  theme(legend.title = element_text (size = 26), axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


psite<-ggplot(metadat, aes(x=Site,y=Elevation_m, fill = Region)) 
psite+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Region"))+
  theme(legend.title = element_text (size = 26), axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


#Litterfall Nutrient Concentration Data distribution by Region

df_nutconc <- nutmeta %>% filter(Raw_Unit == "mg/g")%>%
  group_by(Region, Hold_life_zone) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df_nutconc

plot_nutconc<-ggplot(df_nutconc, aes(x=Region,y=counts, fill = Hold_life_zone)) 
plot_nutconc+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24))+guides(fill = guide_legend(title = "Holdridge Life Zone"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.50, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


#Litterfall Mass Flux Data distribution by Disturbance Type

dfa <- metadat %>%
  group_by(Dist_type) %>%
  summarise(counts = n())
dfa

df_a <- metadat %>%
  group_by(Region, Dist_type) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_a

plot_Dtype<-ggplot(df_a, aes(x=reorder(Dist_type, -counts),y=counts, fill = Region)) 
plot_Dtype+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=25, hjust=1))+guides(fill = guide_legend(title = "Region"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


pa<-ggplot(dfa, aes(x=reorder(Dist_type,-counts), y=counts)) 
pa+geom_bar(fill = "#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20))

#Litterfall Mass Flux Data distribution by Disturbance Name

dfb <- metadat %>%
  group_by(Disturbance) %>%
  summarise(counts = n())%>%
  arrange(desc(counts))
dfb

df_b <- metadat %>%
  group_by(Region, Disturbance) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_b

plot_Dname<-ggplot(df_b, aes(x=reorder(Disturbance, -counts),y=counts, fill = Region)) 
plot_Dname+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=40, hjust=1))+guides(fill = guide_legend(title = "Region"))+labs(x ="", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.75, .85),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

ggplot(dfb, aes(x=Disturbance, y=counts)) +
  geom_bar(fill = "#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20),axis.text.x = element_text(angle=45))

### trying to reorder the disturbance names

p<-ggplot(dfb, aes(x=reorder(Disturbance, -counts),y=counts))
p+geom_bar(fill = "#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=20),axis.text.x = element_text(angle=45))


#Litterfall Mass Flux Data distribution by Site

df_c <- metadat %>%
  group_by(Region, Site) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_c

plot_Site<-ggplot(df_c, aes(x=reorder(Site, -counts),y=counts, fill = Region)) 
plot_Site+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=45, hjust=1))+guides(fill = guide_legend(title = "Region"))+labs(y = "Number of observations", x = "")+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.75, .85),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


dfc <- metadat %>%
  group_by(Region, Site) %>%
  summarise(counts = n())%>%
  arrange(desc(counts))
dfc

pc<-ggplot(dfc, aes(x=reorder(Site, -counts),y=counts))
pc+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))

#Litterfall Mass Flux Data distribution by Year

dfc3 <- metadat %>%
  group_by(Year) %>%
  summarise(counts = n())%>%
  arrange(desc(counts))
dfc3

pc3<-ggplot(dfc3, aes(x=reorder(Year,-counts), y=counts))
pc3+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))

pc3a<-ggplot(dfc3, aes(x=Year,y=counts))
pc3a+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))

#Litterfall mass by Holdridge Life Zone
levels(metadat$Holdridge_life_zone)

df_1 <- metadat%>%
  group_by(Holdridge_life_zone) %>%
  summarise(counts = n())
df_1

dfex <- df_1 %>%
  arrange(desc(Holdridge_life_zone)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfex, 4)

pd<-ggplot(dfex, aes(x=reorder(Holdridge_life_zone,-counts), y=counts))
pd+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))

levels(metadat$Rock_type)
levels(metadat$Rock_type)=c("Acid volcanic", "Basalt", "BI volcanic", "Carbonate","Sandstone", "Shale", "Shield")

df_rock <- metadat %>%
  group_by(Rock_type) %>%
  dplyr::summarise(counts = dplyr::n())%>% dplyr::mutate(percent = (counts/sum(counts))*100)%>%arrange(desc(Rock_type))
df_rock

df_rock2 <- metadat %>%
  group_by(Rock_type, USDA_Soil_Order) %>%
  dplyr::summarise(counts = dplyr::n())%>% dplyr::mutate(percent = (counts/sum(counts))*100)%>%arrange(desc(Rock_type))
df_rock2

plot_rock<-ggplot(df_rock2, aes(y=counts, x = reorder(Rock_type,-counts), fill=USDA_Soil_Order)) 
plot_rock+geom_bar(stat="identity")+scale_fill_brewer(palette = "Set1")+theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "Soil order"))+labs(x = "Parent material", y = "Number of observations")+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 32), legend.position = c(.85, .96),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

names(metadat)
levels(metadat$RockPClass)=c("High", "Intermediate", "Low")
levels(metadat$RockPClass)
df_rock3 <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>%
  group_by(RockPClass,Rock_type, USDA_Soil_Order) %>%
  dplyr::summarise(counts = dplyr::n())%>% dplyr::mutate(percent = (counts/sum(counts))*100)%>%arrange(desc(Rock_type))
df_rock3

plot_rock3<-ggplot(df_rock3, aes(y=counts, x = RockPClass, fill=USDA_Soil_Order)) 
plot_rock3+geom_bar(stat="identity")+scale_fill_brewer(palette = "Set1")+theme_bw()+theme(axis.title=element_text(size=28),
                                                                                               axis.text=element_text(size=26),axis.text.x = element_text(angle=0, hjust=0.5))+guides(fill = guide_legend(title = "Soil order"))+labs(x = "Parent material P", y = "Number of observations")+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 32), legend.position = c(.95, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

plot_rock4<-ggplot(df_rock3, aes(y=counts, x = RockPClass, fill=Rock_type)) 
plot_rock4+geom_bar(stat="identity")+scale_fill_brewer(palette = "Set1")+theme_bw()+theme(axis.title=element_text(size=28),
                                                                                                axis.text=element_text(size=26),axis.text.x = element_text(angle=0, hjust=0.5))+guides(fill = guide_legend(title = "Parent material"))+labs(x = "Parent material P", y = "Number of observations")+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 32), legend.position = c(.95, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))



plot_ParMat+geom_bar(stat="identity")+scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "USDA soil order"))+labs(x = "Soil parent material", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.35, .96),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(0,200,400,600,800))


#Litterfall mass by TSD levels

library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

levels(metadat$Treatment)
str(metadat)

df_TSD <- metadat %>% filter(TSD_months < 61) %>%
  group_by(Country,Holdridge_life_zone, TSD_months)  %>%
  dplyr::summarise(counts = dplyr::n())
df_TSD

#Heatmap observations litterfall mass

df_TSD %>% ggplot(aes(x = Holdridge_life_zone, y = Region)) +geom_tile(aes(fill = counts))+theme_pubr()+
  #facet_grid(~ PlantDiv,switch = "x", scales = "free_x", space = "free_x")+
  scale_fill_gradient(name = "Observations",
                      low = "blue",
                      high = "orange")+coord_flip()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.x = element_text(angle=20, hjust=1,vjust = 1,size=18),axis.text.y=element_text(size=16),axis.title=element_blank(),
        legend.position="right",legend.text =  element_text(size=16,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.box.background = element_rect(colour = "black")) # use theme() to change font sizes)+
  #ggtitle(label = "servations")
  #scale_y_discrete(limits = rev(levels(as.factor(sbg$Site))))

levels(metadat$Raw_Unit)
df_TSDb <- metadat %>% filter(TSD_months < 37) %>%
  group_by(Region,Holdridge_life_zone) %>%
  dplyr::summarise(counts = dplyr::n())
df_TSDb

#Heatmap observations litterfall mass

df_TSDb %>% ggplot(aes(x = Holdridge_life_zone, y = Region)) +geom_tile(aes(fill = counts))+theme_pubr()+
  #facet_grid(~ PlantDiv,switch = "x", scales = "free_x", space = "free_x")+
  scale_fill_gradient(name = "Observations",
                      low = "blue",
                      high = "orange")+coord_flip()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.x = element_text(angle=20, hjust=1,vjust = 1,size=18),axis.text.y=element_text(size=16),axis.title=element_blank(),
        legend.position="right",legend.text =  element_text(size=16,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.box.background = element_rect(colour = "black")) # use theme() to change font sizes)+
#ggtitle(label = "servations")
#scale_y_discrete(limits = rev(levels(as.factor(sbg$Site))))

levels(nutmeta$Variable)

df_TSD_nut_flux <- nutmeta %>% filter(TSD_months < 61) %>% filter(Raw_Unit=="mg/m2/day") %>% filter (Variable=="N"|Variable=="P") %>% 
  group_by(Country,Holdridge_life_zone,TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())
df_TSD_nut_flux

df_TSD_nut_conc <- nutmeta %>% filter(TSD_months < 61) %>% filter(Raw_Unit=="mg/g") %>% filter (Variable=="N"|Variable=="P") %>% 
  group_by(Country,Holdridge_life_zone,TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())
df_TSD_nut_conc

#Heatmap observations litterfall nutrients

df_TSDnut %>% ggplot(aes(x = Holdridge_life_zone, y = Region)) +geom_tile(aes(fill = counts))+theme_pubr()+
  #facet_grid(~ PlantDiv,switch = "x", scales = "free_x", space = "free_x")+
  scale_fill_gradient(name = "Observations",
                      low = "blue",
                      high = "orange")+coord_flip()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.x = element_text(angle=20, hjust=1,vjust = 1,size=18),axis.text.y=element_text(size=16),axis.title=element_blank(),
        legend.position="right",legend.text =  element_text(size=16,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.box.background = element_rect(colour = "black")) # use theme() to change font sizes)+
#ggtitle(label = "servations")
#scale_y_discrete(limits = rev(levels(as.factor(sbg$Site))))

#Heatmap observations litterfall nutrients

df_TSDnut2 %>% ggplot(aes(x = Holdridge_life_zone, y = Region)) +geom_tile(aes(fill = counts))+theme_pubr()+
  #facet_grid(~ PlantDiv,switch = "x", scales = "free_x", space = "free_x")+
  scale_fill_gradient(name = "Observations",
                      low = "blue",
                      high = "orange")+coord_flip()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.x = element_text(angle=20, hjust=1,vjust = 1,size=18),axis.text.y=element_text(size=16),axis.title=element_blank(),
        legend.position="right",legend.text =  element_text(size=16,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.box.background = element_rect(colour = "black")) # use theme() to change font sizes)+
#ggtitle(label = "servations")
#scale_y_discrete(limits = rev(levels(as.factor(sbg$Site))))

####Figure 2 Study Site and Data Distribution####

mu <- ddply(df_TSD, "Country", summarise, grp.mean=counts(TSD_months))
head(mu)

mu2 <- ddply(df_TSD, "Country", summarise, grp.mean=counts(TSD_months))
head(mu2)

p1 <- ggplot(data=df_TSD, aes(x=TSD_months, fill=Basin)) +geom_density(adjust=1.5,alpha=.6) 
p1+theme_classic2()+theme(axis.title=element_text(size=32),axis.text.y=element_text(size=26),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=26))+guides(fill = guide_legend(title = "Basin"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 36), legend.text = element_text (size = 34), legend.position = c(.7, .85),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")+scale_color_brewer(palette="BrBG")

df_TSD <- metadat %>% filter(TSD_months < 61) %>%
  group_by(Country,Holdridge_life_zone,TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())
df_TSD

####Figure 2b-d####

#New colors +scale_fill_manual(values=c("#1dabe6","#1c366a","#c3ced0","#ffa600","#665191","#af060f"))

p1a <- ggplot(data=df_TSD, aes(x=TSD_months, fill=Country)) +geom_density(adjust=1.5,alpha=.8) 
p1a
pmass<-p1a+theme_classic2()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=19))+guides(color = guide_legend(title = "Country"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 12, 24, 36, 48, 60, 72, 84))+
  theme(legend.title = element_blank(), legend.text = element_text (size = 20), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = bquote('Density Litterfall mass'~(g/m^2/day)),x = "Months since disturbance")+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "b", face="bold",size=8,colour="black")
pmass

p_nut_flux <- ggplot(data=df_TSD_nut_flux, aes(x=TSD_months, fill=Country)) +geom_density(adjust=1.5,alpha=.8)
p_nut_flux
p_elem_flux<-p_nut_flux+theme_classic2()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=19))+guides(fill = guide_legend(title = "Country"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 12, 24, 36, 48, 60, 72, 84))+
  theme(legend.title = element_blank(), legend.text = element_text (size = 20), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y =  bquote('Density Litterfall N and P'~(mg/m^2/day)),x = "Months since disturbance")+scale_fill_manual(values=c("#1dabe6","#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "c", face="bold",size=8,colour="black")
p_elem_flux#+ annotate("text", x = 0, y = 0.15, label = "c", face="bold",size=8,colour="black")

df_TSD_nut_conc
p_nut_conc <- ggplot(data=df_TSD_nut_conc, aes(x=TSD_months, fill=Country)) +geom_density(adjust=1.5,alpha=.8) 
p_elem_conc<-p_nut_conc+theme_classic2()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=20))+guides(fill = guide_legend(title = "Country"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 12, 24, 36, 48, 60, 72, 84))+
  theme(legend.title = element_blank(), legend.text = element_text (size = 20), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y =  bquote('Density Litterfall N and P'~(mg/g)),x = "Months since disturbance")+scale_fill_manual(values=c("#c3ced0","#ffa600","#665191","#af060f"))+ annotate("text", x = 0, y = 0.15, label = "d", type="bold",size=8,colour="black")
p_elem_conc

library(patchwork)
Fig2b<-pmass+p_elem_flux+p_elem_conc
Fig2b
ggsave(filename = "Fig2b-d-final.png",
       plot = Fig2b, width = 19, height = 8, units = 'cm',
       scale = 2, dpi = 600)

names(rt)
df_TSD2 <- rt %>% filter(Fraction=="TotLitfall")%>%
  group_by(Basin,Region,Area,Site,Latitude,StormsPerYear,DisturbanceClass,Disturbanceme,WMO_wind_kts,Distance_to_Disturb_km,NYearsCoveredStorms,NStormsRegion) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(StormsPerYear))
df_TSD2

p2_b <- ggplot(data=df_TSD2, aes(x=StormsPerYear, group=Region, fill=Region)) +geom_density(adjust=1.5,alpha=0.5) 
p2_b+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=24),axis.text.x = element_text(angle=0, hjust=1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+scale_x_continuous(name="StormsPerYear")+
  theme(legend.title = element_text (size = 34), legend.text = element_text (size = 30), legend.position = c(.97, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Storms per year")

p2_c <- ggplot(data=df_TSD2, aes(x=StormsPerYear, group=Area, fill=Area)) +geom_bar(alpha=0.5) +scale_x_binned()
p2_c+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 30), legend.text = element_text (size = 28), legend.position = c(.97, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "Storms per year")

p2_d <- ggplot(data=df_TSD2, aes(y=StormsPerYear, x=reorder(Site,-StormsPerYear), color=Region)) +geom_point(stat="identity", size = 5)
p2_d+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=45, hjust=1,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 30), legend.text = element_text (size = 28), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(x = "Site",y = "Storms per year")

p2_e <- ggplot(data=df_TSD2, aes(x=WMO_wind_kts, group=Region, fill=Region)) +geom_bar(alpha=0.5) +scale_x_binned()
p2_e+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 30), legend.text = element_text (size = 28), legend.position = c(.3, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "WMO wind speed (kts)")

p2_f <- ggplot(data=df_TSD2, aes(x=WMO_wind_kts, group=Area, fill=Area)) +geom_bar(alpha=0.5) +scale_x_binned()
p2_f+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 24), legend.position = c(.3, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "WMO wind speed (kts)")

p2_g <- ggplot(data=df_TSD2, aes(x=Distance_to_Disturb_km, group=Region, fill=Region)) +geom_bar(alpha=0.5) +scale_x_binned()
p2_g+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "Distance to track (km)")

p2_h <- ggplot(data=df_TSD2, aes(x=Distance_to_Disturb_km, group=Area, fill=Area)) +geom_bar(alpha=0.5) +scale_x_binned()
p2_h+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "Distance to track (km)")

p2_i <- ggplot(data=df_TSD2, aes(y=StormsPerYear, group=Area, color=Area, x=Latitude)) +geom_point(alpha=0.7, size = 4)
p2_i+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.39, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Storms per year",x = "Latitude")

p2_j <- ggplot(data=df_TSD2, aes(x=Disturbanceme, group=Region, fill=Region)) +geom_bar(alpha=0.5)
p2_j+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=35, hjust=1,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.89, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "Cyclone name")

p2_k <- ggplot(data=df_TSD2, aes(x=Disturbanceme, group=Area, fill=Area)) +geom_bar(alpha=0.5)
p2_k+theme_grey()+theme(axis.title=element_text(size=22),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=35, hjust=1,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22), legend.position = c(.89, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Count",x = "Cyclone name")

legend_size <- c(1,3,6, 9,12)

cols <- c("#d73027",
          "darkgrey",
          "#4575b4")

df_TSD2$NYearsCoveredStorms

p2_l <- ggplot(data=df_TSD2, aes(y=NStormsRegion, x=reorder(Site,-NStormsRegion),size = NYearsCoveredStorms)) +geom_point(stat="identity",aes(size = NYearsCoveredStorms, colour = Region), alpha=0.7)
p2_l+theme_grey()+theme(axis.title=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=40, hjust=1,size=20))+
  theme(legend.title = element_text (size = 22), legend.text = element_text (size = 20), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right")+labs(x = "Site",y = "Number of storms")+ 
  theme(panel.background = element_rect(fill = "linen"))+scale_size(range = c(1,10), breaks =c(50,75,100,125,150))

names(recov)

df_TSD5 <- recov %>%
  group_by(Region,Site,Disturbanceme,TSD_months,recove) %>% filter (TSD_months<37) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_TSD5

prec1 <- ggplot(data=df_TSD5, aes(y=recove, colour=Site, x=TSD_months)) +geom_point(size = 1.4,alpha=0.8)
prec1+theme_bw()+theme(axis.title=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+scale_x_continuous(name="TSD_months")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18), legend.position = "right",legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Recovery (%)",x = "Months since disturbance")

prec1 <- ggplot(data=df_TSD5, aes(y=recove, colour=Site, x=TSD_months)) +geom_point(size = 1.4,alpha=0.8)
prec1+theme_bw()+theme(axis.title=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=22)+scale_x_continuous(name="TSD_months")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18), legend.position = "right",legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Recovery (%)",x = "Months since disturbance")

prec2<-ggplot(df_TSD5, aes(x=TSD_months, y=recove, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "black")+geom_line(linetype = 2)+
  theme_bw()+theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))+ labs(x="Months since disturbance", y="Recovery (%, post/cyclone)")+
  guides(legend.key.width=22) + scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.justification = c("right", "top"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
prec2 + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

####New Resilience Data####

df_Res <- rec %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 61) %>%
  group_by(Region,Area,Site,TSD_months,Resilience) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_Res

p1_Res <- ggplot(data=df_Res, aes(x=TSD_months, fill=Region)) +geom_density(adjust=1.5,alpha=.6) 
p1_Res+theme_classic2()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=22),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")+scale_color_brewer(palette="BrBG")

p1a_Res <- ggplot(data=df_Res, aes(x=TSD_months, fill=Site)) +geom_freqpoly(aes(colour = Site), binwidth = 0.1, adjust=1.5,alpha=.2, na.rm = TRUE) 
p1a_Res+theme_classic2()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=22),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")+scale_color_brewer(palette="BrBG")

library(tidytext)

levels(recmeta$Fraction)
df_Res2 <- recmeta %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 37)%>% filter(Resilience!="NA") %>% filter(Site!="San Felipe") %>% filter(Site!="Wooroonooran Basalt") %>% filter(Site!="Wooroonooran Schist") %>%filter(Disturbanceme!="Luis")%>%filter(Disturbanceme!="Ivor")%>%
  group_by(Region,Area,Site,TSD_months,Disturbanceme,Resilience, ResilienceS)%>% 
  mutate (Res = log(Resilience), ResS = log(ResilienceS))
df_Res2

p2_Res<-ggplot(df_Res2, aes(x=TSD_months, y=Res, color = Other_soil_P))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "darkgray")+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14))+ labs(x="Months since disturbance", y="Resilience (ln post/pre)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1, alpha = .8)+ scale_color_brewer(palette = "Paired")+
   theme(legend.title = element_text (size = 20), legend.position=c(0.92,.2),legend.text = element_text (size = 18),legend.justification = c("right", "top"),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p2_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p2a_Res<-ggplot(df_Res2, aes(x=TSD_months, y=Res, color = Other_soil_P))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=T, color = "black", alpha=.7)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14))+ labs(x="Months since disturbance", y="Resilience (ln post/pre)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.position=c(0.6,.22),legend.text = element_text (size = 16),legend.justification = c("right", "top"),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p2a_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p2a_ResS<-ggplot(df_Res2, aes(x=TSD_months, y=ResS, color = Other_soil_P))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'ga',se=T, color = "black", alpha=.7)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14))+ labs(x="Months since disturbance", y="Resilience season (ln post/pre same month)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.position=c(0.6,.22),legend.text = element_text (size = 16),legend.justification = c("right", "top"),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p2a_ResS + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")


p1a_Res <- ggplot(data=df_Res2, aes(x=TSD_months, fill=Region)) +geom_density(adjust=1.5,alpha=.6) 
p1a_Res+theme_classic2()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=22),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=22))+guides(fill = guide_legend(title = "Region"),legend.key.width=32)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.9, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")+scale_color_brewer(palette="BrBG")

df_Res3 <- recmeta %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 37)%>% filter(ResilienceS!="NA") %>% filter(Site!="San Felipe") %>% filter(Site!="Wooroonooran Basalt") %>% filter(Site!="Wooroonooran Schist") %>%filter(Disturbanceme!="Luis")%>%filter(Disturbanceme!="Ivor")%>%
  group_by(Region,Area,Site,TSD_months,Disturbanceme,Resilience, ResilienceS,Other_soil_P) %>%
  summarise(counts = n())%>% mutate (ResS = log(ResilienceS), Res=log(Resilience), DeltaRes=Res-ResS)
df_Res3

#Resilience plot 36 months

p3_Res<-ggplot(df_Res3, aes(x=TSD_months, y=ResS, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "black")+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience (%, post/pre same month)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=100), lty=2, color = "darkgreen", cex=1, alpha = .8)+ scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.justification = c("right", "top"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p3_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

p3_ResD<-ggplot(df_Res3, aes(x=TSD_months, y=DeltaRes))+geom_point(size = 3, aes(colour=Other_soil_P))+geom_line(linetype = 2)+
  geom_area(mapping = aes(y = ifelse(DeltaRes<0 , DeltaRes, 0)), fill = "orange", alpha=0.5)+
  geom_area(mapping = aes(y = ifelse(DeltaRes>0 , DeltaRes, 0)), fill = "green", alpha=0.5)+
  #geom_area(data = pos, fill = "#3DA4AB", alpha=0.3) +
  #geom_area(data = neg, fill = "tomato", alpha=0.3)+#geom_smooth(method = 'glm',se=FALSE, color = "darkgray", alpha=0.3)
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta Resilience (No season - Season)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.7, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.5,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
p3_ResD + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")


df_Res4 <- df_Res3 %>% filter(TSD_months < 13) %>% mutate(Res = log(Resilience), ResS=log(ResilienceS), DeltaRes = Res - ResS) %>% filter(Site!="El Verde FF" & Site!="El Verde DR")
str(df_Res4)

p4_Res<-ggplot(df_Res4, aes(x=TSD_months, y=Res, colour = Disturbanceme, shape = ))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "darkgray")+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience (ln post/pre)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+ scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="right")+scale_x_continuous(breaks=c(0,2,4,6,8,10,12))
#scale_y_log10()
p4_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

p5_Res<-ggplot(df_Res4, aes(x=TSD_months, y=ResS, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "darkgray")+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience (ln post/pre same month)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+ scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="right")+scale_x_continuous(breaks=c(0,2,4,6,8,10,12))
#scale_y_log10()
p5_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

p6_Res<-ggplot(df_Res4, aes(x=TSD_months, y=DeltaRes, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=FALSE, color = "darkgray")+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta resilience (Resilience - ResilienceS)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+ scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="right")+scale_x_continuous(breaks=c(0,2,4,6,8,10,12))
#scale_y_log10()
p6_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

df_Res5 <- recmeta %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 13) %>% filter(ResilienceS!="NA") %>% filter(Site!="San Felipe") %>% filter(Site!="Wooroonooran Basalt" & Site!="Wooroonooran Schist") %>%filter(Disturbanceme!="Luis" & Disturbanceme!="Ivor")%>%
  mutate(Res = log(Resilience), ResS=log(ResilienceS), DeltaRes = Res - ResS) %>% filter(Site!="El Verde FF" & Site!="El Verde DR")
df_Res5

p7_Res<-ggplot(df_Res5, aes(x=TSD_months, y=DeltaRes, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_smooth(method = 'glm',se=T, color = "darkgray", alpha=.7)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta resilience (Resilience - ResilienceS)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="right")+scale_x_continuous(breaks=c(0,2,4,6,8,10,12))
#scale_y_log10()
p7_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))

df_Res5$Other_soil_P

p8_Res<-ggplot(df_Res5, aes(x=TSD_months, y=DeltaRes, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.7)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta resilience (Resilience - ResilienceS)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="right")+scale_x_continuous(breaks=c(0,2,4,6,8,10,12))
#scale_y_log10()
p8_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p9_Res<-ggplot(df_Res5, aes(x=TSD_months, y=Res, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience (ln post/pre)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.88,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(2,4,6,8,10,12))
#scale_y_log10()
p9_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p10_Res<-ggplot(df_Res5, aes(x=TSD_months, y=ResS, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience season (ln post/pre same month)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.88,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(2,4,6,8,10,12))
#scale_y_log10()
p10_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p11_ResD<-ggplot(df_Res5, aes(x=TSD_months, y=DeltaRes, color = Other_soil_P))+geom_point(size = 3)+#geom_line(linetype = 2)+
  geom_area(fill="#69b3a2", alpha=0.4)+#geom_smooth(method = 'glm',se=FALSE, color = "darkgray", alpha=0.3)
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta Resilience (Resilience - ResilienceS)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.88,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(2,4,6,8,10,12))
#scale_y_log10()
p11_ResD + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

pos <- mutate(df_Res5, DeltaRes = ifelse(DeltaRes >= 0, DeltaRes, 0.0001))
neg <- mutate(df_Res5, DeltaRes = ifelse(DeltaRes < 0, DeltaRes, -0.0001))

p12_ResD<-ggplot(df_Res5, aes(x=TSD_months, y=DeltaRes))+geom_point(size = 3, aes(colour=Other_soil_P))+geom_line(linetype = 2)+
  geom_area(mapping = aes(y = ifelse(DeltaRes<0 , DeltaRes, 0)), fill = "orange", alpha=0.5)+
  geom_area(mapping = aes(y = ifelse(DeltaRes>0 , DeltaRes, 0)), fill = "green", alpha=0.5)+
#geom_area(data = pos, fill = "#3DA4AB", alpha=0.3) +
  #geom_area(data = neg, fill = "tomato", alpha=0.3)+#geom_smooth(method = 'glm',se=FALSE, color = "darkgray", alpha=0.3)
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta Resilience (No season - Season)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.7, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.88,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(2,4,6,8,10,12))
p12_ResD + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

df_Res6 <- recmeta %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 23) %>% filter(ResilienceS!="NA") %>% filter(Site!="San Felipe") %>% filter(Site!="Wooroonooran Basalt" & Site!="Wooroonooran Schist") %>%filter(Disturbanceme!="Luis" & Disturbanceme!="Ivor")%>%
  mutate(Res = log(Resilience), ResS=log(ResilienceS), DeltaRes = Res - ResS)
df_Res6

p9a_Res<-ggplot(df_Res6, aes(x=TSD_months, y=Res, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience (ln post/pre)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.5,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,22))
#scale_y_log10()
p9a_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p10a_Res<-ggplot(df_Res6, aes(x=TSD_months, y=ResS, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Resilience season (ln post/pre same month)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.5,.1),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,22))
p10a_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p12a_ResD<-ggplot(df_Res6, aes(x=TSD_months, y=DeltaRes))+geom_point(size = 3, aes(colour=Other_soil_P))+geom_line(linetype = 2)+
  geom_area(mapping = aes(y = ifelse(DeltaRes<0 , DeltaRes, 0)), fill = "orange", alpha=0.5)+
  geom_area(mapping = aes(y = ifelse(DeltaRes>0 , DeltaRes, 0)), fill = "green", alpha=0.5)+
  #geom_area(data = pos, fill = "#3DA4AB", alpha=0.3) +
  #geom_area(data = neg, fill = "tomato", alpha=0.3)+#geom_smooth(method = 'glm',se=FALSE, color = "darkgray", alpha=0.3)
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta Resilience (No season - Season)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.7, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position=c(0.52,.12),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,22))
p12a_ResD + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")


prega<-ggplot(metaregplot, aes(x=pred.pred, y=data_es0ia.yi,color=data_es0ia.Other_soil_P))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=20))+ labs(x="Predicted Effect Size", y="Observed Effect Size")

prega + geom_label_repel(aes(label=data_es0ia$Basin,fill = factor(data_es0ia$Basin)), color = 'black', size = 5, alpha=0.8) +scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))


##Differences between Resilience and ResilienceSeasonal

df_Res4 <- rec %>% filter(Fraction=="TotLitfall") %>% filter (TSD_months < 37)%>% filter(ResilienceS!="NA") %>% filter(Site!="San Felipe") %>% filter(Site!="Wooroonooran Basalt") %>% filter(Site!="Wooroonooran Schist") %>%filter(Disturbanceme!="Luis")%>%filter(Disturbanceme!="Ivor")%>%
  group_by(Region,Area,Site,TSD_months,Disturbanceme,Resilience, ResilienceS) %>%
  summarise(counts = n())%>% mutate (Res = (Resilience)*100, ResS = (ResilienceS)*100, deltaRes = Res - ResS, deltaRes2 = ResS - Res)
df_Res4

p5_Res<-ggplot(df_Res4, aes(x=TSD_months, y=deltaRes, colour = Disturbanceme))+geom_point(size = 3, alpha = 1)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Delta Resilience (%)")+
  guides(legend.key.width=22)+geom_hline(aes(yintercept=0), lty=2, color = "darkgreen", cex=1, alpha = .8)+ scale_color_brewer(palette = "Paired")+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.justification = c("right", "top"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p5_Res + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))


####New Recovery Data####

names(recmeta)
#Recovery = post/cyclone

df_Rec <- recmeta %>%filter(Recovery!="NA") %>% mutate(Cyclone=Disturbanceme)%>%
  group_by(Region,Area,Site,Cyclone,TSD_months,Recovery) %>%mutate(rec = Recovery*100, rec2=log(Recovery)) %>% arrange(desc(TSD_months))
str(df_Rec)
levels(df_Rec$Site)

p1_Rec<-ggplot(df_Rec, aes(x=TSD_months, y=rec, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Recovery (%, post/cyclone)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+#geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="top",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
#scale_y_log10()
p1_Rec + facet_wrap(~Site+Region, scales = "free_y")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p2_Rec<-ggplot(df_Rec, aes(x=TSD_months, y=rec2, color = Other_soil_P))+geom_point(size = 3)+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Recovery (ln post/cyclone)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+#geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18),legend.position="top",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))

p2_Rec + facet_wrap(~Site+Region, scales = "free")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

p3_Rec<-ggplot(df_Rec, aes(x=TSD_months, y=rec2, fill = Cyclone))+geom_point(size = 3, aes(color=Other_soil_P))+geom_smooth(method = 'glm',se=T, color = "black", alpha=0.5)+geom_line(linetype = 2)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=13))+ labs(x="Months since disturbance", y="Recovery (ln post/cyclone)")+
  guides(legend.key.width=22,color = guide_legend(title = "Total soil P (mg/kg)"))+#geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=0.9, alpha = .8)+
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 18),legend.position="top",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
p3_Rec + facet_wrap(~Site+Region, scales = "free")+theme(strip.text.x = element_text(size=12))+scale_color_gradient(low="blue", high="red")

df_Rec2 <- recmeta %>%filter(Recovery!="NA") %>% mutate(Cyclone=Disturbanceme)%>% filter(TSD_months<23) %>%
  group_by(Region,Area,Site,Cyclone,TSD_months,Recovery) %>%mutate(rec = Recovery*100, rec2=log(Recovery)) %>% arrange(desc(TSD_months))
str(df_Rec2)

df_Rec_B<-df_Rec2 %>% filter(Site=="Bisley") %>% mutate(logrec=log10(Recovery))

lmB<-lm(rec~TSD_months, data=df_Rec_B)
summary(lmB)#intercept = 0.0419289, slope = 0.0039659 for 36 months
            #intercept = 0.016761, slope = 0.006151 for 22 months
fitB <- lm(rec~poly(TSD_months,2,raw=TRUE), data=df_Rec_B)
anova(lmB,fitB)

df_Rec_CCJ<-df_Rec2 %>% filter(Site=="Chamela-Cuixmala") %>% filter(Cyclone=="Jova")

lmCCJ<-lm(rec~TSD_months, data=df_Rec_CCJ)
summary(lmCCJ)#intercept = 12.1429, slope = 0.1071 22 months

df_Rec_CCP<-df_Rec2 %>% filter(Site=="Chamela-Cuixmala") %>% filter(Cyclone=="Patricia")

lmCCP<-lm(rec~TSD_months, data=df_Rec_CCP)
summary(lmCCP)#intercept = 2.8571, slope = 0.3929 22 months

df_Rec_EP<-df_Rec2 %>% filter(Site=="East Peak")

lmEP<-lm(rec~TSD_months, data=df_Rec_EP)
summary(lmEP)#intercept = 0.031667, slope = 0.003095 for 22 months

df_Rec_EV<-df_Rec2 %>% filter(Site=="El Verde")

lmEV<-lm(rec~TSD_months, data=df_Rec_EV)
summary(lmEV)#intercept = 0.093197, slope = 0.006269
              #intercept = 0.053224, slope = 0.010431 22 months

df_Rec_EVTD<-df_Rec2 %>% filter(Site=="El Verde TrDeb")

lmEVTD<-lm(rec~TSD_months, data=df_Rec_EVTD)
summary(lmEVTD)#intercept = 0.0204762, slope = 0.0008752
              #intercept = 0.0206494, slope = 0.0008131 22 months

df_Rec_EVT<-df_Rec2 %>% filter(Site=="El Verde Trim")

lmEVT<-lm(rec~TSD_months, data=df_Rec_EVT)
summary(lmEVT)#intercept = 0.0194444, slope = 0.0004955
              #int = 0.0211688, slope = 0.0002541

####Recovery Slope and Soil Data####

recslope<-rbind(data.frame(group = "Bisley", region="Caribbean",slope = 0.0039659, soilP = 320,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.25, soilP=330,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde", region="Caribbean", slope =0.006269, soilP = 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="East Peak", region="Caribbean", slope = 0.0015385, soilP= 260,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="El Verde TrDeb", region="Caribbean", slope = 0.0008752, soilP= 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde Trim", region="Caribbean", slope = 0.0004955, soilP=340,
                            row.names=FALSE, stringsAsFactors=TRUE))

recslope

recslope_22m<-rbind(data.frame(group = "Bisley", region="Caribbean",slope = 0.006151, soilP = 320,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.1071, soilP=330,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.3929, soilP=330,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde", region="Caribbean", slope =0.010431, soilP = 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="East Peak", region="Caribbean", slope = 0.003095, soilP= 260,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="El Verde TrDeb", region="Caribbean", slope = 0.0008131, soilP= 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde Trim", region="Caribbean", slope = 0.0002541, soilP=340,
                           row.names=FALSE, stringsAsFactors=TRUE))

recslope_22m

precslope<-ggplot(recslope, aes(x=soilP, y=slope))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=21),
        axis.text=element_text(size=20))+ labs(y="Recovery slope", x="Total soil P (mg/kg)")

precslope + geom_label_repel(aes(label=group,fill = group), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

precslope + geom_label_repel(aes(label=region,fill = region), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

precslope_22m<-ggplot(recslope_22m, aes(x=soilP, y=slope))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=20))+ labs(y="Recovery slope 1-22 months", x="Total soil P (mg/kg)")

precslope_22m + geom_label_repel(aes(label=group,fill = group), color = 'black', size = 4.5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

precslope_22m + geom_label_repel(aes(label=region,fill = region), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

sl1<-lm(slope~soilP, data=recslope)
summary(sl1)

sl2<-lm(slope~soilP, data=recslope_22m)
summary(sl2)

####Resilience slope soil P####

str(df_Res5)

df_Rec_B2<-df_Res6 %>% filter(Site=="Bisley")

lmB2<-lm(Res~TSD_months, data=df_Rec_B2)
summary(lmB2)#intercept = -1.94907, slope = 0.09234
            #int = -1.75940, slope = 0.06107

lmB2a<-gam(ResS~TSD_months, data=df_Rec_B2)
summary(lmB2a)#intercept = -1.52254, slope = 0.04568

acf(residuals(lmB2), lag=40)
acf(residuals(lmB2a), lag=40)

df_Rec_CC2<-df_Res6 %>% filter(Site=="Chamela-Cuixmala") %>% filter(Disturbanceme=="Jova")

lmCC2<-lm(Res~TSD_months, data=df_Rec_CC2)#Jova
summary(lmCC2)#intercept = -0.07526, slope = -0.01875

lmCC2a<-lm(ResS~TSD_months, data=df_Rec_CC2)#Jova
summary(lmCC2a)#intercept = -0.80956, slope = 0.13854

df_Rec_CC2P<-df_Res6 %>% filter(Site=="Chamela-Cuixmala") %>% filter(Disturbanceme=="Patricia")

lmCC2P<-lm(Res~TSD_months, data=df_Rec_CC2a)#Patricia
summary(lmCC2P)#intercept = -1.34606, slope = 0.12901

lmCC2Pa<-lm(ResS~TSD_months, data=df_Rec_CC2a)#Patricia
summary(lmCC2Pa)#intercept = -1.9671  slope = 0.2864

df_Rec_EP2<-df_Res6 %>% filter(Site=="East Peak")

lmEP2<-lm(Res~TSD_months, data=df_Rec_EP2)
summary(lmEP2)#intercept = -1.416716, slope = 0.049183

lmEP2a<-lm(ResS~TSD_months, data=df_Rec_EP2)
summary(lmEP2a) #intercept = -1.57952, slope = 0.06507

df_Rec_EV2<-df_Res6 %>% filter(Site=="El Verde")

lmEV2<-lm(Res~TSD_months, data=df_Rec_EV2)
summary(lmEV2)#intercept = -1.27713, slope = 0.04396

lmEV2a<-lm(ResS~TSD_months, data=df_Rec_EV2)
summary(lmEV2a) #intercept = -0.81070 slope = 0.02338

df_Rec_EVTD2<-df_Res6 %>% filter(Site=="El Verde TrDeb")

lmEVTD2<-lm(Res~TSD_months, data=df_Rec_EVTD2)
summary(lmEVTD2)#intercept = -1.13272, slope = 0.02891

lmEVTD2a<-lm(ResS~TSD_months, data=df_Rec_EVTD2)
summary(lmEVTD2a) #intercept = -1.67272 slope = 0.04633

df_Rec_EVT2<-df_Res6 %>% filter(Site=="El Verde Trim")

lmEVT2<-lm(Res~TSD_months, data=df_Rec_EVT2)
summary(lmEVT2)#intercept = -1.56594, slope = 0.01462

lmEVT2a<-lm(ResS~TSD_months, data=df_Rec_EVT2)
summary(lmEVT2a) #intercept = -1.75528 slope = 0.03164

df_Rec_BV<-df_Res6 %>% filter(Site=="Bisley_Verde")

lmBV<-lm(Res~TSD_months, data=df_Rec_BV)
summary(lmBV)#intercept = 0.89031, slope = -0.07711

lmBVa<-lm(ResS~TSD_months, data=df_Rec_BV)
summary(lmBVa) #intercept = 1.07116 slope = -0.06350

df_Rec_EPFF<-df_Res6 %>% filter(Site=="East Peak FF")

lmEPFF<-lm(Res~TSD_months, data=df_Rec_EPFF)
summary(lmEPFF)#intercept = -1.57907, slope = 0.06179

lmEPFFa<-lm(ResS~TSD_months, data=df_Rec_EPFF)
summary(lmEPFFa) #int = -1.73037 slope = 0.07709

df_Rec_Kok<-df_Res6 %>% filter(Site=="Kokee")

lmKok<-lm(Res~TSD_months, data=df_Rec_Kok)
summary(lmKok)#intercept = 0.130955, slope = -0.009229

lmKoka<-lm(ResS~TSD_months, data=df_Rec_Kok)
summary(lmKoka) #int = 0.20003 slope = -0.01432

df_Rec_KokN<-df_Res6 %>% filter(Site=="Kokee N")

lmKokN<-lm(Res~TSD_months, data=df_Rec_KokN)
summary(lmKokN)#intercept = 0.077622, slope = 0.001136

df_Rec_KokNP<-df_Res6 %>% filter(Site=="Kokee NP")

lmKokNP<-lm(Res~TSD_months, data=df_Rec_KokNP)
summary(lmKokNP)#intercept = 0.018685, slope = -0.004072

df_Rec_KokP<-df_Res6 %>% filter(Site=="Kokee P")

lmKokP<-lm(Res~TSD_months, data=df_Rec_KokP)
summary(lmKokP)#intercept = 0.060387, slope = -0.005339

df_Rec_MS<-df_Res6 %>% filter(Site=="Mt Spec")

lmMS<-lm(Res~TSD_months, data=df_Rec_MS)
summary(lmMS)#intercept = -0.91462, slope = 0.05686

lmMSa<-lm(ResS~TSD_months, data=df_Rec_MS)
summary(lmMSa) #int = -0.55641 slope = 0.04319

df_Rec_MSD<-df_Res6 %>% filter(Site=="Mt Spec Disturbed")

lmMSD<-lm(Res~TSD_months, data=df_Rec_MSD)
summary(lmMSD)#intercept = -0.75374, slope = 0.05199

lmMSDa<-lm(ResS~TSD_months, data=df_Rec_MSD)
summary(lmMSDa)#intercept = -0.50775, slope = 0.04467

df_Rec_U<-df_Res6 %>% filter(Site=="Utuado")

lmU<-lm(Res~TSD_months, data=df_Rec_U)
summary(lmU)#intercept = -1.72254, slope = 0.06316

lmUa<-lm(ResS~TSD_months, data=df_Rec_U)
summary(lmUa)#intercept = -1.64998, slope = 0.08504

Res_slope<-rbind(data.frame(group = "Bisley", region="Caribbean",slope = 0.09234, soilP = 320,
                           row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group = "Bisley-El Verde", region="Caribbean",slope = -0.17221, soilP = 320,
                            row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.05513, soilP=330,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde", region="Caribbean", slope =0.009316, soilP = 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="East Peak", region="Caribbean", slope = 0.03509, soilP= 260,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="East Peak FF", region="Caribbean", slope = 0.05587, soilP= 260,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group ="El Verde TrDeb", region="Caribbean", slope = -0.05338, soilP= 340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="El Verde Trim", region="Caribbean", slope = -0.02489, soilP=340,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Kokee", region="Hawaii", slope = 0.02317, soilP=220,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Kokee N", region="Hawaii", slope = -0.005228, soilP=220,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Kokee NP", region="Hawaii", slope = -0.00451, soilP=220,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Kokee P", region="Hawaii", slope = 0.01868, soilP=220,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Mt Spec", region="Australia", slope = -0.06399, soilP=130,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Mt Spec Disturbed", region="Australia", slope = 0.02146, soilP=130,
                           row.names=FALSE, stringsAsFactors=TRUE),
                data.frame(group="Utuado", region="Caribbean", slope = 0.2043, soilP=330,
                           row.names=FALSE, stringsAsFactors=TRUE))

Res_slope

Res_slope2<-rbind(data.frame(group = "Bisley", region="Caribbean",slope = 0.06107, soilP = 320,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Chamela-Cuixmala", region="Mexico", slope =-0.01875, soilP=330,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.12901, soilP=330,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="El Verde", region="Caribbean", slope =0.04396, soilP = 340,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group ="East Peak", region="Caribbean", slope = 0.049183, soilP= 260,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group ="East Peak FF", region="Caribbean", slope = 0.06179, soilP= 260,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group ="El Verde TrDeb", region="Caribbean", slope = 0.02891, soilP= 340,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="El Verde Trim", region="Caribbean", slope = 0.01462, soilP=340,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Kokee", region="Hawaii", slope = -0.009229, soilP=220,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Mt Spec", region="Australia", slope = 0.05686, soilP=130,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Mt Spec Disturbed", region="Australia", slope = 0.05199, soilP=130,
                            row.names=FALSE, stringsAsFactors=TRUE),
                 data.frame(group="Utuado", region="Caribbean", slope = 0.06316, soilP=330,
                            row.names=FALSE, stringsAsFactors=TRUE))

Res_slope2
#data.frame(group = "Bisley-El Verde", region="Caribbean",slope = -0.07711, soilP = 320,
           #row.names=FALSE, stringsAsFactors=TRUE),data.frame(group="Kokee N", region="Hawaii", slope = 0.001136, soilP=220,row.names=FALSE, stringsAsFactors=TRUE),
#data.frame(group="Kokee NP", region="Hawaii", slope = -0.004072, soilP=220,
           #row.names=FALSE, stringsAsFactors=TRUE),
#data.frame(group="Kokee P", region="Hawaii", slope = -0.005339, soilP=220,
           #row.names=FALSE, stringsAsFactors=TRUE),

Res_slope2S<-rbind(data.frame(group = "Bisley", region="Caribbean",slope = 0.04568, soilP = 320,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.13854, soilP=330,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Chamela-Cuixmala", region="Mexico", slope =0.2864, soilP=330,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="El Verde", region="Caribbean", slope =0.02338, soilP = 340,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group ="East Peak", region="Caribbean", slope = 0.06507, soilP= 260,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group ="East Peak FF", region="Caribbean", slope = 0.07709, soilP= 260,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group ="El Verde TrDeb", region="Caribbean", slope = 0.04633, soilP= 340,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="El Verde Trim", region="Caribbean", slope = 0.03164, soilP=340,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Kokee", region="Hawaii", slope = -0.009229, soilP=220,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Mt Spec", region="Australia", slope = 0.04319, soilP=130,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Mt Spec Disturbed", region="Australia", slope = 0.04467, soilP=130,
                             row.names=FALSE, stringsAsFactors=TRUE),
                  data.frame(group="Utuado", region="Caribbean", slope = 0.08504, soilP=330,
                             row.names=FALSE, stringsAsFactors=TRUE))

Res_slope2S

#data.frame(group="Kokee N", region="Hawaii", slope = 0.001136, soilP=220,row.names=FALSE, stringsAsFactors=TRUE),
#data.frame(group="Kokee NP", region="Hawaii", slope = -0.004072, soilP=220,
           #row.names=FALSE, stringsAsFactors=TRUE),
#data.frame(group="Kokee P", region="Hawaii", slope = -0.005339, soilP=220,
           #row.names=FALSE, stringsAsFactors=TRUE),

Reslm<-lm(slope~soilP+region, data=Res_slope2)
summary(Reslm)

pRes_slope<-ggplot(Res_slope, aes(x=soilP, y=slope))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=21),
        axis.text=element_text(size=20))+ labs(y="Resilience slope", x="Total soil P (mg/kg)")

pRes_slope + geom_label_repel(aes(label=region,fill = region), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

Reslm<-lm(slope~soilP+region, data=Res_slope2)
summary(Reslm)

pRes_slope2<-ggplot(Res_slope2, aes(x=soilP, y=slope))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=21),
        axis.text=element_text(size=20))+ labs(y="Resilience slope", x="Total soil P (mg/kg)")

pRes_slope2 + geom_label_repel(aes(label=region,fill = region), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

pRes_slope2 + geom_label_repel(aes(label=group,fill = group), color = 'black', size = 4.5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

pRes_slope2S<-ggplot(Res_slope2S, aes(x=soilP, y=slope))+geom_point(size=3.5)+geom_smooth(method = 'lm', col="black", alpha=0.3)+
  theme_bw()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=21),
        axis.text=element_text(size=20))+ labs(y="Resilience season slope", x="Total soil P (mg/kg)")

pRes_slope2S + geom_label_repel(aes(label=region,fill = region), color = 'black', size = 5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))

pRes_slope2S + geom_label_repel(aes(label=group,fill = group), color = 'black', size = 4.5, alpha=0.8) +#scale_color_gradient(low="blue", high="red")+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 18), legend.text = element_text (size = 16),legend.position = "right",legend.justification = c("right", "top"))


lmRes<-lm(slope~soilP, data=Res_slope2)
summary(lmRes)

lmResS<-lm(slope~soilP, data=Res_slope2S)
summary(lmResS)

p2_3 <- ggplot(data=df_TSD3, aes(x=TSD_months, fill=Holdridge_life_zone)) +geom_density(adjust=1.5,alpha=0.3) 
p2_3+ylim(0,0.075)+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+guides(fill = guide_legend(title = "Holdridge life zone"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 30), legend.position = c(.95, .80),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")

metadat$Dist_type
df_TSD4 <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Disturbance,TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_TSD4

df_TSD4a <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Dist_type,TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_TSD4a

p2_4 <- ggplot(data=df_TSD4, aes(x=TSD_months, fill=Disturbance)) +geom_density(adjust=1.5,alpha=0.2) 
p2_4+ylim(0,0.1)+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=24),axis.text.x = element_text(angle=0, hjust=1,size=24))+guides(fill = guide_legend(title = "Wind disturbance"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 38), legend.text = element_text (size = 34), legend.position = c(.75, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))

p2_4a <- ggplot(data=df_TSD4a, aes(x=TSD_months, fill=Dist_type)) +geom_density(adjust=1.5,alpha=0.3) 
p2_4a+ylim(0,0.1)+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+guides(fill = guide_legend(title = "Type of wind disturbance"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 30), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")


df_TSD5 <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Site,TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_TSD5

p2_5 <- ggplot(data=df_TSD5, aes(x=TSD_months, fill=Site)) +geom_density(adjust=1.5,alpha=0.3) 
p2_5+ ylim(0,0.1) + theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=24),axis.text.x = element_text(angle=0, hjust=1,size=24))+guides(fill = guide_legend(title = "Site"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 30), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")


p2<-ggplot(df_TSD2, aes(x=TSD_months, y=counts))
p2+geom_bar(fill="#0073C2FF",stat="identity") +
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))

p2 <- ggplot(data=df_TSD2, aes(x=TSD_months, fill=TSD_months)) +geom_density(adjust=1.5,alpha=.7) 
p2+theme_pubclean()+theme(axis.title=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+guides(fill = guide_legend(title = "Region"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

#Nutrient flux by time since disturbance

df_nutTSD <- nutmeta %>% filter(Variable=="P"|Variable=="N")%>% filter(Raw_Unit=="mg/g")%>%
  group_by(Region, TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df_nutTSD

pPmgg <- ggplot(data=df_nutTSD, aes(x=TSD_months, group=Region, fill=Region)) +geom_density(adjust=1.5,alpha=0.7) 
pPmgg+theme_pubr()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=0.5,size=20))+guides(fill = guide_legend(title = "Region"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 24), legend.position = c(.85, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")

df_nutTSD2 <- nutmeta %>% filter(Raw_Unit == "mg/m2/day"|Raw_Unit=="mg/g")%>%
  group_by(Region, TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)
df_nutTSD2

pnut2_2 <- ggplot(data=df_nutTSD2, aes(x=TSD_months, group=Region, fill=Region)) +geom_density(adjust=1.5,alpha=.3) 
pnut2_2+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=24),axis.text.x = element_text(angle=0, hjust=1,size=24))+guides(fill = guide_legend(title = "Region"),legend.key.width=20)+ scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 34), legend.text = element_text (size = 30), legend.position = c(.75, .70),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")


df_TSD_0.5 <- metadat %>%
  group_by(Region,TSD_months) %>%
  filter(Cat_TSD_months=="0-0.5") %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(TSD_months))
df_TSD_0.5

p_TSD_0.5 <- ggplot(data=df_TSD_0.5, aes(x=TSD_months, group = Region, fill=Region)) +geom_density(adjust=1.5,alpha=.7) 
p_TSD_0.5+theme_pubclean()+theme(axis.title=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+guides(fill = guide_legend(title = "Region"))+
  theme(legend.title = element_text (size = 26), legend.text = element_text (size = 22), legend.position = c(.25, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

df <- df %>%
  arrange(desc(cut)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)

pk<- model_kdat_Amb_k %>% group_by(group, studies) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(prop = round(counts*100/sum(counts), 1),
                                    lab.ypos = cumsum(prop) - 0.5*prop)
pk
pk %>% ggplot(mapping = aes(x=studies)) + geom_density(adjust=1.5,alpha=.7)+theme_bw()
ggdensity(pk, x = "studies", 
          fill = "#0073C2FF", color = "#0073C2FF",rug = TRUE)

 df_TSD2%>%ggplot( aes(x=TSD_months)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
 
 ggplot(pk, aes(group, studies)) +
   geom_linerange(
     aes(x = group, ymin = 0, ymax = studies), 
     color = "darkgray", size = 1.5)+geom_point(aes(color = group), size = 3)+labs(x = "Number of studies", y = "")+
   theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+
   theme(legend.title = element_blank (), legend.text = element_blank(), legend.key=element_blank())+scale_color_brewer(palette="Dark2")
 
 ggplot(df_TSD, aes(TSD_months, counts)) +
  geom_linerange(
    aes(x = TSD_months, ymin = 0, ymax = counts), 
    color = "lightgray", size = 1.5
  )+
  geom_point(aes(color = TSD_months), size = 2)+
  ggpubr::color_palette("jco")+
  theme_pubclean()

dfex2 <- df_2 %>%
  arrange(desc(Trf_TSD_months)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfex2, 4)

pd2<-ggplot(dfex2, aes(x=Trf_TSD_months, y=counts))
pd2+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))


#Litterfall mass by SSHWS Category
levels(metadat$Trf_TSD_months)
metadat$SSHWS_Cat=as.factor(metadat$SSHWS_Cat)
names(metadat)

df_3 <- metadat%>%
  group_by(Site,SSHWS_Cat, Distance_to_Disturb_km,Disturb_Rainfall_mm) %>%
  summarise(counts = n())%>% arrange(desc(counts))
df_3
df_3$SSHWS_Cat=as.factor(df_3$SSHWS_Cat)

plot_3<-ggplot(df_3, aes(x=reorder(Site,-counts),y=counts, fill = SSHWS_Cat)) 
plot_3+geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=45, hjust=1))+guides(fill = guide_legend(title = "SSHWS"))+labs(x = "Saffir-Simpson Hurricane Wind Scale", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


df_SS <- metadat %>%
  group_by(Region, SSHWS_Cat) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_SS

plot_SSCAT<-ggplot(df_SS, aes(x=SSHWS_Cat,y=counts, fill = Region)) 
plot_SSCAT+geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=0, hjust=1))+guides(fill = guide_legend(title = "Region"))+labs(x = "Saffir-Simpson Hurricane Wind Scale", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

names(metadat)
df_TSD5 <- metadat %>% filter(Fraction=="TotLitfall")%>%filter(TSD_months<61)%>%
  group_by(Country,TSD_months) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_TSD5
rlang::last_error()

p2_5 <- ggplot(data=df_TSD5, aes(x=TSD_months, fill=Site)) +geom_density(adjust=1.5,alpha=0.3) 
p2_5+ylim(0,0.1)+theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=20),axis.text.x = element_text(angle=0, hjust=1,size=20))+guides(fill = guide_legend(title = "Site"),legend.key.width=20)+scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_text (size = 32), legend.text = element_text (size = 30), legend.position = c(.95, .95),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+labs(y = "Density",x = "Months since disturbance")

str(metadat$SSHWS_Cat)

df_SS_2 <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Disturbance, SSHWS_Cat) %>%
  summarise(counts = n())%>% mutate(percent = (counts/sum(counts))*100)%>%
  arrange(desc(counts))
df_SS_2

nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

plot_SSCAT_2<-ggplot(df_SS_2, aes(x=SSHWS_Cat,y=counts, fill = Disturbance)) 
plot_SSCAT_2+geom_bar(stat="identity")+
  scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=0, hjust=1))+guides(fill = guide_legend(title = "Disturbance"))+labs(x = "Saffir-Simpson Hurricane Wind Scale", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+ scale_x_continuous(breaks=c(0, 1, 2, 3, 4,5))



dfex3 <- df_3 %>%
  arrange(desc(SSHWS_Category)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfex3, 4)

pe<-ggplot(dfex3, aes(x=SSHWS_Category, y=counts))
pe+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=0, hjust=1))

#Litterfall mass by Parent Material

levels(metadat$Par_Mat)
str(metadat)

df_4 <- metadat %>% filter(Site!="Tanguro")%>% group_by(RockPClass, Rock_type)%>%dplyr::summarise(counts = dplyr::n())%>%mutate(percent = (counts/sum(counts))*100)
df_4
rlang::last_error()

levels(metadat$USDA_Soil_Order)
metadat$USDA_Soil_Order= factor(metadat$USDA_Soil_Order,levels(metadat$USDA_Soil_Order)[c(8,2,5,7,6,3,4,1)])

levels(metadat$Par_Mat)
metadat$Par_Mat= factor(metadat$Par_Mat,levels(metadat$Par_Mat)[c(8,3,11,6,9,7,2,5,10,13,4,14,1,12)])

library(dplyr)
df_4 <- metadat %>%filter(Site!="Tanguro")%>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Other_soil_P,RockPClass, Rock_type,Par_Mat) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(percent = (counts/sum(counts))*100)
df_4

# Install
install.packages("wesanderson")
# Load
library(wesanderson)
library(RColorBrewer)
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(10, "Dark2"))(nb.cols)

plot_ParMat<-ggplot(df_4, aes(y=counts, x = USDA_Soil_Order, fill=Holdridge_life_zone)) 
plot_ParMat+geom_bar(stat="identity")+scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "Holdridge life zone"))+labs(x = "USDA Soil order", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.5, .96),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

plot_ParMat<-ggplot(df_4, aes(y=counts, x = Par_Mat, fill=USDA_Soil_Order)) 
plot_ParMat+geom_bar(stat="identity")+scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "USDA soil order"))+labs(x = "Soil parent material", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.35, .96),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(0,200,400,600,800))

df_rock <- metadat %>%
  group_by(Rock_type) %>%
  dplyr::summarise(counts = dplyr::n())%>% dplyr::mutate(percent = (counts/sum(counts))*100)%>%arrange(desc(Rock_type))
df_rock

df_4a <- metadat %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb"|Treatment=="Trim")%>% filter(Fraction=="TotLitfall")%>%
  group_by(Holdridge_life_zone, RockPClass,USDA_Soil_Order,TempPrec_ratio_x10.2_Cmmyr,Other_soil_P,Site) %>%
  dplyr::summarise(counts = dplyr::n()) %>% dplyr::mutate(percent = (counts/sum(counts))*100) %>%
  arrange(desc(counts))
df_4a<-df_4a %>% rename(Parent.material = Par_Mat)
df_4a

plot_ParMata<-ggplot(df_4a, aes(y=counts, x = RockPClass, fill=USDA_Soil_Order)) 
plot_ParMata+geom_bar(stat="identity")+scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=60, hjust=1))+guides(fill = guide_legend(title = "USDA soil order"))+labs(x = "Soil parent material", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.24, .97),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

plot_ParMatb<-ggplot(df_4a, aes(y=counts, x = USDA_Soil_Order, fill=Holdridge_life_zone)) 
plot_ParMatb+geom_bar(stat="identity")+scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = mycolors) +
  theme_pubclean()+theme(axis.title=element_text(size=28),
                         axis.text=element_text(size=26),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "Holdridge life zone"))+labs(x = "USDA Soil order", y = "Number of observations")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.4, .96),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))


plot_4a<-ggplot(df_4a, aes(x=Other_soil_P,y=TempPrec_ratio_x10.2_Cmmyr, color = Holdridge_life_zone)) 
plot_4a+geom_point(size=4)+scale_fill_brewer(palette = "Set3")+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=40, hjust=1))+guides(fill = guide_legend(title = "Holdridge life zone"))+labs(x="Total soil P", y="MAT / MAP x 100")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.8,.9),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

plot_4a+geom_point(size=4)+scale_fill_brewer(palette = "Set3")+xlim(0,1000)+
  theme_pubclean()+theme(axis.title=element_text(size=24),
                         axis.text=element_text(size=24),axis.text.x = element_text(angle=30, hjust=1))+guides(fill = guide_legend(title = "Holdridge life zone"))+labs(x="Total soil P", y="MAT / MAP x 100")+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 26), legend.position = c(.8,.9),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))

#Soil Phosphorus plot

ggplot(df_4a, aes(x=RockPClass, y=Other_soil_P,color=RockPClass)) +
geom_linerange(aes(x = RockPClass, ymin = 0, ymax = Other_soil_P), color = "lightgray", size = 1.5)+
  geom_point(aes(color = RockPClass), size = 4)+labs(x = "Parent material P category", y = "Total soil P (mg/kg)")+
  theme_pubclean()+theme(axis.title=element_text(size=28),axis.text.y=element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,size=24))+guides(fill = guide_legend(title = "Soil order"))+
  theme(legend.title = element_text (size = 28), legend.text = element_text (size = 28), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.key = element_blank())+
  scale_y_continuous(breaks=c(0,200,400,600,800,1000,1200, 1400, 1600,1800,2000,2200,2400,2600,2800,3000))

ggplot(df_4a, aes(x=RockPClass, y=Other_soil_P,color=RockPClass)) +geom_boxplot()
  
  #geom_linerange(aes(x = RockPClass, ymin = 0, ymax = Other_soil_P), color = "lightgray", size = 1.5)+
  #geom_point(aes(color = RockPClass), size = 4)+labs(x = "Parent material P category", y = "Total soil P (mg/kg)")+
  #theme_pubclean()+theme(axis.title=element_text(size=28),axis.text.y=element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,size=24))+guides(fill = guide_legend(title = "Soil order"))+
  #theme(legend.title = element_text (size = 28), legend.text = element_text (size = 28), legend.position = c(.99, .99),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.key = element_blank())+
  #scale_y_continuous(breaks=c(0,200,400,600,800,1000,1200, 1400, 1600,1800,2000,2200,2400,2600,2800,3000))


ggplot(df_4a, aes(x=USDA_Soil_Order, y=Other_soil_P,color=Parent.material)) +
  geom_linerange(aes(x = USDA_Soil_Order, ymin = 0, ymax = Other_soil_P), color = "darkgray", size = 1.5)+
  geom_point(aes(color = Parent.material), size = 4)+labs(x = "Soil order", y = "Total soil P (mg/kg)")+
  theme_pubclean()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=24),axis.text.x = element_text(angle=0, hjust=0.5,size=24))+guides(fill = guide_legend(title = "Parent Material"))+
  theme(legend.title = element_text (size = 30), legend.text = element_text (size = 28), legend.position = c(.4, .97),legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.key = element_blank())+ylim(0,100)


dfex4 <- df_4 %>%
  arrange(desc(Par_Mat)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfex4, 4)

pg<-ggplot(dfex4, aes(x=reorder(Par_Mat,-counts), y=counts))
pg+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = counts), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=45, hjust=1))



#Litterfall mass by Soil P
levels(metadat$Trf_TSD_months)
str(metadat)
names(metadat)
metadat$Total_Soil_P..Bicarb....NaOH.Pt.=as.numeric(metadat$Total_Soil_P..Bicarb....NaOH.Pt.)

df_P <- metadat%>%
  group_by(Total_Soil_P..Bicarb....NaOH.Pt.) %>%
  summarise(counts = n())
df_P

dfP<- df_P %>%
  arrange(desc(Total_Soil_P..Bicarb....NaOH.Pt.)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfP, 4)

pf<-ggplot(dfP, aes(x=Total_Soil_P..Bicarb....NaOH.Pt., y=prop))
pf+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = prop), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=0, hjust=1))

#Litterfall mass by Total Soil P
levels(metadat$Trf_TSD_months)
str(metadat)
names(metadat)
metadat$Other_Soil_P=as.numeric(metadat$Other_Soil_P)

df_Pt <- metadat%>%
  group_by(Other_Soil_P) %>%
  summarise(counts = n())
df_Pt

dfPt<- df_Pt %>%
  arrange(desc(Other_Soil_P)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfPt, 4)

pg<-ggplot(dfPt, aes(x=Other_Soil_P, y=prop))
pg+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = prop), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=0, hjust=1))


#Litterfall mass by Distance to disturbance
levels(metadat$Trf_TSD_months)
str(metadat)
names(metadat)
metadat$Other_Soil_P=as.numeric(metadat$Other_Soil_P)

dfD <- metadat%>%
  group_by(Distance_to_Dist_km) %>%
  summarise(counts = n())
dfD

dfDex<- dfD %>%
  arrange(desc(Distance_to_Dist_km)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)
head(dfPt, 4)

pg<-ggplot(dfPt, aes(x=Other_Soil_P, y=prop))
pg+geom_bar(fill="#0073C2FF",stat="identity") +
  geom_text(aes(label = prop), vjust = -0.3, size=7, color="purple")+
  theme_pubclean()+theme(axis.title=element_text(size=20),
                         axis.text=element_text(size=16),axis.text.x = element_text(angle=0, hjust=1))

####PR baseline litterfall plots across forests####

pr<-metadat %>% filter (Area == "Puerto Rico")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde Trim")%>%filter(Site!="El Verde TrDeb")
names(pr)

pr_lit<-metadat %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde DR")
names(pr_lit)

pr_lit2<-metadat %>% filter (Area == "Puerto Rico")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde DR")
names(pr_lit2)

levels(nutmeta$Raw_Unit)
pr_nut<-nutmeta %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable =="N"|Variable=="P")
str(pr_nut)

pr_nut2<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable =="N"|Variable=="P")
str(pr_nut2)

pr_nut2a<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable=="P")
str(pr_nut2a)

pr_nut2b<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/g")%>% filter(Variable=="P")
str(pr_nut2b)

pr_nut3<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable=="N")
str(pr_nut3)

pr_nut3b<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/g")%>% filter(Variable=="N")
str(pr_nut3b)

plot_pr<-ggplot(pr, aes(x=reorder(Site,-Other_soil_P), y=Other_soil_P))+geom_point(alpha=0.5,aes(size=MAT_MAP_x100))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Total soil P (mg/kg)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.8,0.85),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr

plot_pr2<-ggplot(pr, aes(x=reorder(Site,-Other_soil_P), y=Other_soil_P))+geom_point(alpha=0.8,aes(size=MAT_MAP_x100,color=Holdridge_life_zone))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Total soil P (mg/kg)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.7,0.72),legend.background = element_rect(fill=alpha('transparent', 0.8)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_brewer(palette="Dark2")+scale_size(range = c(2, 8))
plot_pr2

plot_pr_lit<-ggplot(pr_lit, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_boxplot(alpha=0.7)+geom_point(alpha=0.5,aes(size=MAT_MAP_x100,color=Other_soil_P))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean total litterfall (g/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.8,0.68),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_lit

plot_pr_lit2<-ggplot(pr_lit2, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.9,aes(color=Other_soil_P,shape=Fraction,size=3))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean litterfall (g/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.85,0.62),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_shape_discrete(solid=F)
plot_pr_lit2

names(pr_nut)
plot_pr_nut<-ggplot(pr_nut, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.5,aes(shape=Variable,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut

plot_pr_nut2<-ggplot(pr_nut2, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,size=Other_soil_P,color=Variable))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_brewer(palette="Dark2")+scale_size(range = c(2, 8))+scale_shape_discrete(solid=F)
plot_pr_nut2

plot_pr_nut2a<-ggplot(pr_nut2a, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean P flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut2a

plot_pr_nut2b<-ggplot(pr_nut2b, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean P in litterfall (mg/g)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut2b

plot_pr_nut3<-ggplot(pr_nut3, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean N flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut3

plot_pr_nut3b<-ggplot(pr_nut3b, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean N in litterfall (mg/g)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))+scale_y_continuous(breaks=c(0,10,15,20,25))
plot_pr_nut3b

pr_C<-nutmeta %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall") %>% filter(Raw_Unit == "mg/m2/day") %>% filter(Variable=="C")
str(pr_C)

plot_pr_C<-ggplot(pr_C, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.5,aes(color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_C

ggplot(metadat, aes(x=reorder(Basin,-StormFrequencyNorm), y=StormFrequencyNorm))+geom_boxplot()
kruskal.test(StormsPerYearSince1955~Basin,data=metadat)
kruskal.test(StormFrequencyNorm~Basin,data=metadat)
kruskal.test(StormFrequencyNorm~Site,data=metadat)

####PR baseline litterfall plots across forests####

pr<-metadat %>% filter (Area == "Puerto Rico")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde Trim")%>%filter(Site!="El Verde TrDeb")
names(pr)

pr_lit<-metadat %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde DR")
names(pr_lit)

pr_lit2<-metadat %>% filter (Area == "Puerto Rico")%>%filter(Site!="East Peak FF")%>%filter(Site!="East Peak FF")%>%filter(Site!="El Verde DR")%>%filter(Site!="El Verde FF")%>%filter(Site!="El Verde DR")
names(pr_lit2)

levels(nutmeta$Raw_Unit)
pr_nut<-nutmeta %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable =="N"|Variable=="P")
str(pr_nut)

pr_nut2<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable =="N"|Variable=="P")
str(pr_nut2)

pr_nut2a<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable=="P")
str(pr_nut2a)

pr_nut2b<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/g")%>% filter(Variable=="P")
str(pr_nut2b)

pr_nut3<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/m2/day")%>% filter(Variable=="N")
str(pr_nut3)

pr_nut3b<-nutmeta %>% filter (Area == "Puerto Rico") %>% filter(Raw_Unit == "mg/g")%>% filter(Variable=="N")
str(pr_nut3b)

plot_pr<-ggplot(pr, aes(x=reorder(Site,-Other_soil_P), y=Other_soil_P))+geom_point(alpha=0.5,aes(size=MAT_MAP_x100))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Total soil P (mg/kg)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.8,0.85),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr

plot_pr2<-ggplot(pr, aes(x=reorder(Site,-Other_soil_P), y=Other_soil_P))+geom_point(alpha=0.8,aes(size=MAT_MAP_x100,color=Holdridge_life_zone))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Total soil P (mg/kg)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.7,0.72),legend.background = element_rect(fill=alpha('transparent', 0.8)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_brewer(palette="Dark2")+scale_size(range = c(2, 8))
plot_pr2

plot_pr_lit<-ggplot(pr_lit, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_boxplot(alpha=0.7)+geom_point(alpha=0.5,aes(size=MAT_MAP_x100,color=Other_soil_P))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean total litterfall (g/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.8,0.68),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_lit

plot_pr_lit2<-ggplot(pr_lit2, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.9,aes(color=Other_soil_P,shape=Fraction,size=3))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean litterfall (g/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position=c(0.85,0.62),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_shape_discrete(solid=F)
plot_pr_lit2

names(pr_nut)
plot_pr_nut<-ggplot(pr_nut, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.5,aes(shape=Variable,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut

plot_pr_nut2<-ggplot(pr_nut2, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,size=Other_soil_P,color=Variable))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_brewer(palette="Dark2")+scale_size(range = c(2, 8))+scale_shape_discrete(solid=F)
plot_pr_nut2

plot_pr_nut2a<-ggplot(pr_nut2a, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean P flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut2a

plot_pr_nut2b<-ggplot(pr_nut2b, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean P in litterfall (mg/g)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut2b

plot_pr_nut3<-ggplot(pr_nut3, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean N flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_nut3

plot_pr_nut3b<-ggplot(pr_nut3b, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.8,aes(shape=Fraction,color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=0.5),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean N in litterfall (mg/g)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))+scale_y_continuous(breaks=c(0,10,15,20,25))
plot_pr_nut3b

pr_C<-nutmeta %>% filter (Area == "Puerto Rico")%>% filter (Fraction=="TotLitfall") %>% filter(Raw_Unit == "mg/m2/day") %>% filter(Variable=="C")
str(pr_C)

plot_pr_C<-ggplot(pr_C, aes(x=reorder(Site,-Pre_Mean), y=Pre_Mean))+geom_point(alpha=0.5,aes(color=Other_soil_P,size=MAT.MAPx100_Cmmyr))+
  theme_bw()+coord_flip()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=17,angle=0, hjust=1),axis.text.y=element_text(size=17, hjust=1))+ labs(x="PR forest", y="Mean element flux in litterfall (mg/m2/day)")+theme(legend.title = element_text (size = 16), legend.text = element_text (size = 14),legend.position="right",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.box.background = element_rect(colour = "black"))+
  scale_color_gradient(low="blue", high="orange")+scale_size(range = c(2, 8))
plot_pr_C

ggplot(metadat, aes(x=reorder(Basin,-StormFrequencyNorm), y=StormFrequencyNorm))+geom_boxplot()
kruskal.test(StormsPerYearSince1955~Basin,data=metadat)
kruskal.test(StormFrequencyNorm~Basin,data=metadat)
kruskal.test(StormFrequencyNorm~Site,data=metadat)
