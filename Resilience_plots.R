##Figure Resilience Observations All Fractions#

#Packages####
library(metafor)
library(tidyverse)
library(ggplot2)
library(patchwork)

####Uploading Litterfall mass data####
metadat<-read.csv(file.choose())#Litterfall_Mass

##Filtering data to retain only points used to calculate the resilience metric
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
str(rec)
unique(levels(as.factor(rec$Case_study)))
#Including litterfall data from 1 to 36 months post-disturbance
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)
str(res_all)
#Including Ambient conditions only
res_amb<-res_all %>% filter(Treatment=="Ambient")
str(res_amb)#948
#Including Ambient and TrimDeb treatment from CTE
res_amb2<-res_all %>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")
str(res_amb2)#1128

#Creating new column named "Case_study"
res_amb$Case_study= paste(res_amb$Site, res_amb$DisturbanceName, sep="|")
res_amb2$Case_study= paste(res_amb2$Site, res_amb2$DisturbanceName, sep="|")

###Calculating Effect sizes (Hedge's g)####
#1 to 36 months post-disturbance - Ambient conditions only
data_esall_amb <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_amb2, measure = "ROM")
str(data_esall_amb)#945 observations including all litterfall mass fractions

#Including CTE TrimDeb
data_esall_amb2 <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = res_amb2, measure = "ROM")
str(data_esall_amb2)#1125 observations

##Total litterfall mass resilience [ln(post/pre)]####
#Months 1 to 36 post-disturbance
tot_amb<-data_esall_amb %>% filter(Fraction=="TotLitfall")
str(tot_amb)#249 observations - ambient conditions only, 285 obs with CTE
unique(levels(as.factor(tot_amb$Case_study)))

#Checking number of case studies per geographic region and country
summary(tot_amb$Region)
summary(tot_amb$Country)
#Number of case studies in Puerto Rico is 108 out of 249 (43%) when CTE is not included

#Data grouping and wrangling to prepare figures
dres_amb_tot <- tot_amb %>%group_by(yi,vi,Case_study,Treatment,TSD_months)  #%>%dplyr::summarise(counts = dplyr::n())
dres_amb_tot
Obs_res_tot<- cbind(data.frame(Obs=tot_amb$yi,Se=sqrt(tot_amb$vi),Months=factor(tot_amb$TSD_months),
                                  SoilP=tot_amb$Other_soil_P,fac_SoilP=factor(tot_amb$Other_soil_P),
                                  Case_study=tot_amb$Case_study,Cyclone_freq=tot_amb$StormFrequencyNorm,Country=tot_amb$Country))
str(Obs_res_tot)
Obs_res_tot$Months
Obs_res_tot$Months<-factor(Obs_res_tot$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
levels(Obs_res_tot$Months)

##Data wrangling Including CTE TrimDeb
data_esall_amb_final2<-data_esall_amb2 %>% filter(Fraction!="Misc fall")
str(data_esall_amb_final2)

Se=sqrt(data_esall_amb_final2$vi)
Se
data_esall_amb_final2$yi
Obs=data_esall_amb_final2$yi
Obs
data_esall_amb_final2$Months<-factor(data_esall_amb_final2$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                                  "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
Obs_res_all2<- cbind(data.frame(Obs=data_esall_amb_final2$yi,Se=sqrt(data_esall_amb_final2$vi),Months=factor(data_esall_amb_final2$TSD_months),
                               SoilP=data_esall_amb_final2$Other_soil_P,Fraction=data_esall_amb_final2$Fraction,fac_SoilP=factor(data_esall_amb_final2$Other_soil_P),
                               Case_study=data_esall_amb_final2$Case_study,Cyclone_freq=data_esall_amb_final2$StormFrequencyNorm,Country=data_esall_amb_final2$Country))

##Calculating Overall Resilience for certain months post-disturbance

#dataframe
data_es_tot_amb<-data_esall_amb_final2 %>% filter(Fraction=="TotLitfall")
str(data_es_tot_amb)#285 observations
summary(data_es_tot_amb$Treatment)
#Number of case studies in Puerto Rico is 144 out of 285 (50%) when CTE is included

#1 month post-disturbance for ambient + TrimDeb CTE####
data_es_tot_amb_1<-data_es_tot_amb %>% filter(TSD_months==1)
names(data_es_tot_amb_1)
#If negative, this is the formula: -(1-exp(random1$b))*100
#If positive, this is the formula: ((exp(random0a$b))-1)*100

#running mixed-effects meta-analysis model for one month post-disturbance
tot_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_1,struct = "HAR",method = "REML")
summary(tot_meta_1)
tot_meta_1$b

#same for new resilience effect size and variance
tot_meta_1_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_1,struct = "HAR",method = "REML")
summary(tot_meta_1_new)
-(1-exp(tot_meta_1_new$b))*100
((exp(tot_meta_1_new$se))-1)*100

#3 months TF####
data_es_tot_amb_3<-data_es_tot_amb %>% filter(TSD_months==3)
tot_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_3,struct = "HAR",method = "REML")
summary(tot_meta_3)
coef(summary(tot_meta_3))

#same for new resilience effect size and variance
tot_meta_3_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_3,struct = "HAR",method = "REML")
summary(tot_meta_3_new)

#5 months TF####
data_es_tot_amb_5<-data_es_tot_amb %>% filter(TSD_months==5)
tot_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_5,struct = "HAR",method = "REML")
                    #mods=~TSD_months)
summary(tot_meta_5)
-(1-exp(tot_meta_5$b))*100
coef(summary(tot_meta_5))
summary(tot_meta_5)

#same for new resilience effect size and variance
tot_meta_5_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_5,struct = "HAR",method = "REML")
#mods=~TSD_months)
summary(tot_meta_5_new)
-(1-exp(tot_meta_5_new$b))*100
((exp(tot_meta_5_new$se))-1)*100

#8 months TF####
data_es_tot_amb_8<-data_es_tot_amb %>% filter(TSD_months==8)
tot_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_8,struct = "HAR",method = "REML")
summary(tot_meta_8)
coef(summary(tot_meta_8))

#same for new resilience effect size and variance
tot_meta_8_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_8,struct = "HAR",method = "REML")
summary(tot_meta_8_new)

#12 months TF####
data_es_tot_amb_12<-data_es_tot_amb %>% filter(TSD_months==12)
tot_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_12,struct = "HAR",method = "REML")
summary(tot_meta_12)

#same for new resilience effect size and variance
tot_meta_12_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_12,struct = "HAR",method = "REML")
summary(tot_meta_12_new)
-(1-exp(tot_meta_12_new$b))*100
((exp(tot_meta_12_new$se))-1)*100

#15 months TF####
data_es_tot_amb_15<-data_es_tot_amb %>% filter(TSD_months==15)
tot_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_15,struct = "HAR",method = "REML")
summary(tot_meta_15)

#same for resilience effect size and variance
tot_meta_15_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_15,struct = "HAR",method = "REML")
summary(tot_meta_15_new)

#18 months TF####
data_es_tot_amb_18<-data_es_tot_amb %>% filter(TSD_months==18)
tot_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_18,struct = "HAR",method = "REML")
summary(tot_meta_18)

#same for new resilience effect size and variance
tot_meta_18_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_18,struct = "HAR",method = "REML")
summary(tot_meta_18_new)

#21 monthsTF####
data_es_tot_amb_21<-data_es_tot_amb %>% filter(TSD_months==21)
tot_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_21,struct = "HAR",method = "REML")
summary(tot_meta_21)

#same for new resilience effect size and variance
tot_meta_21_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_21,struct = "HAR",method = "REML")
summary(tot_meta_21_new)

#data frame
tot_pant_res<-rbind(data.frame(Months="1", estimate=tot_meta_1$b,se=tot_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
  data.frame(Months="3", estimate=tot_meta_3$b,se=tot_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                                  data.frame(Months="5", estimate=tot_meta_5$b,se=tot_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                                  data.frame(Months="8", estimate=tot_meta_8$b,se=tot_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                                  data.frame(Months="12", estimate=tot_meta_12$b,se=tot_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=tot_meta_15$b,se=tot_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="18", estimate=tot_meta_18$b,se=tot_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="21", estimate=tot_meta_21$b,se=tot_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))
           

tot_pant_res

#data frame for new resilience metrics (ln(post/pre)/number of months)
tot_pant_res_new<-rbind(data.frame(Months="1", estimate=tot_meta_1_new$b,se=tot_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="3", estimate=tot_meta_3_new$b,se=tot_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="5", estimate=tot_meta_5_new$b,se=tot_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="8", estimate=tot_meta_8_new$b,se=tot_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="12", estimate=tot_meta_12_new$b,se=tot_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=tot_meta_15_new$b,se=tot_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="18", estimate=tot_meta_18_new$b,se=tot_meta_18_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="21", estimate=tot_meta_21_new$b,se=tot_meta_21_new$se,row.names=FALSE, stringsAsFactors=TRUE))


tot_pant_res_new

##Figure6a####
names(data_es_tot_amb)
Fig_res_amb_tot <- ggplot(data_es_tot_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_tot<-Fig_res_amb_tot+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_tot<-Fig_res_amb_tot+geom_pointrange(data=tot_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_tot<-Fig_res_amb_tot+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","black","#665191","#af060f"))
Fig_res_amb_tot<-Fig_res_amb_tot+labs(color="")+guides(size=FALSE)
Fig_res_amb_tot
Fig_res_amb_tot<-Fig_res_amb_tot+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+
   theme(axis.title.x =element_text(vjust = 0.5,size=24),axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.3)))
Fig_res_amb_tot<-Fig_res_amb_tot+annotate("text", x =7, y = 2, fontface="bold",label = "a Total litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_tot<-Fig_res_amb_tot+ylab("Resilience")#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))
Fig_res_amb_tot

#New Figure 6a TF####
Fig_res_tot_new <- ggplot(data_es_tot_amb, aes(y=yi_new, x=Months,group=Case_study))
Fig_res_tot_new<-Fig_res_tot_new+geom_point(aes(group=Case_study,col=Country,size=vi_new),shape=21,stroke=1.3,alpha=0.5)
Fig_res_tot_new<-Fig_res_tot_new+geom_pointrange(data=tot_pant_res_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_tot_new<-Fig_res_tot_new+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","black","#665191","#af060f"))
Fig_res_tot_new<-Fig_res_tot_new+labs(color="")+guides(size=FALSE)+ylim(-1.5,1.5)
Fig_res_tot_new
Fig_res_tot_new<-Fig_res_tot_new+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-1.5,1.5)+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.3)))
Fig_res_tot_new<-Fig_res_tot_new+annotate("text", x =7, y = 1.5, fontface="bold",label = "a Total litterfall", size=8,colour="black")+labs(x="")
Fig_res_tot_new<-Fig_res_tot_new+ylab("Resilience")#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))
Fig_res_tot_new

####Leaf litterfall (LF)####
##Calculating Overall Resilience for certain points in time since disturbance
data_es_leaf_amb<-data_esall_amb_final2 %>% filter(Fraction=="Leaf fall")
str(data_es_leaf_amb)
leaf_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb,struct = "HAR",method = "REML",mods=~TSD_months)
summary(leaf_meta_all)

#1 month LF####
data_es_leaf_amb_1<-data_es_leaf_amb %>% filter(TSD_months==1)
leaf_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_1,struct = "HAR",method = "REML")
summary(leaf_meta_1)
-(1-exp(leaf_meta_1$b))*100 #value in %

#same for new resilience effect size and variance
leaf_meta_1_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leaf_amb_1,struct = "HAR",method = "REML")
summary(leaf_meta_1_new)
-(1-exp(leaf_meta_1_new$b))*100
((exp(leaf_meta_1_new$se))-1)*100

#3 months LF####
data_es_leaf_amb_3<-data_es_leaf_amb %>% filter(TSD_months==3)
leaf_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_3,struct = "HAR",method = "REML")
summary(leaf_meta_3)
-(1-exp(leaf_meta_3$b))*100
#same for new resilience effect size and variance
leaf_meta_3_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_3,struct = "HAR",method = "REML")
summary(leaf_meta_3_new)

#5 months LF####
data_es_leaf_amb_5<-data_es_leaf_amb %>% filter(TSD_months==5)
leaf_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_5,struct = "HAR",method = "REML")
summary(leaf_meta_5)
-(1-exp(leaf_meta_5$b))*100
#same for new resilience effect size and variance
leaf_meta_5_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_5,struct = "HAR",method = "REML")
summary(leaf_meta_5_new)

#8 months LF####
data_es_leaf_amb_8<-data_es_leaf_amb %>% filter(TSD_months==8)
leaf_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_8,struct = "HAR",method = "REML")
summary(leaf_meta_8)
#same for new resilience effect size and variance
leaf_meta_8_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_8,struct = "HAR",method = "REML")
summary(leaf_meta_8_new)

#12 months LF####
data_es_leaf_amb_12<-data_es_leaf_amb %>% filter(TSD_months==12)
leaf_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_12,struct = "HAR",method = "REML")
summary(leaf_meta_12)
#same for new resilience effect size and variance
leaf_meta_12_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb_12,struct = "HAR",method = "REML")
summary(leaf_meta_12_new)
-(1-exp(leaf_meta_12_new$b))*100
((exp(leaf_meta_12_new$se))-1)*100

#15 months LF####
data_es_leaf_amb_15<-data_es_leaf_amb %>% filter(TSD_months==15)
leaf_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_15,struct = "HAR",method = "REML")
summary(leaf_meta_15)
#same for new resilience effect size and variance
leaf_meta_15_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb_15,struct = "HAR",method = "REML")
summary(leaf_meta_15_new)

#18 months LF####
data_es_leaf_amb_18<-data_es_leaf_amb %>% filter(TSD_months==18)
leaf_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_18,struct = "HAR",method = "REML")
summary(leaf_meta_18)

leaf_meta_18_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb_18,struct = "HAR",method = "REML")
summary(leaf_meta_18_new)

#21 months LF####
data_es_leaf_amb_21<-data_es_leaf_amb %>% filter(TSD_months==21)
leaf_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_21,struct = "HAR",method = "REML")
summary(leaf_meta_21)

leaf_meta_21_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb_21,struct = "HAR",method = "REML")
summary(leaf_meta_21_new)

#data frames for figures
leaf_pant_res<-rbind(data.frame(Months="1", estimate=leaf_meta_1$b,se=leaf_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="3", estimate=leaf_meta_3$b,se=leaf_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="5", estimate=leaf_meta_5$b,se=leaf_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="8", estimate=leaf_meta_8$b,se=leaf_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="12", estimate=leaf_meta_12$b,se=leaf_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=leaf_meta_15$b,se=leaf_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="18", estimate=leaf_meta_18$b,se=leaf_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="21", estimate=leaf_meta_21$b,se=leaf_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


leaf_pant_res
#new data frame for new resilience metric
leaf_pant_res_new<-rbind(data.frame(Months="1", estimate=leaf_meta_1_new$b,se=leaf_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=leaf_meta_3_new$b,se=leaf_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=leaf_meta_5_new$b,se=leaf_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=leaf_meta_8_new$b,se=leaf_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=leaf_meta_12_new$b,se=leaf_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=leaf_meta_15_new$b,se=leaf_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=leaf_meta_18_new$b,se=leaf_meta_18_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=leaf_meta_21_new$b,se=leaf_meta_21_new$se,row.names=FALSE, stringsAsFactors=TRUE))


leaf_pant_res_new

##Figure6b
Fig_res_amb_leaf <- ggplot(data_es_leaf_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_leaf<-Fig_res_amb_leaf+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_leaf<-Fig_res_amb_leaf+geom_pointrange(data=leaf_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_leaf<-Fig_res_amb_leaf+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_leaf<-Fig_res_amb_leaf+labs(color="")+guides(size=FALSE)
Fig_res_amb_leaf<-Fig_res_amb_leaf+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+ylab("Resilience")+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
annotate("text", x =6.5, y = 2, fontface="bold",label = "b Leaf litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_leaf<-Fig_res_amb_leaf+xlab("Time since cyclone (Months)")
Fig_res_amb_leaf

#New Figure6b LF####
Fig_res_amb_leaf_new <- ggplot(data_es_leaf_amb, aes(y=yi_new, x=Months,group=Case_study))
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+geom_point(aes(group=Case_study,col=Country,size=vi_new),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+geom_pointrange(data=leaf_pant_res_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+labs(color="")+guides(size=FALSE)
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-2.5,1)+
  ylab("Resilience")+#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =6.5, y = 1, fontface="bold",label = "b Leaf litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_leaf_new<-Fig_res_amb_leaf_new+xlab("Time since cyclone (Months)")
Fig_res_amb_leaf_new

##Wood litterfall (WF)####
##Calculating Overall Resilience for certain points in time since disturbance
data_es_wood_amb<-data_esall_amb_final2 %>% filter(Fraction=="Wood fall")
wood_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                       tdist = TRUE,
                       data = data_es_wood_amb,struct = "HAR",method = "REML",mods=~TSD_months)
summary(wood_meta_all)

#1 month WF####
data_es_wood_amb_1<-data_es_wood_amb %>% filter(TSD_months==1)
wood_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_1,struct = "HAR",method = "REML")
summary(wood_meta_1)

wood_meta_1_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_1,struct = "HAR",method = "REML")
summary(wood_meta_1_new)
-(1-exp(wood_meta_1_new$b))*100
((exp(wood_meta_1_new$se))-1)*100

#3 months WF####
data_es_wood_amb_3<-data_es_wood_amb %>% filter(TSD_months==3)
wood_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_3,struct = "HAR",method = "REML")
summary(wood_meta_3)
-(1-exp(wood_meta_3$b))*100

wood_meta_3_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_3,struct = "HAR",method = "REML")
summary(wood_meta_3_new)

#5 months WF####
data_es_wood_amb_5<-data_es_wood_amb %>% filter(TSD_months==5)
wood_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_5,struct = "HAR",method = "REML")
summary(wood_meta_5)

wood_meta_5_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_5,struct = "HAR",method = "REML")
summary(wood_meta_5_new)

#8 months WF####
data_es_wood_amb_8<-data_es_wood_amb %>% filter(TSD_months==8)
wood_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_8,struct = "HAR",method = "REML")
summary(wood_meta_8)

wood_meta_8_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_8,struct = "HAR",method = "REML")
summary(wood_meta_8_new)

#12 months WF####
data_es_wood_amb_12<-data_es_wood_amb %>% filter(TSD_months==12)
wood_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_12,struct = "HAR",method = "REML")
summary(wood_meta_12)

wood_meta_12_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_12,struct = "HAR",method = "REML")
summary(wood_meta_12_new)
((exp(wood_meta_12_new$se))-1)*100
-(1-exp(wood_meta_12_new$b))*100

#15 months WF####
data_es_wood_amb_15<-data_es_wood_amb %>% filter(TSD_months==15)
wood_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_15,struct = "HAR",method = "REML")
summary(wood_meta_15)

wood_meta_15_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_15,struct = "HAR",method = "REML")
summary(wood_meta_15_new)

#18 months WF####
data_es_wood_amb_18<-data_es_wood_amb %>% filter(TSD_months==18)
wood_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_18,struct = "HAR",method = "REML")
summary(wood_meta_18)

wood_meta_18_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_18,struct = "HAR",method = "REML")
summary(wood_meta_18_new)

#21 months WF####
data_es_wood_amb_21<-data_es_wood_amb %>% filter(TSD_months==21)
wood_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_21,struct = "HAR",method = "REML")
summary(wood_meta_21)

wood_meta_21_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_21,struct = "HAR",method = "REML")
summary(wood_meta_21_new)

#data frame for figures
wood_pant_res<-rbind(data.frame(Months="1", estimate=wood_meta_1$b,se=wood_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=wood_meta_3$b,se=wood_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=wood_meta_5$b,se=wood_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=wood_meta_8$b,se=wood_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=wood_meta_12$b,se=wood_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=wood_meta_15$b,se=wood_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=wood_meta_18$b,se=wood_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=wood_meta_21$b,se=wood_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


wood_pant_res

wood_pant_res_new<-rbind(data.frame(Months="1", estimate=wood_meta_1_new$b,se=wood_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=wood_meta_3_new$b,se=wood_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=wood_meta_5_new$b,se=wood_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=wood_meta_8_new$b,se=wood_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=wood_meta_12_new$b,se=wood_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=wood_meta_15_new$b,se=wood_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=wood_meta_18_new$b,se=wood_meta_18_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=wood_meta_21_new$b,se=wood_meta_21_new$se,row.names=FALSE, stringsAsFactors=TRUE))


wood_pant_res_new

##Figure 6c 
names(data_es_wood_amb)
Fig_res_amb_wood <- ggplot(data_es_wood_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_wood<-Fig_res_amb_wood+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_wood<-Fig_res_amb_wood+geom_pointrange(data=wood_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_wood<-Fig_res_amb_wood+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_wood<-Fig_res_amb_wood+labs(color="")+guides(size=FALSE)
Fig_res_amb_wood
Fig_res_amb_wood<-Fig_res_amb_wood+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
   scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+ylab("Resilience")+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("Time since cyclone (Months)")+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =6.8, y = 2, fontface="bold",label = "c Wood litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_wood

#New Figure6c WF####
Fig_res_amb_wood_new <- ggplot(data_es_wood_amb, aes(y=yi_new, x=Months,group=Case_study))
Fig_res_amb_wood_new<-Fig_res_amb_wood_new+geom_point(aes(group=Case_study,col=Country,size=vi_new),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_wood_new<-Fig_res_amb_wood_new+geom_pointrange(data=wood_pant_res_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_wood_new<-Fig_res_amb_wood_new+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_wood_new<-Fig_res_amb_wood_new+labs(color="")+guides(size=FALSE)
Fig_res_amb_wood_new
Fig_res_amb_wood_new<-Fig_res_amb_wood_new+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-2,1.6)+
  ylab("Resilience")+#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("Time since disturbance (Months)")+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =6, y = 1.5, fontface="bold",label = "c Wood fall", size=8,colour="black")+labs(x="")
Fig_res_amb_wood_new

##FFS fall####

##Calculating Overall Resilience for certain points in time since disturbance
data_es_ffs_amb<-data_esall_amb_final2 %>% filter(Fraction=="FFS fall")
str(data_es_ffs_amb)#173
ffs_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                       tdist = TRUE,
                       data = data_es_ffs_amb,struct = "HAR",method = "REML",mods=~TSD_months)
summary(ffs_meta_all)

#1 month FFS####
data_es_ffs_amb_1<-data_es_ffs_amb %>% filter(TSD_months==1)
ffs_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_1,struct = "HAR",method = "REML")
summary(ffs_meta_1)

ffs_meta_1_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_ffs_amb_1,struct = "HAR",method = "REML")
summary(ffs_meta_1_new)
-(1-exp(ffs_meta_1_new$b))*100
((exp(ffs_meta_1_new$se))-1)*100

#3 months FFS####
data_es_ffs_amb_3<-data_es_ffs_amb %>% filter(TSD_months==3)
ffs_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_3,struct = "HAR",method = "REML")
summary(ffs_meta_3)

ffs_meta_3_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_ffs_amb_3,struct = "HAR",method = "REML")
summary(ffs_meta_3_new)

#5 months FFS####
data_es_ffs_amb_5<-data_es_ffs_amb %>% filter(TSD_months==5)
ffs_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_5,struct = "HAR",method = "REML")
summary(ffs_meta_5)

ffs_meta_5_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_ffs_amb_5,struct = "HAR",method = "REML")
summary(ffs_meta_5_new)

#8 months FFS####
data_es_ffs_amb_8<-data_es_ffs_amb %>% filter(TSD_months==8)
ffs_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_8,struct = "HAR",method = "REML")
summary(ffs_meta_8)

ffs_meta_8_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_ffs_amb_8,struct = "HAR",method = "REML")
summary(ffs_meta_8_new)

#12 months FFS####
data_es_ffs_amb_12<-data_es_ffs_amb %>% filter(TSD_months==12)
ffs_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_12,struct = "HAR",method = "REML")
summary(ffs_meta_12)

ffs_meta_12_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_12,struct = "HAR",method = "REML")
summary(ffs_meta_12_new)
-(1-exp(ffs_meta_12_new$b))*100
((exp(ffs_meta_12_new$se))-1)*100

#15 months FFS####
data_es_ffs_amb_15<-data_es_ffs_amb %>% filter(TSD_months==15)
ffs_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_15,struct = "HAR",method = "REML")
summary(ffs_meta_15)

ffs_meta_15_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_15,struct = "HAR",method = "REML")
summary(ffs_meta_15_new)

#18 months FFS####
data_es_ffs_amb_18<-data_es_ffs_amb %>% filter(TSD_months==18)
ffs_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_18,struct = "HAR",method = "REML")
summary(ffs_meta_18)

ffs_meta_18_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_18,struct = "HAR",method = "REML")
summary(ffs_meta_18_new)

#21 months FFS####
data_es_ffs_amb_21<-data_es_ffs_amb %>% filter(TSD_months==21)
ffs_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_21,struct = "HAR",method = "REML")
summary(ffs_meta_21)

ffs_meta_21_new<- rma.mv(yi_new,vi_new,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_21,struct = "HAR",method = "REML")
summary(ffs_meta_21_new)

#Data frame for figures
ffs_pant_res<-rbind(data.frame(Months="1", estimate=ffs_meta_1$b,se=ffs_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=ffs_meta_3$b,se=ffs_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=ffs_meta_5$b,se=ffs_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=ffs_meta_8$b,se=ffs_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=ffs_meta_12$b,se=ffs_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=ffs_meta_15$b,se=ffs_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=ffs_meta_18$b,se=ffs_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=ffs_meta_21$b,se=ffs_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


ffs_pant_res
#data frame for figure with new resilience metric
ffs_pant_res_new<-rbind(data.frame(Months="1", estimate=ffs_meta_1_new$b,se=ffs_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="3", estimate=ffs_meta_3_new$b,se=ffs_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="5", estimate=ffs_meta_5_new$b,se=ffs_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="8", estimate=ffs_meta_8_new$b,se=ffs_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="12", estimate=ffs_meta_12_new$b,se=ffs_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=ffs_meta_15_new$b,se=ffs_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="18", estimate=ffs_meta_18_new$b,se=ffs_meta_18_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="21", estimate=ffs_meta_21_new$b,se=ffs_meta_21_new$se,row.names=FALSE, stringsAsFactors=TRUE))


ffs_pant_res_new

##Plotting points and Pantropical effect sizes
Fig_res_amb_ffs <- ggplot(data_es_ffs_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_ffs<-Fig_res_amb_ffs+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_ffs<-Fig_res_amb_ffs+geom_pointrange(data=ffs_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_ffs<-Fig_res_amb_ffs+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","black","#665191"))
Fig_res_amb_ffs<-Fig_res_amb_ffs+labs(color="")+guides(size=FALSE)
Fig_res_amb_ffs
Fig_res_amb_ffs<-Fig_res_amb_ffs+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+ylab("Resilience")+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("Time since disturbance (Months)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =6.5, y = 2, fontface="bold",label = "d FFS litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_ffs<-Fig_res_amb_ffs+xlab("Time since cyclone (Months)")
Fig_res_amb_ffs

#New Figure6d FFS####
Fig_res_amb_ffs_new <- ggplot(data_es_ffs_amb, aes(y=yi_new, x=Months,group=Case_study))
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+geom_point(aes(group=Case_study,col=Country,size=vi_new),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+geom_pointrange(data=ffs_pant_res_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","black","#665191"))
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+labs(color="")+guides(size=FALSE)
Fig_res_amb_ffs_new
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-2.5,1)+
  ylab("Resilience")+#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =5.5, y = 1, fontface="bold",label = "d FFS fall", size=8,colour="black")+labs(x="")
Fig_res_amb_ffs_new<-Fig_res_amb_ffs_new+xlab("Time since cyclone (Months)")
Fig_res_amb_ffs_new

#Final Figure6####
Fig6<-Fig_res_amb_tot+Fig_res_amb_wood+Fig_res_amb_leaf+Fig_res_amb_ffs+plot_layout(ncol=2,heights=c(1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig6

Fig6_new<-Fig_res_tot_new+Fig_res_amb_wood_new+Fig_res_amb_leaf_new+Fig_res_amb_ffs_new+plot_layout(ncol=2,heights=c(1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig6_new

#Saving in High Res
ggsave(filename = "Fig6_Resilience-mass_final.png",
       plot = Fig6, width = 23, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

##END###
