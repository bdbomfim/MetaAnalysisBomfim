##Nutrient flux and concentration Resilience##

#Packages####
library(metafor)
library(tidyverse)
library(ggplot2)
library(patchwork)

##Uploading Nutrient data####
nutmeta<-read.csv(file.choose())#Litterfall_Nutrients
str(nutmeta) #2551 obs. of  78 variables
#create Case study column
nutmeta$Case_study= paste(nutmeta$Site, nutmeta$DisturbanceName, sep="| ")
unique(levels(as.factor(nutmeta$Case_study)))

#Nutrient Resilience Plots
nutrec <- nutmeta %>% filter(Cat_TSD_months == "Rec")
str(nutrec)#2430
#changing Fraction names
levels(nutrec$Fraction)
levels(nutrec$Fraction)<-c("FFS","Leaf","Misc.","Total","Wood")

##Nutrient fluxes less than 25 months due to lack of observations for all regions

## P flux
#Filtering the data
res_pf<-nutrec %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>%filter(Fraction!="Misc.")%>%filter(Variable=="P")%>% filter(Raw_Unit=="mg/m2/day")%>% filter (TSD_months < 37)
str(res_pf)#316 obs - this excludes Miscellaneous

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36 Ambient ####
data_esall_amb_Pflux <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = res_pf, measure = "ROM")
summary(data_esall_amb_Pflux$yi)#316 obs. including all litterfall mass fractions

#Creating new columns yi_new as the result of yi divided by duration####
names(data_esall_amb_Pflux)
data_esall_amb_Pflux$yi_new <- data_esall_amb_Pflux$yi / data_esall_amb_Pflux$TSD_months
#checking new column yi_new
summary(data_esall_amb_Pflux$yi_new)
#Same for the variance vi_new
data_esall_amb_Pflux$vi_new <- data_esall_amb_Pflux$vi / data_esall_amb_Pflux$TSD_months
#checking new column vi_new
summary(data_esall_amb_Pflux$vi_new)

##P flux All Fractions####

##Creating new dataframes with grouped data
dres_amb_Pflux <- data_esall_amb_Pflux %>%group_by(Fraction,Case_study,Treatment,TSD_months)  #%>%dplyr::summarise(counts = dplyr::n())
Obs_res_Pflux<- cbind(data.frame(Obs=data_esall_amb_Pflux$yi,Se=sqrt(data_esall_amb_Pflux$vi),Months=factor(data_esall_amb_Pflux$TSD_months),
                               Fraction=data_esall_amb_Pflux$Fraction,SoilP=data_esall_amb_Pflux$Other_soil_P,fac_SoilP=factor(data_esall_amb_Pflux$Other_soil_P),
                               Case_study=data_esall_amb_Pflux$Case_study,Cyclone_freq=data_esall_amb_Pflux$StormFrequencyNorm,Country=data_esall_amb_Pflux$Country))
#new column Months treated as a factor to anable graphics
Obs_res_Pflux$Months<-factor(Obs_res_Pflux$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                              "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))

#New dataframe including yi_new and vi_new
Obs_res_Pflux_new<- cbind(data.frame(Obs=data_esall_amb_Pflux$yi_new,Se=sqrt(data_esall_amb_Pflux$vi_new),Months=factor(data_esall_amb_Pflux$TSD_months),
                                 Fraction=data_esall_amb_Pflux$Fraction,SoilP=data_esall_amb_Pflux$Other_soil_P,fac_SoilP=factor(data_esall_amb_Pflux$Other_soil_P),
                                 Case_study=data_esall_amb_Pflux$Case_study,Cyclone_freq=data_esall_amb_Pflux$StormFrequencyNorm,Country=data_esall_amb_Pflux$Country))
Obs_res_Pflux_new$Months<-factor(Obs_res_Pflux$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                          "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))

##Calculating Pantropical P flux for Total litterfall
tot_Pflux<-data_esall_amb_Pflux %>% filter(Fraction=="Total")
str(tot_Pflux)#92
#If negative, this is the formula: -(1-exp(random1$b))*100
#If positive, this is the formula: ((exp(random0a$b))-1)*100

#1 and 2 months Pflux-TF####
data_es_totPflux_1<-tot_Pflux %>% filter(TSD_months<3)
totPflux_meta_1<- rma.mv(yi,vi,
                    tdist = TRUE,
                    data = data_es_totPflux_1,struct = "HAR",method = "REML")
summary(totPflux_meta_1)
-(1-exp(totPflux_meta_1$b))*100 #76.9% below
((exp(totPflux_meta_1$se))-1)*100

#with new resilience metric
totPflux_meta_1_new<- rma.mv(yi_new,vi_new,
                         tdist = TRUE,
                         data = data_es_totPflux_1,struct = "HAR",method = "REML")
summary(totPflux_meta_1_new)

#3 and 4 months Pflux-TF####
data_es_totPflux_3<-tot_Pflux %>% filter(TSD_months==3|TSD_months==4)
totPflux_meta_3<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_totPflux_3,struct = "HAR",method = "REML")
summary(totPflux_meta_3)

totPflux_meta_3_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totPflux_3,struct = "HAR",method = "REML")
summary(totPflux_meta_3_new)

#5 and 6 months Pflux_TF####
data_es_totPflux_5<-tot_Pflux %>% filter(TSD_months==5|TSD_months==6)
totPflux_meta_5<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_totPflux_5,struct = "HAR",method = "REML")
summary(totPflux_meta_5)

totPflux_meta_5_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totPflux_5,struct = "HAR",method = "REML")
summary(totPflux_meta_5_new)

#8 and 9 months Pflux-TF####
data_es_totPflux_8<-tot_Pflux %>% filter(TSD_months==8|TSD_months==9)
totPflux_meta_8<- rma.mv(yi,vi,#random = ~(1|Region),
                    tdist = TRUE,
                    data = data_es_totPflux_8,struct = "HAR",method = "REML")
summary(totPflux_meta_8)

totPflux_meta_8_new<- rma.mv(yi_new,vi_new,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_totPflux_8,struct = "HAR",method = "REML")
summary(totPflux_meta_8_new)

#11 and 12 months Pflux-TF####
data_es_totPflux_12<-tot_Pflux %>% filter(TSD_months==11|TSD_months==12)
totPflux_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_totPflux_12,struct = "HAR",method = "REML")
summary(totPflux_meta_12)

totPflux_meta_12_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totPflux_12,struct = "HAR",method = "REML")
summary(totPflux_meta_12_new)

#14 and 15 months Pflux-TF####
data_es_totPflux_15<-tot_Pflux %>% filter(TSD_months==14|TSD_months==15)
totPflux_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_totPflux_15,struct = "HAR",method = "REML")
summary(totPflux_meta_15)

totPflux_meta_15_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totPflux_15,struct = "HAR",method = "REML")
summary(totPflux_meta_15_new)

#Data frame to prepare Pflux figure
tot_pant_Pflux<-rbind(data.frame(Months="1", estimate=totPflux_meta_1$b,se=totPflux_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="3", estimate=totPflux_meta_3$b,se=totPflux_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="5", estimate=totPflux_meta_5$b,se=totPflux_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="8", estimate=totPflux_meta_8$b,se=totPflux_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="12", estimate=totPflux_meta_12$b,se=totPflux_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=totPflux_meta_15$b,se=totPflux_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE))

tot_pant_Pflux

#New dataframe for new resilience metric
tot_pant_Pflux_new<-rbind(data.frame(Months="1", estimate=totPflux_meta_1_new$b,se=totPflux_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="3", estimate=totPflux_meta_3_new$b,se=totPflux_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="5", estimate=totPflux_meta_5_new$b,se=totPflux_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="8", estimate=totPflux_meta_8_new$b,se=totPflux_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="12", estimate=totPflux_meta_12_new$b,se=totPflux_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="15", estimate=totPflux_meta_15_new$b,se=totPflux_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_Pflux_new

#Figure 7a####
names(Obs_res_Pflux_new)
Fig_res_amb_Pflux <- ggplot(Obs_res_Pflux, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Pflux<-Fig_res_amb_Pflux+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#1dabe6","#ffa600","black","#665191","#af060f"))##167923
Fig_res_amb_Pflux<-Fig_res_amb_Pflux+geom_pointrange(data=tot_pant_Pflux,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=22)+theme_bw()
Fig_res_amb_Pflux+scale_shape_discrete(solid=F)
Fig_res_amb_Pflux<-Fig_res_amb_Pflux+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 5, 8,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3.5,3)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+scale_shape_discrete(solid=F)+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  labs(x="Time since disturbance (Months)",y="Resilience")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =4, y = 3, fontface="bold",label = "a P flux", size=8,colour="black")+labs(x="")
Fig_res_amb_Pflux<-Fig_res_amb_Pflux+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=4))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Pflux

#New Figure7a####
Fig_res_amb_Pflux_new <- ggplot(Obs_res_Pflux_new, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Pflux_new<-Fig_res_amb_Pflux_new+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#1dabe6","#ffa600","black","#665191","#af060f"))##167923
Fig_res_amb_Pflux_new<-Fig_res_amb_Pflux_new+geom_pointrange(data=tot_pant_Pflux_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=22)+theme_bw()
Fig_res_amb_Pflux_new+scale_shape_discrete(solid=F)
Fig_res_amb_Pflux_new<-Fig_res_amb_Pflux_new+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 5, 8,12, 15, 18,21,24,27,30,33,36))+
  ylim(-2,2)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+scale_shape_discrete(solid=F)+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  labs(x="Time since disturbance (Months)",y="Resilience")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =4, y = 2, fontface="bold",label = "a P flux", size=8,colour="black")+labs(x="")
Fig_res_amb_Pflux_new<-Fig_res_amb_Pflux_new+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=4))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Pflux_new

## N flux####
#Subsetting dataset
res_nf<-nutrec %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>%filter(Fraction!="Misc.")%>%filter(Variable=="N")%>% filter(Raw_Unit=="mg/m2/day")%>% filter (TSD_months < 37)
str(res_nf)#313
res_nf$Case_study= paste(res_nf$Site, res_nf$Disturbance, sep="|")

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36 Ambient ####
data_esall_amb_Nflux <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                               sd1i = Post_SD, sd2i = Pre_SD, data = res_nf, measure = "ROM")
str(data_esall_amb_Nflux)#313 obs. including all litterfall mass fractions

#Creating new columns yi_new as the result of yi divided by duration####
data_esall_amb_Nflux$yi_new <- data_esall_amb_Nflux$yi / data_esall_amb_Nflux$TSD_months
#checking new column yi_new
summary(data_esall_amb_Nflux$yi_new)
#Same for the variance vi_new
data_esall_amb_Nflux$vi_new <- data_esall_amb_Nflux$vi / data_esall_amb_Nflux$TSD_months
#checking new column vi_new
summary(data_esall_amb_Nflux$vi_new)

##P flux All Fractions

##Preparing the data
dres_amb_Nflux <- data_esall_amb_Nflux %>%group_by(Fraction,Case_study,Treatment,TSD_months)  #%>%dplyr::summarise(counts = dplyr::n())
Obs_res_Nflux<- cbind(data.frame(Obs=data_esall_amb_Nflux$yi,Se=sqrt(data_esall_amb_Nflux$vi),Months=factor(data_esall_amb_Nflux$TSD_months),
                                 Fraction=data_esall_amb_Nflux$Fraction,SoilP=data_esall_amb_Nflux$Other_soil_P,fac_SoilP=factor(data_esall_amb_Nflux$Other_soil_P),
                                 Case_study=data_esall_amb_Nflux$Case_study,Cyclone_freq=data_esall_amb_Nflux$StormFrequencyNorm,Country=data_esall_amb_Nflux$Country))
Obs_res_Nflux$Months<-factor(Obs_res_Nflux$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                              "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
#Including new resilience metric and variance
Obs_res_Nflux_new<- cbind(data.frame(Obs=data_esall_amb_Nflux$yi_new,Se=sqrt(data_esall_amb_Nflux$vi_new),Months=factor(data_esall_amb_Nflux$TSD_months),
                                 Fraction=data_esall_amb_Nflux$Fraction,SoilP=data_esall_amb_Nflux$Other_soil_P,fac_SoilP=factor(data_esall_amb_Nflux$Other_soil_P),
                                 Case_study=data_esall_amb_Nflux$Case_study,Cyclone_freq=data_esall_amb_Nflux$StormFrequencyNorm,Country=data_esall_amb_Nflux$Country))
Obs_res_Nflux_new$Months<-factor(Obs_res_Nflux$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                              "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))

##Calculating Pantropical P flux for Total litterfall###
#Total litterfall N flux Nflux-TF####
tot_Nflux<-data_esall_amb_Nflux %>% filter(Fraction=="Total")
#1 and 2 months Nflux-TF####
data_es_totNflux_1<-tot_Nflux %>% filter(TSD_months<3)
totNflux_meta_1<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_1,struct = "HAR",method = "REML")
summary(totNflux_meta_1)
-(1-exp(totNflux_meta_1$b))*100
((exp(totNflux_meta_1$se))-1)*100

#with new resilience metric
totNflux_meta_1_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                        tdist = TRUE,
                        data = data_es_totNflux_1,struct = "HAR",method = "REML")
summary(totNflux_meta_1_new)

#3 and 4 months Nflux-TF####
data_es_totNflux_3<-tot_Nflux %>% filter(TSD_months==3|TSD_months==4)
totNflux_meta_3<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_3,struct = "HAR",method = "REML")
summary(totNflux_meta_3)
-(1-exp(totNflux_meta_3$b))*100
#with new resilience metric
totNflux_meta_3_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_3,struct = "HAR",method = "REML")
summary(totNflux_meta_3_new)

#5 and 6 months Nflux_TF####
data_es_totNflux_5<-tot_Nflux %>% filter(TSD_months==5|TSD_months==6)
totNflux_meta_5<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_5,struct = "HAR",method = "REML")
summary(totNflux_meta_5)
#with new resilience metric
totNflux_meta_5_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_5,struct = "HAR",method = "REML")
summary(totNflux_meta_5_new)

#8 and 9 months Nflux-TF####
data_es_totNflux_8<-tot_Nflux %>% filter(TSD_months==8|TSD_months==9)
totNflux_meta_8<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_8,struct = "HAR",method = "REML")
summary(totNflux_meta_8)
#resilience metric
totNflux_meta_8_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_totNflux_8,struct = "HAR",method = "REML")
summary(totNflux_meta_8_new)

#11 and 12 months Nflux-TF####
data_es_totNflux_12<-tot_Nflux %>% filter(TSD_months==11|TSD_months==12)
totNflux_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totNflux_12,struct = "HAR",method = "REML")
summary(totNflux_meta_12)
#new resilience metric
totNflux_meta_12_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totNflux_12,struct = "HAR",method = "REML")
summary(totNflux_meta_12_new)

#14 and 15 months Nflux-TF####
data_es_totNflux_15<-tot_Nflux %>% filter(TSD_months==14|TSD_months==15)
totNflux_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totNflux_15,struct = "HAR",method = "REML")
summary(totNflux_meta_15)
#resilience metric
totNflux_meta_15_new<- rma.mv(yi_new,vi_new,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_totNflux_15,struct = "HAR",method = "REML")
summary(totNflux_meta_15_new)

#Data frame for figure
tot_pant_Nflux<-rbind(data.frame(Months="1", estimate=totNflux_meta_1$b,se=totNflux_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="3", estimate=totNflux_meta_3$b,se=totNflux_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="5", estimate=totNflux_meta_5$b,se=totNflux_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="8", estimate=totNflux_meta_8$b,se=totNflux_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="12", estimate=totNflux_meta_12$b,se=totNflux_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="15", estimate=totNflux_meta_15$b,se=totNflux_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_Nflux
#New data frame for new resilience metric and variance
tot_pant_Nflux_new<-rbind(data.frame(Months="1", estimate=totNflux_meta_1_new$b,se=totNflux_meta_1_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="3", estimate=totNflux_meta_3_new$b,se=totNflux_meta_3_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="5", estimate=totNflux_meta_5_new$b,se=totNflux_meta_5_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="8", estimate=totNflux_meta_8_new$b,se=totNflux_meta_8_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="12", estimate=totNflux_meta_12_new$b,se=totNflux_meta_12_new$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="15", estimate=totNflux_meta_15_new$b,se=totNflux_meta_15_new$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_Nflux_new

#Figure 7b####
Fig_res_amb_Nflux <- ggplot(Obs_res_Nflux, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Nflux<-Fig_res_amb_Nflux+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#1dabe6","#ffa600","black","#665191","#af060f"))##1C39A8
Fig_res_amb_Nflux<-Fig_res_amb_Nflux+geom_pointrange(data=tot_pant_Nflux,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=22)+theme_bw()
Fig_res_amb_Nflux<-Fig_res_amb_Nflux+theme_bw()+scale_shape_discrete(solid=F)+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 5, 8,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3.5,3)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+ #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  ylab("Resilience")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right",legend.justification="left")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =4, y = 3, fontface="bold",label = "b N flux", size=8,colour="black")+labs(x="")
Fig_res_amb_Nflux<-Fig_res_amb_Nflux+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=4))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Nflux<-Fig_res_amb_Nflux+labs(x="Time since disturbance (Months)")
Fig_res_amb_Nflux

#New Figure7b####
Fig_res_amb_Nflux_new <- ggplot(Obs_res_Nflux_new, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Nflux_new<-Fig_res_amb_Nflux_new+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#1dabe6","#ffa600","black","#665191","#af060f"))##1C39A8
Fig_res_amb_Nflux_new<-Fig_res_amb_Nflux_new+geom_pointrange(data=tot_pant_Nflux_new,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=22)+theme_bw()
Fig_res_amb_Nflux_new<-Fig_res_amb_Nflux_new+theme_bw()+scale_shape_discrete(solid=F)+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 5, 8,12, 15, 18,21,24,27,30,33,36))+
  ylim(-2.5,2)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+ #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  ylab("Resilience")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right",legend.justification="left")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =4, y = 2, fontface="bold",label = "b N flux", size=8,colour="black")+labs(x="")
Fig_res_amb_Nflux_new<-Fig_res_amb_Nflux_new+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=4))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Nflux_new<-Fig_res_amb_Nflux_new+labs(x="Time since disturbance (Months)")
Fig_res_amb_Nflux_new

##P concentration in Leaf + Wood fall ####
#Subsetting dataset
res_pc<-nutrec %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>%filter(Fraction!="Misc.")%>%filter(Variable=="P")%>% filter(Raw_Unit=="mg/g")%>% filter (TSD_months < 37)
str(res_pc)#128
res_pf$Fraction

res_pc$Case_study= paste(res_pc$Site, res_pc$Disturbance, sep="|")

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36 Ambient ####
data_esall_amb_Pc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                               sd1i = Post_SD, sd2i = Pre_SD, data = res_pc, measure = "ROM")
str(data_esall_amb_Pc)#316 obs. including all litterfall mass fractions

data_esall_amb_Pc$Fraction

##Preparing the data
dres_amb_Pc <- data_esall_amb_Pc %>%group_by(Fraction,Case_study,Treatment,TSD_months)  #%>%dplyr::summarise(counts = dplyr::n())
dres_amb_Pc

Obs_res_Pc<- cbind(data.frame(Obs=data_esall_amb_Pc$yi,Se=sqrt(data_esall_amb_Pc$vi),Months=factor(data_esall_amb_Pc$TSD_months),
                                 Fraction=data_esall_amb_Pc$Fraction,SoilP=data_esall_amb_Pc$Other_soil_P,fac_SoilP=factor(data_esall_amb_Pc$Other_soil_P),
                                 Case_study=data_esall_amb_Pc$Case_study,Cyclone_freq=data_esall_amb_Pc$StormFrequencyNorm,Country=data_esall_amb_Pc$Country))
str(Obs_res_Pc)

Obs_res_Pc$Months<-factor(Obs_res_Pc$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                              "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
levels(Obs_res_Pc$Months)

Obs_res_Pc$Fraction

##Calculating Leaf fall P concentration
leaf_Pc<-data_esall_amb_Pc %>% filter(Fraction=="Leaf")
str(leaf_Pc)#67

#1 and 2 months
data_es_lpc_1<-leaf_Pc %>% filter(TSD_months<3)
str(data_es_lpc_1)

lpc_meta_1<- rma.mv(yi,vi,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_lpc_1,struct = "HAR",method = "REML")
summary(lpc_meta_1)
((exp(lpc_meta_1$b))-1)*100

#3 and 4 months
data_es_lpc_3<-leaf_Pc %>% filter(TSD_months==3|TSD_months==4)

lpc_meta_3<- rma.mv(yi,vi,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_lpc_3,struct = "HAR",method = "REML")
summary(lpc_meta_3)
((exp(lpc_meta_3$b))-1)*100
((exp(lpc_meta_3$se))-1)*100

#5 and 6 months
data_es_lpc_5<-leaf_Pc %>% filter(TSD_months==5|TSD_months==6)
str(data_es_lpc_5)

lpc_meta_5<- rma.mv(yi,vi,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_lpc_5,struct = "HAR",method = "REML")
summary(lpc_meta_5)

#8 and 9 months
data_es_lpc_8<-leaf_Pc %>% filter(TSD_months==8|TSD_months==9)
str(data_es_lpc_8)

lpc_meta_8<- rma.mv(yi,vi,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_lpc_8,struct = "HAR",method = "REML")
summary(lpc_meta_8)

#11 and 12 months
data_es_lpc_12<-leaf_Pc %>% filter(TSD_months==11|TSD_months==12)
str(data_es_lpc_12)

lpc_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_lpc_12,struct = "HAR",method = "REML")
summary(lpc_meta_12)

#14 and 15 months
data_es_lpc_15<-leaf_Pc %>% filter(TSD_months==14|TSD_months==15)
str(data_es_lpc_15)

lpc_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_lpc_15,struct = "HAR",method = "REML")
summary(lpc_meta_15)

#Data frame
tot_pant_Pc<-rbind(data.frame(Months="1", estimate=lpc_meta_1$b,se=lpc_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="3", estimate=lpc_meta_3$b,se=lpc_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="5", estimate=lpc_meta_5$b,se=lpc_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="8", estimate=lpc_meta_8$b,se=lpc_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="12", estimate=lpc_meta_12$b,se=lpc_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                      data.frame(Months="15", estimate=lpc_meta_15$b,se=lpc_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_Pc

#Figure
Fig_res_amb_Pc <- ggplot(Obs_res_Pc, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Pc<-Fig_res_amb_Pc+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#ffa600","black","#665191","#af060f"))##167923 "#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"        #"#ffa600","black","#665191","#af060f"
Fig_res_amb_Pc<-Fig_res_amb_Pc+geom_pointrange(data=tot_pant_Pc,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=24)+theme_bw()
Fig_res_amb_Pc<-Fig_res_amb_Pc+theme_bw()+scale_shape_discrete(solid=F)+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-1,1.5)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+ 
  ylab("Resilience")+#ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("")+
  theme(axis.title.x =element_text(vjust = 0.5,size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =7.6, y = 1.5, fontface="bold",label = "c P concentration", size=8,colour="black")+labs(x="")
Fig_res_amb_Pc<-Fig_res_amb_Pc+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=3))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Pc

#patchwork to produce final figure
FigNP<-Fig_res_amb_Pflux+Fig_res_amb_Pc+Fig_res_amb_Nflux+plot_layout(ncol=1,heights=c(1,1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
FigNP

ggsave(filename = "Fig_Res_NPflux.png",
       plot = FigNPflux, width = 20, height = 16, units = 'cm',
       scale = 2, dpi = 600)

##N concentration Leaf ####

#Subsetting dataset
res_nc<-nutrec %>%filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>%filter(Fraction!="Misc.")%>%filter(Variable=="N")%>% filter(Raw_Unit=="mg/g")%>% filter (TSD_months < 37)
str(res_nc)#126 obs

res_nc$Case_study= paste(res_nc$Site, res_nc$Disturbance, sep="|")

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36 Ambient ####
data_esall_amb_Nc <- escalc(n1i = Post_n, n2i = Pre_n, m1i = Post_Mean, m2i = Pre_Mean, 
                            sd1i = Post_SD, sd2i = Pre_SD, data = res_nc, measure = "ROM")
str(data_esall_amb_Nc)#316 obs. including all litterfall mass fractions

data_esall_amb_Nc$Fraction

##Preparing the data
dres_amb_Nc <- data_esall_amb_Nc %>%group_by(Fraction,Case_study,Treatment,TSD_months)  #%>%dplyr::summarise(counts = dplyr::n())
dres_amb_Nc

Obs_res_Nc<- cbind(data.frame(Obs=data_esall_amb_Nc$yi,Se=sqrt(data_esall_amb_Nc$vi),Months=factor(data_esall_amb_Nc$TSD_months),
                              Fraction=data_esall_amb_Nc$Fraction,SoilP=data_esall_amb_Nc$Other_soil_P,fac_SoilP=factor(data_esall_amb_Nc$Other_soil_P),
                              Case_study=data_esall_amb_Nc$Case_study,Cyclone_freq=data_esall_amb_Nc$StormFrequencyNorm,Country=data_esall_amb_Nc$Country))
str(Obs_res_Nc)

Obs_res_Nc$Months<-factor(Obs_res_Nc$Months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                        "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
levels(Obs_res_Nc$Months)

##Calculating Pantropical Leaf N concentration
leaf_Nc<-data_esall_amb_Nc %>% filter(Fraction=="Leaf")
str(leaf_Nc)#67

#1 and 2 months
data_es_lnc_1<-leaf_Nc %>% filter(TSD_months<3)

lnc_meta_1<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_lnc_1,struct = "HAR",method = "REML")
summary(lnc_meta_1)
((exp(lnc_meta_1$b))-1)*100
summary(lpc_meta_1)

#3 and 4 months
data_es_lnc_3<-leaf_Nc %>% filter(TSD_months==3|TSD_months==4)
str(data_es_lnc_3)

lnc_meta_3<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_lnc_3,struct = "HAR",method = "REML")
summary(lnc_meta_3)
((exp(lnc_meta_3$b))-1)*100
((exp(lnc_meta_3$se))-1)*100

#5 and 6 months
data_es_lnc_5<-leaf_Nc %>% filter(TSD_months==5|TSD_months==6)
str(data_es_lnc_5)

lnc_meta_5<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_lnc_5,struct = "HAR",method = "REML")
summary(lnc_meta_5)
((exp(lnc_meta_5$b))-1)*100
((exp(lnc_meta_5$se))-1)*100

#8 and 9 months
data_es_lnc_8<-leaf_Nc %>% filter(TSD_months==8|TSD_months==9)
str(data_es_lnc_8)

lnc_meta_8<- rma.mv(yi,vi,#random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_lnc_8,struct = "HAR",method = "REML")
summary(lnc_meta_8)

#11 and 12 months
data_es_lnc_12<-leaf_Nc %>% filter(TSD_months==11|TSD_months==12)
str(data_es_lnc_12)

lnc_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_lnc_12,struct = "HAR",method = "REML")
summary(lnc_meta_12)

#14 and 15 months
data_es_lnc_15<-leaf_Nc %>% filter(TSD_months==14|TSD_months==15)
str(data_es_lnc_15)

lnc_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_lnc_15,struct = "HAR",method = "REML")
summary(lnc_meta_15)

#Data frame
tot_pant_Nc<-rbind(data.frame(Months="1", estimate=lnc_meta_1$b,se=lnc_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                   data.frame(Months="3", estimate=lnc_meta_3$b,se=lnc_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                   data.frame(Months="5", estimate=lnc_meta_5$b,se=lnc_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                   data.frame(Months="8", estimate=lnc_meta_8$b,se=lnc_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                   data.frame(Months="12", estimate=lnc_meta_12$b,se=lnc_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                   data.frame(Months="15", estimate=lnc_meta_15$b,se=lnc_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_Nc

#Figure
Fig_res_amb_Nc <- ggplot(Obs_res_Nc, aes(y=Obs, x=Months,group=Case_study))
Fig_res_amb_Nc<-Fig_res_amb_Nc+geom_point(aes(group=Case_study,col=Country,size=Se,shape=Fraction),stroke=1.5,alpha=0.5)+scale_color_manual(values=c("#ffa600","black","#665191","#af060f"))##1C39A8 "#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"        #"#ffa600","black","#665191","#af060f"
Fig_res_amb_Nc<-Fig_res_amb_Nc+geom_pointrange(data=tot_pant_Nc,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1.4, stroke=1.3,shape=24)+theme_bw()
Fig_res_amb_Nc<-Fig_res_amb_Nc+theme_bw()+scale_shape_discrete(solid=F)+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-1,1.5)+ scale_fill_discrete(breaks=c("Total","Leaf","Wood","FFS"))+
  #ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(size=28),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=28),#strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text = element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 5)))+
  annotate("text", x =7.6, y = 1.5, fontface="bold",label = "d N concentration", size=8,colour="black")+labs(x="")
Fig_res_amb_Nc<-Fig_res_amb_Nc+
  guides(colour = guide_legend(nrow = 5, byrow = T, override.aes=list(size=3))) +
  guides(shape = guide_legend(nrow = 4, byrow = T, override.aes=list(size=4)),size=FALSE) +
  theme(legend.direction = "vertical", legend.box = "vertical")+labs(color="",shape="")
Fig_res_amb_Nc<-Fig_res_amb_Nc+labs(x="Time since disturbance (Months)",y="Resilience")
Fig_res_amb_Nc

#Figure 7a-d Nutrient Resilience##
FigNP_v2<-Fig_res_amb_Pflux+Fig_res_amb_Pc+Fig_res_amb_Nflux+Fig_res_amb_Nc+plot_layout(ncol=1,heights=c(1,1,1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
FigNP_v2

#Saving in high res
ggsave(filename = "Fig7_Resilience_NP.png",
       plot = FigNP_v2, width = 20, height = 23, units = 'cm',
       scale = 2, dpi = 1200)

#Final Figure 7a-d Nutrient Resilience####
FigNP_v8<-Fig_res_amb_Pflux+Fig_res_amb_Pc+Fig_res_amb_Nflux+Fig_res_amb_Nc+plot_layout(ncol=2,heights=c(1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
FigNP_v8

ggsave(filename = "Fig7_Resilience_NP.png",
       plot = FigNP_v8, width = 24, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

##END###
