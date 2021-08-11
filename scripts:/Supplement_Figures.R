##Discussion Figures####
library(metafor)
library(ggplot2)
library(patchwork)
library(tidyverse)

##Upload data####
metadat<-read.csv(file.choose())#20210520_Litterfall_Mass
attach(metadat)
str(metadat)#2370 obs of 77 variables
names(metadat)
metadat$Case_study= paste(metadat$Site, metadat$DisturbanceName,metadat$Treatment, sep=" | ")

####FigureS1 Response annual vs subannual####
pfrac_S1<-ggplot(data_frac_S1, aes(x=group,y=-1*estimate,ymax=-1*ci_up,ymin=-1*ci_low, shape = variable))+scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2))+
  geom_pointrange(mapping=aes(color=variable),size=1, position=position_dodge(width=c(0.8, 1.2)))+coord_flip()+
  geom_hline(aes(yintercept=0), lty=2,size=1,col="magenta",alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="Pantropical resistance", x="Litterfall fraction") +scale_shape_discrete(solid=F)+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),axis.title=element_text(size=28),
        axis.text=element_text(size=26),legend.text =  element_text(size=20),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.position = c(0.3,0.89),legend.box.background = element_rect(colour = "gray"))
Fig_S1<-pfrac_S1+ scale_color_manual(values=c("#141212","#5E3FBA","#A81C38"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#scale_color_grey(start=0.65, end=0.25)
Fig_S1<-Fig_S1+ annotate("text", y = -6.5, x = 5.2, fontface="bold",label = "(7)", size=5,colour="black")+ annotate("text", y = -6.5, x = 5, fontface="bold",label = "(7)", size=5,colour="black")+ annotate("text", y = -6.5, x = 4.7, fontface="bold",label = "(11)", size=5,colour="black")+
  annotate("text", y = -6.5, x = 4.22, fontface="bold",label = "(4)", size=5,colour="black")+ annotate("text", y = -6.5, x = 4, fontface="bold",label = "(4)", size=5,colour="black")+ annotate("text", y = -6.5, x = 3.7, fontface="bold",label = "(9)", size=5,colour="black")+
  annotate("text", y = -6.5, x = 3.22, fontface="bold",label = "(16)", size=5,colour="black")+ annotate("text", y = -6.5, x = 3, fontface="bold",label = "(16)", size=5,colour="black")+ annotate("text", y = -6.5, x = 2.7, fontface="bold",label = "(30)", size=5,colour="black")+
  annotate("text", y = -6.5, x = 2.22, fontface="bold",label = "(14)", size=5,colour="black")+ annotate("text", y = -6.5, x = 2, fontface="bold",label = "(14)", size=5,colour="black")+ annotate("text", y = -6.5, x = 1.7, fontface="bold",label = "(29)", size=5,colour="black")+
  annotate("text", y = -6.5, x = 1.22, fontface="bold",label = "(23)", size=5,colour="black")+ annotate("text", y = -6.5, x = 1, fontface="bold",label = "(23)", size=5,colour="black")+ annotate("text", y = -6.5, x = 0.7, fontface="bold",label = "(48)", size=5,colour="black")
Fig_S1
#Saving Figure S1
ggsave(filename = "FigS1_Fractions_Mass_Resistance.png",
       plot = Fig_S1, width = 12, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

#FigureS2 Resilience annual and subannual####

##Data wrangling####
#subset including observations 1 month and on post-cyclone#
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
rec

#subset including 1 to 36 months post-cyclone and excluding duplicated studies and obs. in Bisley
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)

res_allS<-res_all %>% filter(Pre_Mean_MonthSpecific!="NA")#subannual data only

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36####
#Annual all data
data_esall <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_all, measure = "ROM")
str(data_esall)#1527

#Subannual reduced data
data_esallS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean_MonthSpecific, 
                      sd1i = Post_SD, sd2i = Pre_SD_MonthSpecific, data = res_allS, measure = "ROM")
str(data_esallS)#1316

#Annual reduced data
data_esallSS <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                       sd1i = Post_SD, sd2i = Pre_SD, data = res_allS, measure = "ROM")
str(data_esallSS)#1316

#Data frames
levels(data_esall$Treatment)
#Annual all data
tot_lit_all<-data_esall %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_all)#213

#Subannual reduced data
tot_lit_allS<-data_esallS %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_allS)#168

#Annual reduced data
tot_lit_allSS<-data_esallSS %>% filter(Fraction=="TotLitfall")%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>% filter(TSD_months<22)%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_allSS)#168

#Wrangling number of observations
pobs_tot_lit_all<- tot_lit_all %>% group_by(TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(prop = round(counts*100/sum(counts), 1),
                                                  lab.ypos = cumsum(prop) - 0.5*prop,Group="Annual (all)")
pobs_tot_lit_all

pobs_tot_lit_allS<- tot_lit_allS %>% group_by(TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(prop = round(counts*100/sum(counts), 1),
                                                  lab.ypos = cumsum(prop) - 0.5*prop,Group="Subnnual (red.)")
pobs_tot_lit_allS

pobs_tot_lit_allSS<- tot_lit_allSS %>% group_by(TSD_months) %>%
  dplyr::summarise(counts = dplyr::n())%>% mutate(prop = round(counts*100/sum(counts), 1),
                                                  lab.ypos = cumsum(prop) - 0.5*prop,Group="Annual (red.)")
pobs_tot_lit_allSS

#Data frame
Res_all_k<-rbind(data.frame(Group="Annual (all)", TSD=pobs_tot_lit_all$TSD_months,n=pobs_tot_lit_all$counts),
                 data.frame(Group="Subannual (red.)", TSD=pobs_tot_lit_allS$TSD_months,n=pobs_tot_lit_allS$counts),
                 data.frame(Group="Annual (red.)", TSD=pobs_tot_lit_allSS$TSD_months,n=pobs_tot_lit_allSS$counts))
str(Res_all_k)

#Figure S2a####
FigS2a<-ggplot(Res_all_k, aes(x=factor(TSD), y=n, group=Group,color=Group,shape=Group)) +
  geom_linerange(aes(x = factor(TSD), ymin = 0, ymax = n), alpha=0.2, 
                 size = 1)+#geom_hline(aes(yintercept=20), lty=4,size=0.8,col="#C8B9B9")+
  geom_point(size = 4,stroke=1.2)+labs(x = "", y = "Number of observations")+scale_shape_discrete(solid=F)+
  theme_pubr()+theme(axis.title=element_text(size=24),axis.text.y=element_text(size=22),axis.text.x = element_text(angle=0, hjust=0.5,size=22))+
  theme(legend.title = element_blank (), legend.position="top",legend.background = element_rect(fill=alpha('transparent', 0.4)),
        legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.text = element_text(size=26))+ scale_color_manual(values=c("#141212","#5E3FBA","#A81C38"))+scale_y_continuous(breaks=c(0,5,10,15,20))+ #scale_x_discrete(breaks = c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36))+ 
  annotate("text", x = 0.65, y = 15, fontface="bold",label = "a", size=8,colour="black")
FigS2a

#Data wrangling for FigS2b - Pantropical resilience annual and subannual over time since disturbance

#Data frame is ready - from Res_plots_New
tot_pant_res<-rbind(data.frame(Group="Annual (all)",Months="1", estimate=tot_meta_1$b,se=tot_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="3", estimate=tot_meta_3$b,se=tot_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="5", estimate=tot_meta_5$b,se=tot_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="8", estimate=tot_meta_8$b,se=tot_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="12", estimate=tot_meta_12$b,se=tot_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="15", estimate=tot_meta_15$b,se=tot_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="18", estimate=tot_meta_18$b,se=tot_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Annual (all)",Months="21", estimate=tot_meta_21$b,se=tot_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


tot_pant_res

#Subannual reduced

#1 month
data_es_tot_amb_1S<-tot_lit_allS %>% filter(TSD_months==1)
str(data_es_tot_amb_1S)#8 obs
tot_meta_1S<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_1S,struct = "HAR",method = "REML")
summary(tot_meta_1S)
coef(summary(tot_meta_1S))

#3 months
data_es_tot_amb_3S<-tot_lit_allS %>% filter(TSD_months==3)
str(data_es_tot_amb_3S)
tot_meta_3S<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_3S,struct = "HAR",method = "REML")
summary(tot_meta_3S)
coef(summary(tot_meta_3S))

#5 months
data_es_tot_amb_5S<-tot_lit_allS %>% filter(TSD_months==5)
str(data_es_tot_amb_5S)
tot_meta_5S<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_5S,struct = "HAR",method = "REML")
#mods=~TSD_months)
summary(tot_meta_5S)

#8 months
data_es_tot_amb_8S<-tot_lit_allS %>% filter(TSD_months==8)
str(data_es_tot_amb_8S)

tot_meta_8S<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_8S,struct = "HAR",method = "REML")
summary(tot_meta_8S)
coef(summary(tot_meta_8S))

#12 months
data_es_tot_amb_12S<-tot_lit_allS %>% filter(TSD_months==12)
str(data_es_tot_amb_12S)
tot_meta_12S<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_12S,struct = "HAR",method = "REML")
summary(tot_meta_12S)
coef(summary(tot_meta_12S))

#15 months
data_es_tot_amb_15S<-tot_lit_allS %>% filter(TSD_months==15)
str(data_es_tot_amb_15S)
tot_meta_15S<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_15S,struct = "HAR",method = "REML")
summary(tot_meta_15S)
coef(summary(tot_meta_15S))

#18 months
data_es_tot_amb_18S<-tot_lit_allS %>% filter(TSD_months==18)
str(data_es_tot_amb_18S)
tot_meta_18S<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_18S,struct = "HAR",method = "REML")
summary(tot_meta_18S)

#21 months
data_es_tot_amb_21S<-tot_lit_allS %>% filter(TSD_months==21)
str(data_es_tot_amb_21S)
tot_meta_21S<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_21S,struct = "HAR",method = "REML")
summary(tot_meta_21S)

#Data frame subannual red.
tot_pant_resS<-rbind(data.frame(Group="Subannual (red.)",Months="1", estimate=tot_meta_1S$b,se=tot_meta_1S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="3", estimate=tot_meta_3S$b,se=tot_meta_3S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="5", estimate=tot_meta_5S$b,se=tot_meta_5S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="8", estimate=tot_meta_8S$b,se=tot_meta_8S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="12", estimate=tot_meta_12S$b,se=tot_meta_12S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="15", estimate=tot_meta_15S$b,se=tot_meta_15S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="18", estimate=tot_meta_18S$b,se=tot_meta_18S$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Group="Subannual (red.)",Months="21", estimate=tot_meta_21S$b,se=tot_meta_21S$se,row.names=FALSE, stringsAsFactors=TRUE))

tot_pant_resS

#Annual reduced

#1 month
data_es_tot_amb_1SS<-tot_lit_allSS %>% filter(TSD_months==1)
str(data_es_tot_amb_1SS)#8 obs
tot_meta_1SS<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_1SS,struct = "HAR",method = "REML")
summary(tot_meta_1SS)
coef(summary(tot_meta_1SS))

#3 months
data_es_tot_amb_3SS<-tot_lit_allSS %>% filter(TSD_months==3)
str(data_es_tot_amb_3SS)
tot_meta_3SS<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_3SS,struct = "HAR",method = "REML")
summary(tot_meta_3SS)
coef(summary(tot_meta_3SS))

#5 months
data_es_tot_amb_5SS<-tot_lit_allSS %>% filter(TSD_months==5)
str(data_es_tot_amb_5SS)
tot_meta_5SS<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_5SS,struct = "HAR",method = "REML")
summary(tot_meta_5SS)

#8 months
data_es_tot_amb_8SS<-tot_lit_allSS %>% filter(TSD_months==8)
str(data_es_tot_amb_8SS)
tot_meta_8SS<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_8SS,struct = "HAR",method = "REML")
summary(tot_meta_8SS)
coef(summary(tot_meta_8SS))

#12 months
data_es_tot_amb_12SS<-tot_lit_allSS %>% filter(TSD_months==12)
str(data_es_tot_amb_12SS)
tot_meta_12SS<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_tot_amb_12SS,struct = "HAR",method = "REML")
summary(tot_meta_12SS)

#15 months
data_es_tot_amb_15SS<-tot_lit_allSS %>% filter(TSD_months==15)
str(data_es_tot_amb_15SS)
tot_meta_15SS<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_tot_amb_15SS,struct = "HAR",method = "REML")
summary(tot_meta_15SS)
coef(summary(tot_meta_15S))

#18 months
data_es_tot_amb_18SS<-tot_lit_allSS %>% filter(TSD_months==18)
str(data_es_tot_amb_18SS)
tot_meta_18SS<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_tot_amb_18SS,struct = "HAR",method = "REML")
summary(tot_meta_18SS)

#21 months
data_es_tot_amb_21SS<-tot_lit_allSS %>% filter(TSD_months==21)
str(data_es_tot_amb_21SS)
tot_meta_21SS<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_tot_amb_21SS,struct = "HAR",method = "REML")
summary(tot_meta_21SS)

#Data frame subannual red.
tot_pant_resSS<-rbind(data.frame(Group="Annual (red.)",Months="1", estimate=tot_meta_1SS$b,se=tot_meta_1SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="3", estimate=tot_meta_3SS$b,se=tot_meta_3SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="5", estimate=tot_meta_5SS$b,se=tot_meta_5SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="8", estimate=tot_meta_8SS$b,se=tot_meta_8SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="12", estimate=tot_meta_12SS$b,se=tot_meta_12SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="15", estimate=tot_meta_15SS$b,se=tot_meta_15SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="18", estimate=tot_meta_18SS$b,se=tot_meta_18SS$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Group="Annual (red.)",Months="21", estimate=tot_meta_21SS$b,se=tot_meta_21SS$se,row.names=FALSE, stringsAsFactors=TRUE))
tot_pant_resSS

#Combining data frames

FigS2b_data <- rbind(tot_pant_res,tot_pant_resS,tot_pant_resSS)
FigS2b_data

#FigureS2b####

FigS2b<-ggplot(FigS2b_data, aes(x=factor(Months),y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col=Group,group=Group))+scale_y_continuous(breaks=c(-1,-0.5,0,0.5))+
  geom_pointrange(mapping=aes(shape=Group),size=1, stroke=1.5,position=position_dodge(width=c(0.5,0.5,0.5)))+
  #geom_hline(aes(yintercept=0), lty=2,size=1.2,col="magenta") + # this adds a dotted line for effect size of 0
  labs(y="", x="") +scale_shape_discrete(solid=F)+ scale_color_manual(values=c("#141212","#5E3FBA","#A81C38"))+
  theme_pubr()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1.2, alpha = 0.8)+
  theme(axis.title.x =element_text(size=24),
        axis.title.y =element_text(vjust = 0.5,size=24),
        axis.text=element_text(size=22),legend.box="horizontal",legend.text =  element_text(size=22),legend.title = element_blank(),
        legend.position = "none",legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+labs(x="Time since disturbance (Months)", y="Pantropical resilience")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))+
annotate("text", x = 0.5, y = 0.7, fontface="bold",label = "b", size=8,colour="black")
FigS2b

#Final FigureS2####
FigS2<-FigS2a+FigS2b+plot_layout(ncol=1,heights=c(0.8,1))
FigS2

#Saving FigureS2
ggsave(filename = "FigS2_Resilience_annual_subannual.png",
       plot = FigS2, width = 24, height = 16, units = 'cm',
       scale = 2, dpi = 1000)

##Figure Total Litterfall Mass Flux Response Predictors variability#####

names(CPFig4b_cor)#these are the variables to be plotted

sup_fig_mod<-data_es0ia%>% filter(HURRECON_wind_ms!="NA")
names(sup_fig_mod)

#Boxplot showing range in data points for the reponse predictors of total litterfall mass flux
#Holdridge
df_h <- sup_fig_mod %>%
  group_by(Holdridge_life_zone,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_h<-df_h %>% ggplot(aes(x=Holdridge_life_zone, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=60, hjust=1,vjust = 1,size=12))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+
  labs(y = "n",x = "Holdridge life zone")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

df_rt <- sup_fig_mod %>%
  group_by(RockType,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_rt<-df_rt %>% ggplot(aes(x=RockType, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=45, hjust=1,vjust = 1,size=12))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(1,3,5,7,9))+
  labs(y = "n",x = "Geological group")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

df_rp <- sup_fig_mod %>%
  group_by(RockPClass,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_rp<-df_rp %>% ggplot(aes(x=RockPClass, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 18), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(1,3,5,7,9,11,13))+
  labs(y = "n",x = "Parent material P class")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

df_pm <- sup_fig_mod %>%
  group_by(Par_Mat,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_pm<-df_pm %>% ggplot(aes(x=Par_Mat, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=60, hjust=1,vjust = 1,size=12))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 18), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(1,3,5,7,9,11,13))+
  labs(y = "n",x = "Parent material")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

df_so <- sup_fig_mod %>%
  group_by(USDA_Soil_Order,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_so<-df_so %>% ggplot(aes(x=USDA_Soil_Order, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=40, hjust=1,vjust = 1,size=12))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(1,3,5,7,9,11,13))+
  labs(y = "n",x = "Soil order (US)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))

#Longitude
df_lo <- sup_fig_mod %>%
  group_by(Longitude,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_lo<-df_lo %>% ggplot(aes(x=Longitude, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_y_continuous(breaks=c(1,2,3,4,5))+scale_x_continuous(breaks=c(-150,-100,-50,0,50,100,150))+
  labs(y = "n",x = "Longitude (UTM)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_lo <- df_lo %>%
  ggplot(aes(x=Longitude)) +
  geom_density(fill="#e9ecef",alpha=0.8,position = 'identity') +
  #scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(fill="")+labs(y = "Density",x = "Longitude (UTM)")+
  geom_label( aes(x=80, y=0.015, label="Longitude"), color="#404080",size=5.5)
p_lo

#Elevation
df_el <- sup_fig_mod %>%
  group_by(Elevation_m,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_el<-df_el %>% ggplot(aes(x=Elevation_m, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_x_continuous(breaks=c(0,100,300,500,700,900,1100,1300))+
  labs(y = "n",x = "Elevation (m)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_el <- df_el %>%
  ggplot(aes(x=Elevation_m)) +
  geom_density(fill="#e9ecef",alpha=0.8,position = 'identity') +
  #scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(fill="")+labs(y = "Density",x = "Elevation (m)")+
  geom_label( aes(y=0.001, x=880, label="Elevation"), color="#404080",size=5.5)
p_el

#MAT/MAP
df_map <- sup_fig_mod %>%
  group_by(MAT_MAP_x100,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_map<-df_map %>% ggplot(aes(x=MAT_MAP_x100, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_x_continuous(breaks=c(0.5,1,1.5,2,2.5,3,3.5))+
  labs(y = "n",x = "MAT/MAP x 100 (C/mm)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_map <- df_map %>%
  ggplot(aes(x=MAT_MAP_x100)) +
  geom_density(fill="#e9ecef",alpha=0.8,position = 'identity') +
  #scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(fill="")+labs(y = "Density",x = "MAT/MAP x 100 (C/mm)")+
  geom_label( aes(x=2.5, y=0.65, label="(MAT/MAP)x100"), color="#404080",size=5.5)
p_map

#Soil P
df_sp <- sup_fig_mod %>%
  group_by(Other_soil_P,Country) %>%
  dplyr::summarise(counts = dplyr::n())
summary(df_sp)

p_df_sp<-df_sp %>% ggplot(aes(x=Other_soil_P, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_x_continuous(breaks=c(0,100,300,500,700,900,1100,1300,1500))+
  labs(y = "n",x = "Total soil P (mg/kg)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_sp <- df_sp %>%
  ggplot(aes(x=log(Other_soil_P))) +
  geom_density(fill="#e9ecef",alpha=0.8,position = 'identity') +
  #scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(fill="")+labs(y = "Density",x = "Total soil P (ln mg/kg)")+
  geom_label( aes(x=3.8, y=0.5, label="Total soil P"), color="#404080",size=5.5)
p_sp

df_sf <- sup_fig_mod %>%
  group_by(StormFrequencyNorm,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_sf<-df_sf %>% ggplot(aes(x=StormFrequencyNorm, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_x_continuous(breaks=c(0.25,0.5,0.75,1))+
  labs(y = "n",x = "Storm frequency (storms/year)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_sf<-df_sf %>% ggplot(aes(x=StormFrequencyNorm)) +
  geom_density(fill="#e9ecef",alpha=0.8) +
  #scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(x="Storm frequency (storms/year)",y="Density")+
  geom_label( aes(x=0.7, y=2.5, label="Storm frequency"), color="#404080",size=5.5)
p_sf

#Time since disturbance
df_tss <- sup_fig_mod %>%
  group_by(YearsSinceLastStorm,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_tss<-df_tss %>% ggplot(aes(x=YearsSinceLastStorm, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+#scale_x_continuous(breaks=c(0.25,0.5,0.75,1))+
  labs(y = "n",x = "Time since last storm (years)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_tss<-df_tss %>% ggplot(aes(x=YearsSinceLastStorm)) +
  geom_density(fill="#e9ecef",alpha=0.8) +
  #scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(x="Time since last storm (years)",y="Density")+
  geom_label( aes(x=8, y=0.09, label="Time since last storm"), color="#404080",size=5.5)
p_tss

#Cyclone rainfall
df_cr <- sup_fig_mod %>%
  group_by(Disturb_Rainfall_mm,Country) %>%
  dplyr::summarise(counts = dplyr::n())
p_df_cr<-df_cr %>% ggplot(aes(x=Disturb_Rainfall_mm, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+#scale_x_continuous(breaks=c(0.25,0.5,0.75,1))+
  labs(y = "n",x = "Cyclone rainfall (mm)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_cr<-df_cr %>% ggplot(aes(x=Disturb_Rainfall_mm)) +
  geom_density(fill="#e9ecef",alpha=0.8) +
  #scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(x="Cyclone rainfall (mm)",y="Density")+
  geom_label( aes(x=600, y=0.0025, label="Cyclone rainfall"), color="#404080",size=5.5)
p_cr

#Wind speed
df_wp <- sup_fig_mod %>%
  group_by(HURRECON_wind_ms,Country) %>%
  dplyr::summarise(counts = dplyr::n())
summary(df_wp)
p_df_wp<-df_wp %>% ggplot(aes(x=HURRECON_wind_ms, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 18), legend.position = "top",
        legend.justification = c("center", "top"),legend.box.just = "right",legend.margin = margin(2, 2, 2, 2))+#scale_x_continuous(breaks=c(0.25,0.5,0.75,1))+
  labs(y = "n",x = "Wind speed (m/s)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_df_wp
p_wp <- df_wp %>%
  ggplot( aes(x=HURRECON_wind_ms)) +
  geom_density(fill="#e9ecef",alpha=0.8) +
  #scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(x="Wind speed (m/s)",y="Density")+
  geom_label( aes(x=23, y=0.025, label="Wind speed"), color="#404080",size=5.5)
p_wp

#Wind duration
df_wd <- sup_fig_mod %>%
  group_by(Gale_wind_duration_minutes,Country) %>%
  dplyr::summarise(counts = dplyr::n())
df_wd
p_df_wd<-df_wd %>% ggplot(aes(x=Gale_wind_duration_minutes, y = counts)) +geom_point(aes(color=Country),shape=21,size=2,stroke=1.5) +theme_classic2()+
  theme(axis.title=element_text(size=18),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=0, hjust=0.5,vjust = 1,size=18))+#guides(fill = guide_legend(title = ""),legend.key.width=32)+ #scale_x_continuous(name="Months since disturbance", breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84))+
  theme(legend.title = element_blank (), legend.text = element_text (size = 24), legend.position = "none",
        legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))+scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))+
  labs(y = "n",x = "Gale wind duration (minutes)")+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
p_wd <- df_wd %>%
  ggplot( aes(x=Gale_wind_duration_minutes)) +
  geom_density(fill="#e9ecef",alpha=0.8) +
  #scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f")) +
  theme_classic2() +labs(x="Wind duration (minutes)",y="Density")+
  geom_label( aes(x=4800, y=0.0002, label="Wind duration"), color="#404080",size=5.5)
p_wd


#sup_fig<-ggplot(CPFig4d_cor, aes(x=Holdridge zone,y=H, fill = meso_cap2.Season)) +
sup_fig_dd_1<-(p_wp|p_lo|p_sf|p_wd)/(p_df_tss|p_df_sp|p_df_map|p_df_cr)/(p_df_el|p_df_so|p_df_rp)/(p_df_rt|p_df_h|p_df_pm)+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 22,face="bold"),plot.tag.position =c(0.02,1))#
sup_fig_dd_1

sup_fig_dd_2<-(p_wp|p_lo|p_sf)/(p_wd|p_tss|p_sp)/(p_map|p_cr|p_el)+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 22,face="bold"),plot.tag.position =c(0.02,1))#
sup_fig_dd_2

#saving as png
ggsave(filename = "Fig_sup_distribution_response.png",
       plot = sup_fig_dd_2, width = 26, height = 26, units = 'cm',scale = 2, dpi = 1000)

##Figure S3 Resilience Correlations####

####Correlation plot of final variables in random forest total litterfall mass flux resilience

#Data
names(datametaforest_restot)
CPFig8b_cor<-datametaforest_restot[,c(2:6,9:18)]
names(CPFig8b_cor)## Final!
#Changing column names
names(CPFig8b_cor)[names(CPFig8b_cor) == "Holdridge_ID"] <- "Holdridge zone"
names(CPFig8b_cor)[names(CPFig8b_cor) == "long"] <- "Longitude"
names(CPFig8b_cor)[names(CPFig8b_cor) == "elev"] <- "Elevation"
names(CPFig8b_cor)[names(CPFig8b_cor) == "mat_map"] <- "MAT/MAP"
names(CPFig8b_cor)[names(CPFig8b_cor) == "soilP"] <- "Soil P"
names(CPFig8b_cor)[names(CPFig8b_cor) == "Soil_ID"] <- "Soil order"
names(CPFig8b_cor)[names(CPFig8b_cor) == "stormfreq"] <- "Storm frequency"
names(CPFig8b_cor)[names(CPFig8b_cor) == "timesincestorm"] <- "Time since last storm"
names(CPFig8b_cor)[names(CPFig8b_cor) == "distrain"] <- "Cyclone rainfall"
names(CPFig8b_cor)[names(CPFig8b_cor) == "hurrwind"] <- "Wind speed"
names(CPFig8b_cor)[names(CPFig8b_cor) == "windur"] <- "Wind duration"
names(CPFig8b_cor)[names(CPFig8b_cor) == "Rocktype_ID"] <- "Geological group"
names(CPFig8b_cor)[names(CPFig8b_cor) == "RockP_ID"] <- "Parent material P"
names(CPFig8b_cor)[names(CPFig8b_cor) == "Par_Mat_ID"] <- "Parent material"
names(CPFig8b_cor)[names(CPFig8b_cor) == "tsd"] <- "Time since cyclone"

#Calculating correlation coefficients and p values
corrf8b <- round(cor(CPFig8b_cor,method="pearson"),2)
p.matf8b <- cor_pmat(CPFig8b_cor)

#FigureS3a correlation
FigS3a<-ggcorrplot(corrf8b, hc.order = TRUE, type = "lower",hc.method = "ward.D2",sig.level = 0.05,
                    outline.col = "white", p.mat = p.matf8b,method="square",ggtheme=ggplot2::theme_classic(),show.legend=TRUE, 
                    legend.title="Pearson's r", lab=TRUE, lab_size=6, tl.cex=28,insig="blank",
                    colors = c("#ABA0A0", "white", "#ffa600",pch.cex=20,nbreaks = 8,legend.text.cex=26))+font("legend.text",size=18)+font("legend.title", size=22)#+theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),axis.text.y = element_text(margin=margin(0,-2,0,0)))
FigS3a

#Saving figure in high res
ggsave(filename = "FigS3a_Resilience_Correlations_TotLit.png",
       plot = FigS3a, width = 16, height = 18, units = 'cm',
       scale = 2, dpi = 800)

#Correlation among leaf fall resilience predictors####

#Selecting columns
names(datametaforest_reslf)
CPFig8d_cor<-datametaforest_reslf[,c(2:6,9:18)]
names(CPFig8d_cor)
#Changing column names
names(CPFig8d_cor)[names(CPFig8d_cor) == "Holdridge_ID"] <- "Holdridge zone"
names(CPFig8d_cor)[names(CPFig8d_cor) == "long"] <- "Longitude"
names(CPFig8d_cor)[names(CPFig8d_cor) == "tsd"] <- "Time since cyclone"
names(CPFig8d_cor)[names(CPFig8d_cor) == "elev"] <- "Elevation"
names(CPFig8d_cor)[names(CPFig8d_cor) == "mat_map"] <- "MAT/MAP"
names(CPFig8d_cor)[names(CPFig8d_cor) == "soilP"] <- "Soil P"
names(CPFig8d_cor)[names(CPFig8d_cor) == "Soil_ID"] <- "Soil order"
names(CPFig8d_cor)[names(CPFig8d_cor) == "stormfreq"] <- "Storm frequency"
names(CPFig8d_cor)[names(CPFig8d_cor) == "timesincestorm"] <- "Time since last storm"
names(CPFig8d_cor)[names(CPFig8d_cor) == "distrain"] <- "Cyclone rainfall"
names(CPFig8d_cor)[names(CPFig8d_cor) == "hurrwind"] <- "Wind speed"
names(CPFig8d_cor)[names(CPFig8d_cor) == "windur"] <- "Wind duration"
names(CPFig8d_cor)[names(CPFig8d_cor) == "Rocktype_ID"] <- "Geological group"
names(CPFig8d_cor)[names(CPFig8d_cor) == "RockP_ID"] <- "Parent material P"
names(CPFig8d_cor)[names(CPFig8d_cor) == "Par_Mat_ID"] <- "Parent material"

#Checking if names are correct
names(CPFig8d_cor)

#Calculating correlation coefficients and p values
corrf8d <- round(cor(CPFig8d_cor,method="pearson"),2)
p.matf8d <- cor_pmat(CPFig8d_cor)

#FigureS3b
FigS3b<-ggcorrplot(corrf8d, hc.order = TRUE, type = "lower",hc.method = "ward.D2",sig.level = 0.05,
                   outline.col = "white", p.mat = p.matf8d,method="square",ggtheme=ggplot2::theme_classic(),show.legend=TRUE, 
                   legend.title="Pearson's r", lab=TRUE, lab_size=6, tl.cex=28,insig="blank",
                   colors = c("#46A332", "white", "#ffa600",pch.cex=20,nbreaks = 8,legend.text.cex=26))+font("legend.text",size=18)+font("legend.title", size=22)#+theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),axis.text.y = element_text(margin=margin(0,-2,0,0)))
FigS3b

#Saving in high res
ggsave(filename = "FigS3b_Resilience_Corr_Leaf.png",
       plot = FigS3b, width = 16, height = 18, units = 'cm',
       scale = 2, dpi = 800)


##Mass flux change with P concentration change####

str(data_es0ia)
data_es0ia$Author
names(data_es0ilf)
str(data_es0ilpc)
str(data_es0ilnc)
data_es0ilpc$Case.ID
head(data_es0ilpc)

data_frac4 <- rbind(data.frame(group="FFS", variable="Annual", estimate=full.model3ff$b, var="Mass flux",
                               ci_low=(full.model3ff$b-(1.96*full.model3ff$se)),ci_up=(full.model3ff$b+(1.96*full.model3ff$se)),
                               row.names=FALSE, stringsAsFactors=TRUE)
dat_nut<-rbind(data.frame(group="Leaf P",Pre_Mean=data_es0ilpc$Pre_Mean, ES=data_es0ilpc$yi,SE=sqrt(data_es0ilpc$vi),Case_study=data_es0ilpc$Case_study, SoilP=data_es0ilpc$Other_soil_P,stringsAsFactors=TRUE),
               data.frame(group="Leaf N",Pre_Mean=data_es0ilnc$Pre_Mean,ES=data_es0ilnc$yi,SE=sqrt(data_es0ilnc$vi),Case_study=data_es0ilnc$Case_study,SoilP=data_es0ilnc$Other_soil_P,stringsAsFactors=TRUE))
dat_nut

##Looking at change in Leaf P concentration by Soil P
pre_nut<-lm(ES~group+log(SoilP), data=dat_nut)
summary(pre_nut)
tab_model(pre_nut)

pre_lnc<-lm(Pre_Mean~log(Other_soil_P), data=data_es0ilnc)
summary(pre_lnc)
tab_model(pre_lnc)

pre_lpc<-lm(Pre_Mean~log(Other_soil_P), data=data_es0ilpc)
summary(pre_lpc)
tab_model(pre_lpc)

##Data viz
#Leaf P
FigS_lpc<-ggplot(data_es0ilpc, aes(x=log(Other_soil_P), y=Pre_Mean))+geom_point(aes(color=log(Other_soil_P)),size=6,alpha=0.9,shape=21,stroke=2)+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=20))+ labs(y = bquote('Leaf fall P'~(mg/g)),x = "Total soil P (ln mg/kg)")
rsq_label.a <- paste('R^2 == 0.54')
FigS_lpc<-FigS_lpc+scale_color_gradient(low="#FCFF00", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,color="Soil P")+ annotate("text", x = 5.4, y = 0.75, label = rsq_label.a, size=6,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS_lpc<-FigS_lpc+ annotate("text", x = 5.4, y = 0.8, label = "c", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS_lnc<-FigS_lnc+ annotate("text", x = 5.4, y = 18, label = "N = 10", size=5,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS_lnc<-FigS_lnc+ annotate("text", x = 3, y = 7.5, label = "a", size=8,hjust=0,fontface="bold",colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS_lpc

#Leaf N
FigS_lnc<-ggplot(data_es0ilnc, aes(x=log(Other_soil_P), y=Pre_Mean))+geom_point(aes(color=log(Other_soil_P)),size=6,alpha=0.9,shape=21,stroke=2)+geom_smooth(method = 'glm', col="#9f8a89", alpha=0.1)+
  theme_pubr()+ # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=20))+ labs(y = bquote('Leaf fall N'~(mg/g)),x = "Total soil P (ln mg/kg)")
rsq_label.b <- paste('R^2 == 0.73')
FigS_lnc<-FigS_lnc+scale_color_gradient(low="#FCFF00", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+#+scale_color_gradient(low="blue", high="red")
  theme(legend.title = element_text (size = 24), legend.text = element_text (size = 22),legend.position = "none",legend.justification = c("right"))+#+scale_size_continuous(range=c(1,20))+
  labs(fill=NULL,color="Soil P")+ annotate("text", x = 5.4, y = 19, label = rsq_label.b, size=6,hjust=0,colour="black",parse=TRUE) #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS_lnc<-FigS_lnc+ annotate("text", x = 5.4, y = 20, label = "d", size=8,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS_lnc<-FigS_lnc+ annotate("text", x = 5.4, y = 17.5, label = "N = 10", size=5,hjust=0,colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
#FigS_lnc<-FigS_lnc+ annotate("text", x = 3, y = 7.5, label = "a", size=8,hjust=0,fontface="bold",colour="black") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))
FigS_lnc

Leaf_nut<-FigS_lpc+FigS_lnc+plot_layout(ncol=2)
Leaf_nut
rlang::last_error()
ggsave(filename = "SupFig_Leaf_Nut_Conc.png",
       plot = Leaf_nut, width = 15, height = 8, units = 'cm',
       scale = 2, dpi = 1000)

FinalFigS3<-FigS3a+FigS3b+ plot_layout(ncol=2)
FinalFigS3

Final_Pre_Sup<-(FigS3a | FigS3b)/(FigS_lpc | FigS_lnc)
Final_Pre_Sup

ggsave(filename = "SupFig_Pre_Mass_Leaf_Nut_SoilP_v3.png",
       plot = Final_Pre_Sup, width = 16, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

#Data wrangliing
unique(levels(as.factor(data_es0ilpc$Case_study_b)))
unique(levels(as.factor(data_es0ia$Case_study2)))

dat_mass_lpc<-data_es0ia %>% filter(Treatment=="Ambient")%>% 
  filter(Case_study2=="Bisley| Hugo| Ambient")%>%filter(Case_study2=="East Peak| Hugo| Ambient")%>%
  filter(Case_study2=="El Verde| Hugo| Ambient")%>%filter(Case_study2=="Kokee| Iniki| Ambient")%>%
  filter(Case_study2=="Lienhuachi| Fungwong| Ambient")%>%filter(Case_study2=="Lienhuachi| Jangmi| Ambient")%>%
  filter(Case_study2=="Lienhuachi| Kalmaegi| Ambient")%>%filter(Case_study2=="Lienhuachi| Sinlaku| Ambient")%>%
  filter(Case_study2=="Bisley| Hugo| Ambient")%>%filter(Case_study2=="East Peak| Hugo| Ambient")%>%
  