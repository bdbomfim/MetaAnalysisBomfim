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
pfrac_S1<-ggplot(data_frac_S1, aes(x=group,y=estimate,ymax=ci_up,ymin=ci_low, shape = variable))+scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8))+
  geom_pointrange(mapping=aes(color=variable),size=1, position=position_dodge(width=c(0.8, 1.2)))+coord_flip()+
  geom_hline(aes(yintercept=0), lty=2,size=1,col="magenta",alpha=0.8) + # this adds a dotted line for effect size of 0
  labs(y="Pantropical response with 95% CI", x="Litterfall fraction") +scale_shape_discrete(solid=F)+
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(axis.title.x =element_text(vjust = -0.5),axis.text.x =element_text(vjust = -0.3),
        axis.title.y =element_text(vjust = 1),axis.title=element_text(size=28),
        axis.text=element_text(size=26),legend.text =  element_text(size=20),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),legend.position = c(0.7,0.89),legend.box.background = element_rect(colour = "gray"))
Fig_S1<-pfrac_S1+ scale_color_manual(values=c("#141212","#5E3FBA","#A81C38"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#scale_color_grey(start=0.65, end=0.25)
Fig_S1<-Fig_S1+ annotate("text", y = 6.5, x = 5.2, fontface="bold",label = "(7)", size=5,colour="black")+ annotate("text", y = 6.5, x = 5, fontface="bold",label = "(7)", size=5,colour="black")+ annotate("text", y = 6.5, x = 4.7, fontface="bold",label = "(11)", size=5,colour="black")+
  annotate("text", y = 6.5, x = 4.22, fontface="bold",label = "(4)", size=5,colour="black")+ annotate("text", y = 6.5, x = 4, fontface="bold",label = "(4)", size=5,colour="black")+ annotate("text", y = 6.5, x = 3.7, fontface="bold",label = "(9)", size=5,colour="black")+
  annotate("text", y = 6.5, x = 3.22, fontface="bold",label = "(16)", size=5,colour="black")+ annotate("text", y = 6.5, x = 3, fontface="bold",label = "(16)", size=5,colour="black")+ annotate("text", y = 6.5, x = 2.7, fontface="bold",label = "(30)", size=5,colour="black")+
  annotate("text", y = 6.5, x = 2.22, fontface="bold",label = "(14)", size=5,colour="black")+ annotate("text", y = 6.5, x = 2, fontface="bold",label = "(14)", size=5,colour="black")+ annotate("text", y = 6.5, x = 1.7, fontface="bold",label = "(29)", size=5,colour="black")+
  annotate("text", y = 6.5, x = 1.22, fontface="bold",label = "(23)", size=5,colour="black")+ annotate("text", y = 6.5, x = 1, fontface="bold",label = "(23)", size=5,colour="black")+ annotate("text", y = 6.5, x = 0.7, fontface="bold",label = "(48)", size=5,colour="black")

#Saving Figure S1
ggsave(filename = "FigS1_Fractions_Mass_Response.png",
       plot = Fig_S1, width = 13, height = 10, units = 'cm',
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
  