##Figure Resilience Observations All Fractions
library(metafor)
library(tidyverse)

####STEP 0 Upload Litterfall mass data####
metadat<-read.csv(file.choose())#20210415_Litterfall_Mass_Flux
attach(metadat)

##Mass
rec <- metadat %>% filter(Cat_TSD_months == "Rec")
rec

#1 to 36 months
res_all<- rec %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% filter (TSD_months < 37)
str(res_all)

res_amb<-res_all %>% filter(Treatment=="Ambient")
str(res_amb)#948

res_amb2<-res_all %>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")
str(res_amb2)#1128

res_amb$Case_study= paste(res_amb$Site, res_amb$DisturbanceName, sep="|")
res_amb2$Case_study= paste(res_amb2$Site, res_amb2$DisturbanceName, sep="|")

#Calculating Effect sizes (i.e. the log transformed ratio of means, or Hedge's g) 1 to 36 Ambient ####
data_esall_amb <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = res_amb, measure = "ROM")
str(data_esall_amb)#948 obs. including all litterfall mass fractions

data_esall_amb2 <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = res_amb2, measure = "ROM")
str(data_esall_amb2)#1128

##Total Litterfall mass resilience [ln(post/pre)] including all Treatments####

#1 to 36 months
tot_amb<-data_esall_amb %>% filter(Fraction=="TotLitfall")
str(tot_amb)#250 obs
levels(tot_amb$Treatment)

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

##Including CTE TrimDeb
data_esall_amb_final2<-data_esall_amb2 %>% filter(Fraction!="Misc fall")
str(data_esall_amb_final2)
unique(levels(as.factor(data_esall_amb_final2$Case_study)))
Se=sqrt(data_esall_amb_final2$vi)
Se
data_esall_amb_final2$yi
Obs=data_esall_amb_final2$yi
Obs
data_esall_amb_final2$Months<-factor(data_esall_amb_final2$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                                  "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
levels(data_esall_amb_final2$Months)


Obs_res_all2<- cbind(data.frame(Obs=data_esall_amb_final2$yi,Se=sqrt(data_esall_amb_final2$vi),Months=factor(data_esall_amb_final2$TSD_months),
                               SoilP=data_esall_amb_final2$Other_soil_P,Fraction=data_esall_amb_final2$Fraction,fac_SoilP=factor(data_esall_amb_final2$Other_soil_P),
                               Case_study=data_esall_amb_final2$Case_study,Cyclone_freq=data_esall_amb_final2$StormFrequencyNorm,Country=data_esall_amb_final2$Country))
str(Obs_res_all2)
data_esall_amb_final2$TSD_months

str(data_esall_amb_final2)
data_esall_amb_final2$Case_study= paste(data_esall_amb_final2$Site, data_esall_amb_final2$DisturbanceName, sep="|")

##Calculating Overall Resilience for certain points in time since disturbance
data_es_tot_amb<-data_esall_amb_final2 %>% filter(Fraction=="TotLitfall")
str(data_es_tot_amb)#286

#1 month
data_es_tot_amb_1<-data_es_tot_amb %>% filter(TSD_months==1)
str(data_es_tot_amb_1)
unique(levels(as.factor(data_es_tot_amb_1$Case_study)))

tot_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_1,struct = "HAR",method = "REML")
summary(tot_meta_1)
coef(summary(tot_meta_1))

#3 months
data_es_tot_amb_3<-data_es_tot_amb %>% filter(TSD_months==3)
str(data_es_tot_amb_3)

tot_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_3,struct = "HAR",method = "REML")
summary(tot_meta_3)
coef(summary(tot_meta_3))

#5 months
data_es_tot_amb_5<-data_es_tot_amb %>% filter(TSD_months==5)
str(data_es_tot_amb_5)

tot_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_5,struct = "HAR",method = "REML")
                    #mods=~TSD_months)
summary(tot_meta_5)
-(1-exp(tot_meta_5$b))*100
coef(summary(tot_meta_5))
summary(tot_meta_5)

#8 months
data_es_tot_amb_8<-data_es_tot_amb %>% filter(TSD_months==8)
str(data_es_tot_amb_8)

tot_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_8,struct = "HAR",method = "REML")
summary(tot_meta_8)
coef(summary(tot_meta_8))

#12 months
data_es_tot_amb_12<-data_es_tot_amb %>% filter(TSD_months==12)
str(data_es_tot_amb_12)

tot_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_tot_amb_12,struct = "HAR",method = "REML")
summary(tot_meta_12)
coef(summary(tot_meta_12))

#15 months
data_es_tot_amb_15<-data_es_tot_amb %>% filter(TSD_months==15)
str(data_es_tot_amb_15)

tot_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_15,struct = "HAR",method = "REML")
summary(tot_meta_15)
coef(summary(tot_meta_15))

#18 months
data_es_tot_amb_18<-data_es_tot_amb %>% filter(TSD_months==18)
str(data_es_tot_amb_18)

tot_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_18,struct = "HAR",method = "REML")
summary(tot_meta_18)

#21 months
data_es_tot_amb_21<-data_es_tot_amb %>% filter(TSD_months==21)
str(data_es_tot_amb_21)

tot_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_tot_amb_21,struct = "HAR",method = "REML")
summary(tot_meta_21)

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

##Plotting points and Pantropical effect sizes
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
  #ylab(expression(Response~(ln~litterfall[ti]/litterfall[t0])))+xlab("Time since disturbance (Months)")+ #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
 theme(axis.title.x =element_text(vjust = 0.5,size=24),axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="vertical",legend.position="right")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.3)))
Fig_res_amb_tot<-Fig_res_amb_tot+annotate("text", x =7, y = 2, fontface="bold",label = "a Total litterfall", size=8,colour="black")+labs(x="")
Fig_res_amb_tot<-Fig_res_amb_tot+ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))
Fig_res_amb_tot

##Leaf fall####

##Calculating Overall Resilience for certain points in time since disturbance

##5 months
data_es_leaf_amb<-data_esall_amb_final2 %>% filter(Fraction=="Leaf fall")
str(data_es_leaf_amb)

leaf_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_leaf_amb,struct = "HAR",method = "REML",mods=~TSD_months)
summary(leaf_meta_all)

mypreds_tot_meta_all<-predict(tot_meta_all,addx=TRUE)
data_es_tot_amb$pred<-mypreds_tot_meta_all$pred
data_es_tot_amb$se<-mypreds_tot_meta_all$se

#1 month
data_es_leaf_amb_1<-data_es_leaf_amb %>% filter(TSD_months==1)
str(data_es_leaf_amb_1)

leaf_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_1,struct = "HAR",method = "REML")
summary(leaf_meta_1)
-(1-exp(leaf_meta_1$b))*100

#3 months
data_es_leaf_amb_3<-data_es_leaf_amb %>% filter(TSD_months==3)
str(data_es_leaf_amb_3)

leaf_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_3,struct = "HAR",method = "REML")
summary(leaf_meta_3)
-(1-exp(leaf_meta_3$b))*100

#5 months
data_es_leaf_amb_5<-data_es_leaf_amb %>% filter(TSD_months==5)
str(data_es_leaf_amb_5)

leaf_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_5,struct = "HAR",method = "REML")
summary(leaf_meta_5)
-(1-exp(leaf_meta_5$b))*100

#8 months
data_es_leaf_amb_8<-data_es_leaf_amb %>% filter(TSD_months==8)
str(data_es_leaf_amb_8)

leaf_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                    tdist = TRUE,
                    data = data_es_leaf_amb_8,struct = "HAR",method = "REML")
summary(leaf_meta_8)

#12 months
data_es_leaf_amb_12<-data_es_leaf_amb %>% filter(TSD_months==12)
str(data_es_leaf_amb_12)

leaf_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_12,struct = "HAR",method = "REML")
summary(leaf_meta_12)

#15 months
data_es_leaf_amb_15<-data_es_leaf_amb %>% filter(TSD_months==15)
str(data_es_leaf_amb_15)

leaf_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_15,struct = "HAR",method = "REML")
summary(leaf_meta_15)

#18 months
data_es_leaf_amb_18<-data_es_leaf_amb %>% filter(TSD_months==18)
str(data_es_leaf_amb_18)

leaf_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_18,struct = "HAR",method = "REML")
summary(leaf_meta_18)

#21 months
data_es_leaf_amb_21<-data_es_leaf_amb %>% filter(TSD_months==21)
str(data_es_leaf_amb_21)

leaf_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_leaf_amb_21,struct = "HAR",method = "REML")
summary(leaf_meta_21)

#data frame
leaf_pant_res<-rbind(data.frame(Months="1", estimate=leaf_meta_1$b,se=leaf_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="3", estimate=leaf_meta_3$b,se=leaf_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="5", estimate=leaf_meta_5$b,se=leaf_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="8", estimate=leaf_meta_8$b,se=leaf_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="12", estimate=leaf_meta_12$b,se=leaf_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="15", estimate=leaf_meta_15$b,se=leaf_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="18", estimate=leaf_meta_18$b,se=leaf_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                    data.frame(Months="21", estimate=leaf_meta_21$b,se=leaf_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


leaf_pant_res

##Plotting points and Pantropical effect sizes
Fig_res_amb_leaf <- ggplot(data_es_leaf_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_leaf<-Fig_res_amb_leaf+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_leaf<-Fig_res_amb_leaf+geom_pointrange(data=leaf_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_leaf<-Fig_res_amb_leaf+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_leaf<-Fig_res_amb_leaf+labs(color="")+guides(size=FALSE)
Fig_res_amb_leaf<-Fig_res_amb_leaf+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
annotate("text", x =5, y = 2, fontface="bold",label = "b Leaf fall", size=8,colour="black")+labs(x="")
Fig_res_amb_leaf<-Fig_res_amb_leaf+xlab("Time since disturbance (Months)")
Fig_res_amb_leaf

##Wood fall####

##Calculating Overall Resilience for certain points in time since disturbance

data_es_wood_amb<-data_esall_amb_final2 %>% filter(Fraction=="Wood fall")
str(data_es_wood_amb)#199

wood_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                       tdist = TRUE,
                       data = data_es_wood_amb,struct = "HAR",method = "REML",mods=~TSD_months)
summary(wood_meta_all)

mypreds_tot_meta_all<-predict(tot_meta_all,addx=TRUE)
data_es_tot_amb$pred<-mypreds_tot_meta_all$pred
data_es_tot_amb$se<-mypreds_tot_meta_all$se

#1 month
data_es_wood_amb_1<-data_es_wood_amb %>% filter(TSD_months==1)
str(data_es_wood_amb_1)

wood_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_1,struct = "HAR",method = "REML")
summary(wood_meta_1)

#3 months
data_es_wood_amb_3<-data_es_wood_amb %>% filter(TSD_months==3)
str(data_es_wood_amb_3)

wood_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_3,struct = "HAR",method = "REML")
summary(wood_meta_3)
-(1-exp(wood_meta_3$b))*100

#5 months
data_es_wood_amb_5<-data_es_wood_amb %>% filter(TSD_months==5)
str(data_es_wood_amb_5)

wood_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_5,struct = "HAR",method = "REML")
summary(wood_meta_5)

#8 months
data_es_wood_amb_8<-data_es_wood_amb %>% filter(TSD_months==8)
str(data_es_wood_amb_8)

wood_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_wood_amb_8,struct = "HAR",method = "REML")
summary(wood_meta_8)

#12 months
data_es_wood_amb_12<-data_es_wood_amb %>% filter(TSD_months==12)
str(data_es_wood_amb_12)

wood_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_12,struct = "HAR",method = "REML")
summary(wood_meta_12)

#15 months
data_es_wood_amb_15<-data_es_wood_amb %>% filter(TSD_months==15)
str(data_es_wood_amb_15)

wood_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_15,struct = "HAR",method = "REML")
summary(wood_meta_15)

#18 months
data_es_wood_amb_18<-data_es_wood_amb %>% filter(TSD_months==18)
str(data_es_wood_amb_18)

wood_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_18,struct = "HAR",method = "REML")
summary(wood_meta_18)

#21 months
data_es_wood_amb_21<-data_es_wood_amb %>% filter(TSD_months==21)
str(data_es_wood_amb_21)

wood_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_wood_amb_21,struct = "HAR",method = "REML")
summary(wood_meta_21)

#data frame
wood_pant_res<-rbind(data.frame(Months="1", estimate=wood_meta_1$b,se=wood_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=wood_meta_3$b,se=wood_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=wood_meta_5$b,se=wood_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=wood_meta_8$b,se=wood_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=wood_meta_12$b,se=wood_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=wood_meta_15$b,se=wood_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=wood_meta_18$b,se=wood_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=wood_meta_21$b,se=wood_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


wood_pant_res

##Plotting points and Pantropical effect sizes

names(data_es_wood_amb)
Fig_res_amb_wood <- ggplot(data_es_wood_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_wood<-Fig_res_amb_wood+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_wood<-Fig_res_amb_wood+geom_pointrange(data=wood_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_wood<-Fig_res_amb_wood+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","#ffa600","black","#665191","#af060f"))
Fig_res_amb_wood<-Fig_res_amb_wood+labs(color="")+guides(size=FALSE)
Fig_res_amb_wood
Fig_res_amb_wood<-Fig_res_amb_wood+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
   scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("Time since disturbance (Months)")+#ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =6, y = 2, fontface="bold",label = "c Wood fall", size=8,colour="black")+labs(x="")
Fig_res_amb_wood

##FFS fall####

##Calculating Overall Resilience for certain points in time since disturbance

data_es_ffs_amb<-data_esall_amb_final2 %>% filter(Fraction=="FFS fall")
str(data_es_ffs_amb)#173

ffs_meta_all<- rma.mv(yi,vi,random = ~(1|Site),
                       tdist = TRUE,
                       data = data_es_ffs_amb,struct = "HAR",method = "REML",mods=~TSD_months)
#mods=~TSD_months)
summary(ffs_meta_all)

mypreds_tot_meta_all<-predict(tot_meta_all,addx=TRUE)
data_es_tot_amb$pred<-mypreds_tot_meta_all$pred
data_es_tot_amb$se<-mypreds_tot_meta_all$se

#1 month
data_es_ffs_amb_1<-data_es_ffs_amb %>% filter(TSD_months==1)
str(data_es_ffs_amb_1)

ffs_meta_1<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_1,struct = "HAR",method = "REML")
summary(ffs_meta_1)

#3 months
data_es_ffs_amb_3<-data_es_ffs_amb %>% filter(TSD_months==3)
str(data_es_ffs_amb_3)

ffs_meta_3<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_3,struct = "HAR",method = "REML")
summary(ffs_meta_3)

#5 months
data_es_ffs_amb_5<-data_es_ffs_amb %>% filter(TSD_months==5)
str(data_es_ffs_amb_5)

ffs_meta_5<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_5,struct = "HAR",method = "REML")
summary(ffs_meta_5)
-(1-exp(ffs_meta_5$b))*100

#8 months
data_es_ffs_amb_8<-data_es_ffs_amb %>% filter(TSD_months==8)
str(data_es_ffs_amb_8)

ffs_meta_8<- rma.mv(yi,vi,random = ~(1|Site),
                     tdist = TRUE,
                     data = data_es_ffs_amb_8,struct = "HAR",method = "REML")
summary(ffs_meta_8)

#12 months
data_es_ffs_amb_12<-data_es_ffs_amb %>% filter(TSD_months==12)
str(data_es_ffs_amb_12)

ffs_meta_12<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_12,struct = "HAR",method = "REML")
summary(ffs_meta_12)

#15 months
data_es_ffs_amb_15<-data_es_ffs_amb %>% filter(TSD_months==15)
str(data_es_ffs_amb_15)

ffs_meta_15<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_15,struct = "HAR",method = "REML")
summary(ffs_meta_15)

#18 months
data_es_ffs_amb_18<-data_es_ffs_amb %>% filter(TSD_months==18)
str(data_es_ffs_amb_18)
ffs_meta_18<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_18,struct = "HAR",method = "REML")
summary(ffs_meta_18)

#21 months
data_es_ffs_amb_21<-data_es_ffs_amb %>% filter(TSD_months==21)
str(data_es_ffs_amb_21)

ffs_meta_21<- rma.mv(yi,vi,random = ~(1|Site),
                      tdist = TRUE,
                      data = data_es_ffs_amb_21,struct = "HAR",method = "REML")
summary(ffs_meta_21)

#data frame
ffs_pant_res<-rbind(data.frame(Months="1", estimate=ffs_meta_1$b,se=ffs_meta_1$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="3", estimate=ffs_meta_3$b,se=ffs_meta_3$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="5", estimate=ffs_meta_5$b,se=ffs_meta_5$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="8", estimate=ffs_meta_8$b,se=ffs_meta_8$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="12", estimate=ffs_meta_12$b,se=ffs_meta_12$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="15", estimate=ffs_meta_15$b,se=ffs_meta_15$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="18", estimate=ffs_meta_18$b,se=ffs_meta_18$se,row.names=FALSE, stringsAsFactors=TRUE),
                     data.frame(Months="21", estimate=ffs_meta_21$b,se=ffs_meta_21$se,row.names=FALSE, stringsAsFactors=TRUE))


ffs_pant_res

##Plotting points and Pantropical effect sizes
Fig_res_amb_ffs <- ggplot(data_es_ffs_amb, aes(y=yi, x=Months,group=Case_study))
Fig_res_amb_ffs<-Fig_res_amb_ffs+geom_point(aes(group=Case_study,col=Country,size=vi),shape=21,stroke=1.3,alpha=0.5)
Fig_res_amb_ffs<-Fig_res_amb_ffs+geom_pointrange(data=ffs_pant_res,mapping=aes(group=Months,x=Months,y=estimate,ymax=estimate+(1.96*se),ymin=estimate-(1.96*se), col="Pantropical"),size=1, stroke=1,shape=21)+theme_bw()
Fig_res_amb_ffs<-Fig_res_amb_ffs+guides(title="Region")+scale_color_manual(values=c("#1dabe6","#b35a2d","black","#665191"))
Fig_res_amb_ffs<-Fig_res_amb_ffs+labs(color="")+guides(size=FALSE)
Fig_res_amb_ffs
Fig_res_amb_ffs<-Fig_res_amb_ffs+theme_bw()+geom_segment(aes(x=1, y=0, xend=36, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21,24,27,30,33,36))+
  ylim(-3,2)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+
  #ylab(expression(Resilience~(ln~Litterfall~t[x]~t[0]^-1)))+xlab("Time since disturbance (Months)")+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text.y=element_text(size=22),legend.text =  element_text(size=26),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=22),legend.box="horizontal",legend.position="none")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+guides(color = guide_legend(override.aes = list(size = 1.5)))+
  annotate("text", x =5, y = 2, fontface="bold",label = "d FFS fall", size=8,colour="black")+labs(x="")
Fig_res_amb_ffs<-Fig_res_amb_ffs+xlab("Time since disturbance (Months)")
Fig_res_amb_ffs

#Final Figure 6####
Fig6v2<-Fig_res_amb_tot+Fig_res_amb_wood+Fig_res_amb_leaf+Fig_res_amb_ffs+plot_layout(ncol=2,heights=c(1,1))+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig6v2

#Saving in High Res
ggsave(filename = "Fig6v4_Resilience_mass.png",
       plot = Fig6v2, width = 23, height = 14, units = 'cm',
       scale = 2, dpi = 1000)

##END###
