##Leaf fall Resilience Analysis

leaf_all<-data_esall %>% filter(Fraction=="Leaf fall")
str(leaf_all)#358 obs
levels(leaf_all$Treatment)

#1 to 36 months - Ambient only
leaf_amb<-leaf_all%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")
str(leaf_amb)#278 obs
leaf_amb$TSD_months=as.numeric(leaf_amb$TSD_months)
leaf_amb$Case_study= paste(leaf_amb$Site, leaf_amb$DisturbanceName, sep="|")

##Subseting the data 1 - 21 months by deleting cyclones Ivor, Jova, Patricia and Gilbert and Study #4
leaf_amb_1to21<- leaf_amb  %>% filter (TSD_months<22) %>% filter(Study_ID!="4")%>% filter(DisturbanceName!="Keith")%>% filter(DisturbanceName!="Ivor")%>% filter(DisturbanceName!="Jova")%>% filter(DisturbanceName!="Patricia")%>% filter(DisturbanceName!="Gilbert")
str(leaf_amb_1to21)#193
unique(levels(as.factor(leaf_amb_1to21$Country)))

##Fitting random-effects meta-analysis model to obtain weights

#Transforming the numerical predictors
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
leaf_amb_1to21$soilP<-z.trans(leaf_amb_1to21$Other_soil_P)
leaf_amb_1to21$hurrwind<-z.trans(leaf_amb_1to21$HURRECON_wind_ms)
leaf_amb_1to21$windur<-z.trans(leaf_amb_1to21$Gale_wind_duration_minutes)
leaf_amb_1to21$tsd<-z.trans(leaf_amb_1to21$TSD_months)

#Weights####
unique(levels(as.factor(leaf_amb_1to21$Case_study)))

#This is how weights are calculated in multilevel models
#wi <- 1 / (sum(model$sigma2) + data$vi)
#wi

#Fitting multilevel meta-analysis models with varying possible random effects
leaf_meta<- rma.mv(yi,vi,random = list(~1|Site,~1|DisturbanceName),
                  tdist = TRUE,
                  data = leaf_amb_1to21,struct = "HAR",method = "REML")
summary(leaf_meta)
weight_lf<-weights(leaf_meta,type="matrix")
leaf_amb_1to21$weight<-colSums(weight_lf)/sum(weight_lf)
leaf_amb_1to21$weight
#Best method
leaf_meta$sigma2
leaf_amb_1to21$weight2<-(1/(leaf_amb_1to21$vi+0.0464996+0.3230058))
leaf_amb_1to21$weight2

leaf_meta_b<- rma.mv(yi,vi,random = ~1|Country/DisturbanceName,
                    tdist = TRUE,
                    data = leaf_amb_1to21,struct = "HAR",method = "REML")
summary(leaf_meta_b)

leaf_amb_1to21$weight_b<-weights(leaf_meta_b)
weight_lfb<-weights(leaf_meta_b,type="matrix")
leaf_amb_1to21$weight_b1<-colSums(weight_lfb)/sum(weight_lfb)
leaf_amb_1to21$weight_b1
#Best method
leaf_meta_b$sigma2
leaf_amb_1to21$weight_b2<-(1/(leaf_amb_1to21$vi+0.02965614+0.19747667))
leaf_amb_1to21$weight_b2

leaf_meta_c<- rma.mv(yi,vi,random = list(~1|Country,~1|DisturbanceName),
                     tdist = TRUE,
                     data = leaf_amb_1to21,struct = "HAR",method = "REML")
summary(leaf_meta_c)
leaf_amb_1to21$weight_c<-weights(leaf_meta_c)
wi <- 1 / (sum(res$sigma2) + dat$vi)

weight_lfc<-weights(leaf_meta_c,type="matrix")
leaf_amb_1to21$weight_c1<-colSums(weight_lfc)/sum(weight_lfc)
leaf_amb_1to21$weight_c1
#Best method
leaf_meta_c$sigma2
leaf_amb_1to21$weight_c2<-1 / (sum(leaf_meta_c$sigma2) + leaf_amb_1to21$vi)
leaf_amb_1to21$weight_c2


leaf_meta_d<- rma.mv(yi,vi,random = ~1|DisturbanceName,
                     tdist = TRUE,
                     data = leaf_amb_1to21,struct = "HAR",method = "REML")
summary(leaf_meta_d)
leaf_amb_1to21$weight_d<-weights(leaf_meta_d)

weight_lfd<-weights(leaf_meta_d,type="matrix")
leaf_amb_1to21$weight_d1<-colSums(weight_lfd)/sum(weight_lfd)
leaf_amb_1to21$weight_d1
#Best method
leaf_meta_d$sigma2
leaf_amb_1to21$weight_d2<-(1/(leaf_amb_1to21$vi+0.1837246))
leaf_amb_1to21$weight_d2

leaf_meta_e<- rma.mv(yi,vi,random = ~1|Country,
                     tdist = TRUE,
                     data = leaf_amb_1to21,struct = "HAR",method = "REML")
summary(leaf_meta_e)
leaf_amb_1to21$weight_e<-weights(leaf_meta_e)

weight_lfe<-weights(leaf_meta_e,type="matrix")
leaf_amb_1to21$weight_e1<-colSums(weight_lfe)/sum(weight_lfe)
leaf_amb_1to21$weight_e1
#Best method
leaf_meta_e$sigma2
leaf_amb_1to21$weight_e2<-(1/(leaf_amb_1to21$vi+0.2090469))
leaf_amb_1to21$weight_e2

##Fitting GAMM model based on the best model for total litterfall
unique(levels(as.factor(leaf_amb_1to21$Case_study)))

names(leaf_amb_1to21)
#This is the best model for interaction between soil P and time since cyclone
gamm_2y_mixed_1l <- gamm4(yi ~ s(tsd,soilP)                                          ,weights=leaf_amb_1to21$weight_c, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
gamm_2y_mixed_1.1l <- gamm4(yi ~ s(tsd,soilP)                                          ,weights=leaf_amb_1to21$weight_c1, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
gamm_2y_mixed_1.1l <- gamm4(yi ~ s(tsd,soilP)                                          ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
mods_1l <- list(mixed_1l = gamm_2y_mixed_1l$mer, mixed_1.1l = gamm_2y_mixed_1.1l$mer,mixed_1.2l = gamm_2y_mixed_1.2l$mer)
atab_1l <- aictab(mods_1l)
atab_1l

gamm_2y_mixed_1l <- gamm4(yi ~ s(tsd,soilP)                                          ,weights=leaf_amb_1to21$weight, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
gamm_2y_mixed_1l <- gamm4(yi ~ s(tsd,soilP)                                          ,weights=leaf_amb_1to21$weight, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)

summary(gamm_2y_mixed_1l$gam)
summary(gamm_2y_mixed_1l$mer)

gamm_2y_mixed_1.1l <- gamm4(yi ~ s(tsd,hurrwind)                           ,weights=leaf_amb_1to21$weight, random = ~(1|Site)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1.1l$gam)
summary(gamm_2y_mixed_1.1l$mer)

gamm_2y_mixed_1.1l <- gamm4(yi ~ HURRECON_wind_ms+TSD_months                          ,weights=leaf_amb_1to21$weight, random = ~(1|Site)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1.1l$gam)#R2=018
summary(gamm_2y_mixed_1.1l$mer)

gamm_2y_mixed_1.2l <- gamm4(yi ~ hurrwind*tsd                          ,weights=leaf_amb_1to21$weight, random = ~(1|Site)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1.2l$gam)#R2=021
summary(gamm_2y_mixed_1.2l$mer)

gamm_2y_mixed_1.3l <- gamm4(yi ~ s(exp(soilP),tsd)                         ,weights=leaf_amb_1to21$weight, random = ~(1|Site)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1.3l$gam)#R2=0
summary(gamm_2y_mixed_1.3l$mer)

gamm_2y_mixed_1bl <- gamm4(yi ~ soilP*tsd                          ,weights=leaf_amb_1to21$weight_b, random = ~(1|Country/DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1bl$gam)
summary(gamm_2y_mixed_1bl$mer)

gamm_2y_mixed_1.1cl <- gamm4(yi ~ s(soilP,tsd)                         ,weights=leaf_amb_1to21$weight_c, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_1.1cl$gam)
summary(gamm_2y_mixed_1.1cl$mer)

gamm_2y_mixed_2cl <- gamm4(yi ~ s(tsd,by=hurrwind)                           ,weights=leaf_amb_1to21$weight_c, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_2cl$gam)#0.33
summary(gamm_2y_mixed_2cl$mer)


gamm_2y_mixed_2cl.1 <- gamm4(yi ~ s(tsd,by=hurrwind)                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.1$gam)
summary(gamm_2y_mixed_2cl.1$mer)

gamm_2y_mixed_2cl.2 <- gamm4(yi ~ s(tsd,hurrwind)                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.2$gam)

gamm_2y_mixed_2cl.2a <- gamm4(yi ~ s(tsd,hurrwind,k=20)                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.2$gam)

gamm_2y_mixed_2cl.3 <- gamm4(yi ~ s(tsd) + hurrwind                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.3$gam)

gamm_2y_mixed_2cl.4 <- gamm4(yi ~ tsd + hurrwind                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.4$gam)

gamm_2y_mixed_2cl.5 <- gamm4(yi ~ tsd*hurrwind                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=F)
summary(gamm_2y_mixed_2cl.5$gam)

#COMPARING MODELS FOR Table 5 Model 2b
mods_2b <- list(mixed_2cl.1 = gamm_2y_mixed_2cl.1$mer, mixed_2cl.2 = gamm_2y_mixed_2cl.2$mer, mixed_2cl.2a = gamm_2y_mixed_2cl.2a$mer,
                mixed_2cl.3 = gamm_2y_mixed_2cl.3$mer, mixed_2cl.4 = gamm_2y_mixed_2cl.4$mer, mixed_2cl.5 = gamm_2y_mixed_2cl.5$mer)
atab_2b <- aictab(mods_2b)
atab_2b

#Table5 Model2b####
gamm_2y_mixed_2cl.2final <- gamm4(yi ~ s(tsd,hurrwind,k=20)                           ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=TRUE)
summary(gamm_2y_mixed_2cl.2final$gam)
summary(gamm_2y_mixed_2cl.2final$mer)

#Table 5 Model 1b####
gamm_2y_mixed_1.1lfinal <- gamm4(yi ~ s(soilP,tsd,k=20)                                          ,weights=leaf_amb_1to21$weight_c2, random = ~(1|Country)+(1|DisturbanceName), data=leaf_amb_1to21, REML=T)
summary(gamm_2y_mixed_1.1lfinal$gam)
summary(gamm_2y_mixed_1.1lfinal$mer)
names(leaf_amb_1to21)

#Data wrangling for plotting - Best Overall Model
mypreds_1to21_l<-predict(gamm_2y_mixed_2cl.2final$gam,newdata=leaf_amb_1to21,se.fit=T)
leaf_amb_1to21$pred_l<-mypreds_1to21_l$fit
leaf_amb_1to21$se_l<-mypreds_1to21_l$se.fit
names(leaf_amb_1to21)

##Data wrangling - Time since disturbance as factor for plotting
leaf_amb_1to21$Months<-factor(leaf_amb_1to21$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                                      "18","19","20","21"))
levels(leaf_amb_1to21$Months)

#Figure with the best overall model for leaf fall resilience
Fig_res1to21_l <- ggplot(leaf_amb_1to21, aes(y=pred_l, x=Months,group=Case_study))
Fig_res1to21_l<-Fig_res1to21_l+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.9,size=1.2)
Fig_res1to21_l<-Fig_res1to21_l+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(leaf_amb_1to21$vi)),alpha=0.9,shape=21,stroke=1.3)#+ scale_size_continuous(range = c(1, 5))#size=1.5)#color="darkgray")
Fig_res1to21_l<-Fig_res1to21_l+geom_ribbon(aes(ymin=pred_l-(1.96*se_l),ymax=pred_l+(1.96*se_l),col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.02,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,5.5,6,6.5))+scale_fill_gradient(low="#FCFF00",high="#D8001F",breaks=c(5,5.5,6,6.5))#+scale_color_manual(values=c("#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#ef9250","#e9774c","#d43d51"))+scale_fill_manual(values=c("#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#ef9250","#e9774c","#d43d51"))#+ scale_color_manual(values=c("#FFD700","#A07717","#F06E05","#268FE1","#265E89","#192DB5","#85EE85","#409E40","#0C620C","#D37EA5","#D14082","#EF0672"))+ scale_fill_manual(values=c("#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#EF0672"))
Fig_res1to21_l<-Fig_res1to21_l+geom_line(aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
Fig_res1to21_l
Fig_res1to21_l<-Fig_res1to21_l+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21))+scale_shape_discrete(solid=F)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("Time since disturbance (Months)")+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text=element_text(size=22),legend.text =  element_text(size=20),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=20),legend.box="horizontal",legend.position="top")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+labs(color="Soil P \n(ln mg/kg)")+guides(fill=FALSE,size=FALSE)+
annotate("text", x = 1, y = 1.2, label = "b Leaf fall", size=8,hjust=0,colour="black",fontface="bold")
Fig_res1to21_l

#Data wrangling for plotting - Best Model with Soil P
mypreds_1to21_l2<-predict(gamm_2y_mixed_1.1lfinal$gam,newdata=leaf_amb_1to21,se.fit=T)
leaf_amb_1to21$pred_l2<-mypreds_1to21_l2$fit
leaf_amb_1to21$se_l2<-mypreds_1to21_l2$se.fit
names(leaf_amb_1to21)

##Data wrangling - Time since disturbance as factor for plotting
leaf_amb_1to21$Months<-factor(leaf_amb_1to21$TSD_months, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                                                    "18","19","20","21"))
levels(leaf_amb_1to21$Months)

#Figure with the best model for leaf fall resilience with soil P and Time since cyclone
Fig_res1to21_l2 <- ggplot(leaf_amb_1to21, aes(y=pred_l2, x=Months,group=Case_study))
Fig_res1to21_l2<-Fig_res1to21_l2+geom_point(aes(group=Case_study,col=log(Other_soil_P)),alpha=0.9,size=1.2)
Fig_res1to21_l2<-Fig_res1to21_l2+geom_point(aes(y=yi,x=Months,group=Case_study,col=log(Other_soil_P),size=sqrt(leaf_amb_1to21$vi)),alpha=0.9,shape=21,stroke=1.3)#+ scale_size_continuous(range = c(1, 5))#size=1.5)#color="darkgray")
Fig_res1to21_l2<-Fig_res1to21_l2+geom_ribbon(aes(ymin=pred_l2-(1.96*se_l2),ymax=pred_l2+(1.96*se_l2),col=log(Other_soil_P),fill=log(Other_soil_P)),alpha=0.05,linetype=3,size=0.4)+scale_color_gradient(low="#FCFF00",high="#D8001F")+scale_fill_gradient(low="#FCFF00",high="#D8001F")#+scale_color_manual(values=c("#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#ef9250","#e9774c","#d43d51"))+scale_fill_manual(values=c("#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#ef9250","#e9774c","#d43d51"))#+ scale_color_manual(values=c("#FFD700","#A07717","#F06E05","#268FE1","#265E89","#192DB5","#85EE85","#409E40","#0C620C","#D37EA5","#D14082","#EF0672"))+ scale_fill_manual(values=c("#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#F7EAEA","#EF0672"))
Fig_res1to21_l2<-Fig_res1to21_l2+geom_line(aes(group=Case_study,col=log(Other_soil_P)),size=1.2)
#Fig_res1to21_l2
Fig_res1to21_l2<-Fig_res1to21_l2+theme_pubr()+geom_segment(aes(x=1, y=0, xend=21, yend=0), lty=2, color = "magenta", cex=1.4)+
  scale_x_discrete(breaks = c(1, 3, 6, 9,12, 15, 18,21))+scale_shape_discrete(solid=F)+
  ylab(expression(Resilience~(ln~litterfall[tx]/litterfall[t0])))+xlab("Time since disturbance (Months)")+guides(color = guide_colourbar(barwidth = 15, barheight = 1.5,nbin=30,ticks.colour="black",ticks.linewidth = 2.5))+
  theme(axis.title.x =element_text(vjust = 0.5,size=24),
        axis.text.x =element_text(vjust = 1,size=22),
        axis.title.y =element_text(vjust = 1,size=24),strip.background = element_rect(color="white", fill="white",linetype="solid"),
        axis.text=element_text(size=22),legend.text =  element_text(size=28),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),
        legend.title = element_text(size=24),legend.box="horizontal",legend.position="top")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+labs(color="Soil P \n(ln mg/kg)")+guides(fill=FALSE,size=FALSE)+
annotate("text", x = 1, y = 1.2, label = "c Leaf fall - soil P, time", size=8,hjust=0,colour="black",fontface="bold")
Fig_res1to21_l2

ggsave(filename = "Fig6b_Pred_Leaf.png",
       plot = Fig_res1to21_l_facet, width = 24, height = 12, units = 'cm',
       scale = 2, dpi = 600)

##Visualizing the total litterfall data 1 to 36 months - Ambient only
dres_amb_l <- leaf_amb %>%group_by(Case_study)  %>%
  dplyr::summarise(counts = dplyr::n())
dres_amb_l
dres_amb_l$counts=as.numeric(dres_amb_l$counts)
dres_amb_l
p_data_clus_l<- dres_amb_l %>% ggplot(aes(x = reorder(Case_study,-counts), y=counts)) +geom_point(stroke=2,color="blue",size=4,shape=21)+theme_bw()+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text.y=element_text(size=18),axis.text.x = element_text(angle=52, hjust=1,vjust = 1,size=18),axis.title=element_text(size=20),
        legend.position="right",legend.text =  element_text(size=20,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))+geom_hline(aes(yintercept=1), lty=2, color = "red", cex=0.9, alpha = .8)+
  #scale_fill_brewer(palette="Paired")+
  labs(x="Case study",y="Number of observations")#+ annotate("text", x = 10, y = 25, label = "3-year resilience of total litterfall",size=6,colour="black")
p_data_clus_l
#saving impage
ggsave(filename = "SupFig_TotAmbRes_Obs.png",
       plot = FigSupTotAmb, width = 22, height = 19, units = 'cm',
       scale = 2, dpi = 600)

##Checking time scale
dres_amb_l2 <- leaf_amb %>%group_by(Case_study,Treatment,TSD_months)  %>%
  dplyr::summarise(counts = dplyr::n())
dres_amb_l2
p_amb_l<- dres_amb_l2 %>% ggplot(aes(x = TSD_months, fill=Treatment)) +geom_histogram(binwidth = 0.5)+theme_bw()+guides(color = guide_legend(title = "Country"),legend.key.width=32)+
  theme(strip.background = element_rect(color="white", fill="white",linetype="solid"),axis.text= element_text(angle=0, hjust=0.5,vjust = 1,size=18),axis.title=element_text(size=22),
        legend.position="right",legend.text =  element_text(size=20,angle=0),legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)),legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))+scale_fill_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))+
  #scale_fill_brewer(palette="Paired")+
  labs(x="Time since disturbance (months)",y="Number of observations")#+ annotate("text", x = 10, y = 25, label = "3-year resilience of total litterfall",size=6,colour="black")
FigSup_Amb_l<-p_amb_l + facet_wrap(~Case_study)+theme(strip.text.x = element_text(size=16))
FigSup_Amb_l

##END##
