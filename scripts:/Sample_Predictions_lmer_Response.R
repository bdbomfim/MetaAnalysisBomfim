##Sample script predictions lmer with response data

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
