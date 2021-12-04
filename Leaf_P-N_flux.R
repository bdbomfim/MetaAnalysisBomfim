##Leaf P flux##

leaf_Pflux<-data_esall_amb_Pflux %>% filter(Fraction=="Leaf")
str(leaf_Pflux)#92
#1 and 2 months time bins
data_es_leafPflux_1<-leaf_Pflux %>% filter(TSD_months<3)

leafPflux_meta_1<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafPflux_1,struct = "HAR",method = "REML")
summary(leafPflux_meta_1)
-(1-exp(leafPflux_meta_1$b))*100 #76.9% below
((exp(leafPflux_meta_1$se))-1)*100

#3 and 4 months
data_es_leafPflux_3<-leaf_Pflux %>% filter(TSD_months==3|TSD_months==4)
str(data_es_leafPflux_3)

leafPflux_meta_3<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafPflux_3,struct = "HAR",method = "REML")
summary(leafPflux_meta_3)

#5 and 6 months
data_es_leafPflux_5<-leaf_Pflux %>% filter(TSD_months==5|TSD_months==6)
str(data_es_leafPflux_5)

leafPflux_meta_5<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafPflux_5,struct = "HAR",method = "REML")
summary(leafPflux_meta_5)

#8 and 9 months
data_es_leafPflux_8<-leaf_Pflux %>% filter(TSD_months==8|TSD_months==9)
str(data_es_leafPflux_8)

leafPflux_meta_8<- rma.mv(yi,vi,#random = ~(1|Region),
                         tdist = TRUE,
                         data = data_es_leafPflux_8,struct = "HAR",method = "REML")
summary(leafPflux_meta_8)

#11 and 12 months
data_es_leafPflux_12<-leaf_Pflux %>% filter(TSD_months==11|TSD_months==12)
str(data_es_leafPflux_12)

leafPflux_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_leafPflux_12,struct = "HAR",method = "REML")
summary(leafPflux_meta_12)

#14 and 15 months
data_es_leafPflux_15<-leaf_Pflux %>% filter(TSD_months==14|TSD_months==15)
str(data_es_leafPflux_15)

leafPflux_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_leafPflux_15,struct = "HAR",method = "REML")
summary(leafPflux_meta_15)

##Leaf N flux##

#Total
leaf_Nflux<-data_esall_amb_Nflux %>% filter(Fraction=="Leaf")
str(leaf_Nflux)#92

#1 and 2 months
data_es_leafNflux_1<-leaf_Nflux %>% filter(TSD_months<3)
str(data_es_totNflux_1)
leafNflux_meta_1<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafNflux_1,struct = "HAR",method = "REML")
summary(leafNflux_meta_1)
-(1-exp(leafNflux_meta_1$b))*100
((exp(leafNflux_meta_1$se))-1)*100

#3 and 4 months
data_es_leafNflux_3<-leaf_Nflux %>% filter(TSD_months==3|TSD_months==4)

leafNflux_meta_3<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data =  data_es_leafNflux_3,struct = "HAR",method = "REML")
summary(leafNflux_meta_3)
-(1-exp(leafNflux_meta_3$b))*100

#5 and 6 months
data_es_leafNflux_5<-leaf_Nflux %>% filter(TSD_months==5|TSD_months==6)
str(data_es_leafNflux_5)

leafNflux_meta_5<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafNflux_5,struct = "HAR",method = "REML")
summary(leafNflux_meta_5)

#8 and 9 months
data_es_leafNflux_8<-leaf_Nflux %>% filter(TSD_months==8|TSD_months==9)
str(data_es_leafNflux_8))

leafNflux_meta_8<- rma.mv(yi,vi,#random = ~(1|Site),
                         tdist = TRUE,
                         data = data_es_leafNflux_8,struct = "HAR",method = "REML")
summary(leafNflux_meta_8)

#11 and 12 months
data_es_leafNflux_12<-leaf_Nflux %>% filter(TSD_months==11|TSD_months==12)
str(data_es_leafNflux_12)

leafNflux_meta_12<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_leafNflux_12,struct = "HAR",method = "REML")
summary(leafNflux_meta_12)

#14 and 15 months
data_es_leafNflux_15<-leaf_Nflux %>% filter(TSD_months==14|TSD_months==15)
str(data_es_leafNflux_15)

leafNflux_meta_15<- rma.mv(yi,vi,#random = ~(1|Site),
                          tdist = TRUE,
                          data = data_es_leafNflux_15,struct = "HAR",method = "REML")
summary(leafNflux_meta_15)
