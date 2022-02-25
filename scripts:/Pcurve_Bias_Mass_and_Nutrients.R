##Resistance and Resilience P-curve (Bias) Analysis###

library(metafor)
library(dmetar)

##Resistance####
#Data wrangling to prepare dataframes
data0a<-metadat %>% filter(Fraction=="TotLitfall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
data_es0ia <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = data0a, measure = "ROM")

#Mass flux - Supplementary Table S3######

#Total litterfall 
str(data_es0ia)
pcurvedata<-cbind(data.frame(TE = data_es0ia$yi, seTE=data_es0ia$vi, studlab=data_es0ia$Study_ID))
pcurvedata
pcurve(pcurvedata)
pcurve(pcurvedata, effect.estimation = TRUE, N = data_es0ia$S_size, dmin = -1, dmax = 7)

#Leaf litterfall
pcurvedatalf<-cbind(data.frame(TE = data_es0ilf$yi, seTE=data_es0ilf$vi, studlab=data_es0ilf$Study_ID))
pcurvedatalf
pcurve(pcurvedatalf)
pcurve(pcurvedatalf, effect.estimation = TRUE, N = data_es0ia$S_size, dmin = -1, dmax = 7)

#Wood litterfall
pcurvedatawf<-cbind(data.frame(TE = data_es0iwf$yi, seTE=data_es0iwf$vi, studlab=data_es0iwf$Study_ID))
pcurvedatawf
pcurve(pcurvedatawf)
pcurve(pcurvedatalf, effect.estimation = TRUE, N = data_es0ia$S_size, dmin = -1, dmax = 7)

#FFS fall
data_es0iff_f<-data_es0iff %>% filter(yi!="NA")
pcurvedataff<-cbind(data.frame(TE = data_es0iff_f$yi, seTE=data_es0iff_f$vi, studlab=data_es0iff_f$Study_ID))
pcurvedataff
pcurve(pcurvedataff)

#Miscellaneous fall
pcurvedatamf<-cbind(data.frame(TE = data_es0imf$yi, seTE=data_es0imf$vi, studlab=data_es0imf$Study_ID))
pcurvedatamf
pcurve(pcurvedatamf)

##P curve analysis for bias in Nutrient resistance####

#P flux - Supplementary Table S4####

#total litterfall
str(data_es0ia)
pcurvedata_tpf<-cbind(data.frame(TE = data_es0itpf$yi, seTE=data_es0itpf$vi, studlab=data_es0itpf$Study_ID))
pcurvedata_tpf
pcurve(pcurvedata_tpf)

##Leaf litterfall
pcurvedata_lpf<-cbind(data.frame(TE = data_es0ilpf$yi, seTE=data_es0ilpf$vi, studlab=data_es0ilpf$Study_ID))
pcurvedata_lpf
pcurve(pcurvedata_lpf)

pcurvedata_wpf<-cbind(data.frame(TE = data_es0iwpf$yi, seTE=data_es0iwpf$vi, studlab=data_es0iwpf$Study_ID))
pcurvedata_wpf
pcurve(pcurvedata_wpf)

#N flux - Supplementary Table S4####

#Total litterfall N flux
pcurvedata_tnf<-cbind(data.frame(TE = data_es0itnf$yi, seTE=data_es0itnf$vi, studlab=data_es0itnf$Study_ID))
pcurvedata_tnf
pcurve(pcurvedata_tnf)

##Leaf litterfall N flux
pcurvedata_lnf<-cbind(data.frame(TE = data_es0ilnf$yi, seTE=data_es0ilnf$vi, studlab=data_es0ilnf$Study_ID))
pcurvedata_lnf
pcurve(pcurvedata_lnf)

#Wood litterfall N flux
pcurvedata_wnf<-cbind(data.frame(TE = data_es0iwnf$yi, seTE=data_es0iwnf$vi, studlab=data_es0iwnf$Study_ID))
pcurvedata_wnf
pcurve(pcurvedata_wnf)

#P concentration - Supplementary Table S4####
#Leaf litterfall
pcurvedata_lpc<-cbind(data.frame(TE = data_es0ilpc$yi, seTE=data_es0ilpc$vi, studlab=data_es0ilpc$Study_ID))
pcurvedata_lpc
pcurve(pcurvedata_lpc)

#N concentration - Supplementary Table S4####
#Leaf litterfall
pcurvedata_lnc<-cbind(data.frame(TE = data_es0ilnc$yi, seTE=data_es0ilnc$vi, studlab=data_es0ilnc$Study_ID))
pcurvedata_lnc
pcurve(pcurvedata_lnc)

##Resilience####
#Data wrangling
res_all_21<- metadat %>% filter(Cat_TSD_months == "Rec") %>% filter(Case_ID!="25.2")%>% filter(Case_ID!="18.1")%>% 
  filter (TSD_months < 22)%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")%>% filter(Fraction!="Misc fall")
str(res_all_21)#731 all fractions
data_esall_21 <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                         sd1i = Post_SD, sd2i = Pre_SD, data = res_all_21, measure = "ROM")
#New resilience metric
data_esall_21$yi_new <- data_esall_21$yi / data_esall_21$TSD_months
#checking new column yi_new
summary(data_esall_21$yi_new)
#Same for the variance vi_new
data_esall_21$vi_new <- data_esall_21$vi / data_esall_21$TSD_months
#checking new column vi_new
summary(data_esall_21$vi_new)

##Mass Flux All Fractions - Supplementary Table S5####

#Total litterfall###
tot_lit_amb_1to21_final<-data_esall_21 %>% filter(Fraction=="TotLitfall")%>% filter(DisturbanceName!="Keith")%>% filter(Site!="San Felipe")%>% filter (Site!="Grande-Terre")%>% filter(DisturbanceName!="Ivor")
str(tot_lit_amb_1to21_final)#213 obs amb+CTE,
pcurve_tot_mass<-cbind(data.frame(TE = tot_lit_amb_1to21_final$yi_new, seTE=tot_lit_amb_1to21_final$vi_new, studlab=tot_lit_amb_1to21_final$Study_ID))
pcurve_tot_mass
pcurve(pcurve_tot_mass)

#Testing without CTE#
tot_lit_amb_1to21_noCTE<-tot_lit_amb_1to21_final %>% filter(Treatment!="TrimDeb")
str(tot_lit_amb_1to21_noCTE)#192 observations
pcurve_tot_mass_noCTE<-cbind(data.frame(TE = tot_lit_amb_1to21_noCTE$yi_new, seTE=tot_lit_amb_1to21_noCTE$vi_new, studlab=tot_lit_amb_1to21_noCTE$Study_ID))
pcurve_tot_mass_noCTE
pcurve(pcurve_tot_mass_noCTE)

##Leaf litterfall###
pcurve_leaf_mass<-cbind(data.frame(TE = leaf_amb_1to21$yi, seTE=leaf_amb_1to21$vi, studlab=leaf_amb_1to21$Study_ID))
pcurve_leaf_mass
pcurve(pcurve_leaf_mass)

leaf_1to21_final<-data_esall_21 %>% filter(Fraction=="Leaf fall")
summary(leaf_1to21_final$Treatment)
pcurve_leaf_mass_new<-cbind(data.frame(TE = leaf_1to21_final$yi_new, seTE=leaf_1to21_final$vi_new, studlab=leaf_1to21_final$Study_ID))
pcurve_leaf_mass_new
pcurve(pcurve_leaf_mass_new)

##Wood litterfall#####

#1 to 21 months subset
wood_amb_1to21<-data_esall %>% filter(Fraction=="Wood fall") %>% filter (TSD_months<22)%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")

pcurve_wood_mass<-cbind(data.frame(TE = wood_amb_1to21$yi, seTE=wood_amb_1to21$vi, studlab=wood_amb_1to21$Study_ID))
pcurve_wood_mass
pcurve(pcurve_wood_mass)

##FFS fall mass flux#####

#1 to 21 months subset
ffs_amb_1to21<-data_esall %>% filter(Fraction=="FFS fall") %>% filter (TSD_months<22)%>% filter(Treatment=="Ambient"|Treatment=="TrimDeb")

pcurve_ffs_mass<-cbind(data.frame(TE = ffs_amb_1to21$yi, seTE=ffs_amb_1to21$vi, studlab=ffs_amb_1to21$Study_ID))
pcurve_ffs_mass
pcurve(pcurve_ffs_mass)

##N and P Flux All Fractions - Supplementary Table S6####

##P flux####

#1 to 15 months
data_esall_amb_Pflux_1to15<-data_esall_amb_Pflux %>% filter(TSD_months<16)
data_esall_amb_Pflux_1to15$Fraction

#Total P flux##
tot_Pflux_1to15<-data_esall_amb_Pflux_1to15 %>% filter(Fraction == "Total")

pcurve_tot_pflux<-cbind(data.frame(TE = tot_Pflux_1to15$yi, seTE=tot_Pflux_1to15$vi, studlab=tot_Pflux_1to15$Study_ID))
pcurve_tot_pflux
pcurve(pcurve_tot_pflux)

#Leaf P flux###
leaf_Pflux_1to15<-data_esall_amb_Pflux_1to15 %>% filter(Fraction == "Leaf")

pcurve_leaf_pflux<-cbind(data.frame(TE = leaf_Pflux_1to15$yi, seTE=leaf_Pflux_1to15$vi, studlab=leaf_Pflux_1to15$Study_ID))
pcurve_leaf_pflux
pcurve(pcurve_leaf_pflux)

#Wood P flux###
wood_Pflux_1to15<-data_esall_amb_Pflux_1to15 %>% filter(Fraction == "Wood")

pcurve_wood_pflux<-cbind(data.frame(TE = wood_Pflux_1to15$yi, seTE=wood_Pflux_1to15$vi, studlab=wood_Pflux_1to15$Study_ID))
pcurve_wood_pflux
pcurve(pcurve_wood_pflux)

#FFS P flux###
ffs_Pflux_1to15<-data_esall_amb_Pflux_1to15 %>% filter(Fraction == "FFS")

pcurve_ffs_pflux<-cbind(data.frame(TE = ffs_Pflux_1to15$yi, seTE=ffs_Pflux_1to15$vi, studlab=ffs_Pflux_1to15$Study_ID))
pcurve_ffs_pflux
pcurve(pcurve_ffs_pflux)

##N flux####

#1 to 15 months
data_esall_amb_Nflux_1to15<-data_esall_amb_Nflux %>% filter(TSD_months<16)
data_esall_amb_Nflux_1to15$Fraction

#Total N flux##

tot_Nflux_1to15<-data_esall_amb_Nflux_1to15 %>% filter(Fraction == "Total")

pcurve_tot_nflux<-cbind(data.frame(TE = tot_Nflux_1to15$yi, seTE=tot_Nflux_1to15$vi, studlab=tot_Nflux_1to15$Study_ID))
pcurve_tot_nflux
pcurve(pcurve_tot_nflux)

#Leaf N flux###
leaf_Nflux_1to15<-data_esall_amb_Nflux_1to15 %>% filter(Fraction == "Leaf")

pcurve_leaf_nflux<-cbind(data.frame(TE = leaf_Nflux_1to15$yi, seTE=leaf_Nflux_1to15$vi, studlab=leaf_Nflux_1to15$Study_ID))
pcurve_leaf_nflux
pcurve(pcurve_leaf_nflux)

#Wood N flux###
wood_Nflux_1to15<-data_esall_amb_Nflux_1to15 %>% filter(Fraction == "Wood")

pcurve_wood_nflux<-cbind(data.frame(TE = wood_Nflux_1to15$yi, seTE=wood_Nflux_1to15$vi, studlab=wood_Nflux_1to15$Study_ID))
pcurve_wood_nflux
pcurve(pcurve_wood_nflux)

#FFS N flux###
ffs_Nflux_1to15<-data_esall_amb_Nflux_1to15 %>% filter(Fraction == "FFS")

pcurve_ffs_nflux<-cbind(data.frame(TE = ffs_Nflux_1to15$yi, seTE=ffs_Nflux_1to15$vi, studlab=ffs_Nflux_1to15$Study_ID))
pcurve_ffs_nflux
pcurve(pcurve_ffs_nflux)

##P and N concentration - Supplementary Table S7####

##P concentration####

#1 to 15 months
data_esall_amb_Pconc_1to15<-data_esall_amb_Pc %>% filter(TSD_months<16)#data from script Resilience_Nutrient_Plots
data_esall_amb_Pconc_1to15$Fraction

##Leaf P concentration###

leaf_Pc_1to15<-data_esall_amb_Pconc_1to15 %>% filter(Fraction == "Leaf")

pcurve_leaf_Pc<-cbind(data.frame(TE = leaf_Pc_1to15$yi, seTE=leaf_Pc_1to15$vi, studlab=leaf_Pc_1to15$Study_ID))
pcurve_leaf_Pc
pcurve(pcurve_leaf_Pc)

##Leaf N concentration##

data_esall_amb_Nconc_1to15<-data_esall_amb_Nc %>% filter(TSD_months<16)#data from script Resilience_Nutrient_Plots
data_esall_amb_Nconc_1to15$Fraction

leaf_Nc_1to15<-data_esall_amb_Nconc_1to15 %>% filter(Fraction == "Leaf")

pcurve_leaf_Nc<-cbind(data.frame(TE = leaf_Nc_1to15$yi, seTE=leaf_Nc_1to15$vi, studlab=leaf_Nc_1to15$Study_ID))
pcurve_leaf_Nc
pcurve(pcurve_leaf_Nc)

##END##