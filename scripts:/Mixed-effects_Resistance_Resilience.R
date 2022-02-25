##Running mixed-effects models for resilience with resistance as random effect


#annualfuture_45_2030_tmean = split.data.frame(future_45_2030_tmean, future_45_2030_tmean$futureYear)
#lapply(names(annualfuture_45_2030_tmean), function(x){
  #write.csv(annualfuture_45_2030_tmean[[x]], paste(x, "_45_2030_tmean.csv", sep = ""))
#})
#AUS_historicS002$TMIN[index_max] <- sapply(index_max, function(i) with(AUS_historicS002, mean(c(TMAX[i-1], TMAX[i+1]))))

#tot_lit_amb_1to21_final$Resistance <-
summary(tot_lit_amb_1to21_final$Treatment)

unique_vals_rt <- lapply(data_es0ia$Case_study, unique)
unique_vals_rt

unique_vals_rs<-lapply(tot_lit_amb_1to21_final$Case_study, unique)
unique_vals_rs

#Merging dataframes by Case_study
vetor <- c(1,2,3)
key <- data.frame(vetor=vetor, mat=c('a', 'b', 'c'))
data <- data.frame(id=c('a', 'b', 'a', 'c', 'a'))
data$vector1 <- key[match(data$id, key$mat), 'vetor']
data

#My code to assign the resistance value to the respective case study in the resilience dataframe
Resistance <- data_es0ia$yi
rt<- data.frame(Case_study = data_es0ia$Case_study, Resistance = Resistance)
rt
#Same for leaf litterfall
Resistance_lf <- data_es0ilf$yi
Resistance_lf
rt_lf<- data.frame(Case_study = data_es0ilf$Case_study, Resistance_lf = Resistance_lf)
rt_lf

#Assigning Resistance values to case studies in Total litterfall data frame
tot_lit_amb_1to21_final$Resistance <- rt[match(tot_lit_amb_1to21_final$Case_study, rt$Case_study), 'Resistance']
head(tot_lit_amb_1to21_final)
summary(tot_lit_amb_1to21_final$Resistance)
tot_lit_amb_1to21_final

#Assigning Resistance values to case studies in Leaf litterfall data frame
leaf_amb_1to21$Resistance_lf <- rt_lf[match(leaf_amb_1to21$Case_study, rt_lf$Case_study), 'Resistance_lf']
head(leaf_amb_1to21)
summary(leaf_amb_1to21$Resistance_lf)
head(Data_leaf_amb_1to21)

z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
tot_lit_amb_1to21_final$soilP=z.trans(tot_lit_amb_1to21_final$Other_soil_P)
tot_lit_amb_1to21_final$hurrwind=z.trans(tot_lit_amb_1to21_final$HURRECON_wind_ms)
tot_lit_amb_1to21_final$tsls=z.trans(tot_lit_amb_1to21_final$YearsSinceLastStorm)

tot_lit_amb_1to21_final_analysis<-tot_lit_amb_1to21_final %>% filter(resist!="NA")%>% filter(hurrwind!="NA")%>% filter(TSD_months<13)%>% filter(TSD_months>10)
str(tot_lit_amb_1to21_final_analysis)

mix_model_12<-rma.mv(yi,vi,random = ~ 1 | resist, 
                        tdist = TRUE,data = tot_lit_amb_1to21_final_analysis,
                        method = "REML",mods = ~soilP)
mix_model_12

mix_model_12b<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|DisturbanceName), 
                     tdist = TRUE,data = tot_lit_amb_1to21_final_analysis,
                     method = "REML",mods = ~soilP)
mix_model_12b

mix_model_12c<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|Site), 
                      tdist = TRUE,data = tot_lit_amb_1to21_final_analysis,
                      method = "REML",mods = ~soilP)
mix_model_12c

mix_model_12d<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|Site, ~1|DisturbanceName), 
                      tdist = TRUE,data = tot_lit_amb_1to21_final_analysis,
                      method = "REML",mods = ~soilP)
mix_model_12d

AIC(mix_model_12,mix_model_12b,mix_model_12c,mix_model_12d)

#Same for leaf litterfall
head(Data_leaf_amb_1to21)
Data_leaf_amb_1to21_analysis<-Data_leaf_amb_1to21 %>% filter(resist!="NA")%>% filter(hurrwind!="NA")%>% filter(TSD_months<13)%>% filter(TSD_months>10)
str(Data_leaf_amb_1to21_analysis)

mix_model_12lf<-rma.mv(yi,vi,random = ~ 1 | resist, 
                     tdist = TRUE,data = Data_leaf_amb_1to21_analysis,
                     method = "REML",mods = ~soilP)
mix_model_12lf

mix_model_12blf<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|DisturbanceName), 
                      tdist = TRUE,data = Data_leaf_amb_1to21_analysis,
                      method = "REML",mods = ~soilP)
mix_model_12blf

mix_model_12clf<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|Site), 
                      tdist = TRUE,data = Data_leaf_amb_1to21_analysis,
                      method = "REML",mods = ~soilP)
mix_model_12clf

mix_model_12dlf<-rma.mv(yi,vi,random = list(~ 1 | resist, ~1|Site, ~1|DisturbanceName), 
                      tdist = TRUE,data = Data_leaf_amb_1to21_analysis,
                      method = "REML",mods = ~soilP)
mix_model_12dlf

AIC(mix_model_12lf,mix_model_12blf,mix_model_12clf,mix_model_12dlf)

#Checking litterfall values
#My code to assign the resistance value to the respective case study in the resilience dataframe
Lit_peak <- data_es0ia$Post_Mean
LP<- data.frame(Case_study = data_es0ia$Case_study, Lit_peak = Lit_peak)
LP

#Assigning values
tot_lit_amb_1to21_final$Lit_peak <- LP[match(tot_lit_amb_1to21_final$Case_study, LP$Case_study), 'Lit_peak']
head(tot_lit_amb_1to21_final)
summary(tot_lit_amb_1to21_final$Lit_peak)
tot_lit_amb_1to21_final

