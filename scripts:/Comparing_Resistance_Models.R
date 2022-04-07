library(AICcmodavg)
library(glmulti)

# Steps to get to R2 calculation
metadat<-read.csv(file.choose())#Litterfall_Mass in data folder
data0a<-metadat %>% filter(Fraction=="TotLitfall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
data_es0ia <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                     sd1i = Post_SD, sd2i = Pre_SD, data = data0a, measure = "ROM")

#Table 2 - Model1a and 2a####
names(data_es0ia)
#standardizing the moderators
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ia$soilP=z.trans(data_es0ia$Other_soil_P)
data_es0ia$hurrwind=z.trans(data_es0ia$HURRECON_wind_ms)
data_es0ia$tsls=z.trans(data_es0ia$YearsSinceLastStorm)
data_es0ia$long=z.trans(data_es0ia$Longitude)
data_es0ia$stofreq=z.trans(data_es0ia$StormFrequencyNorm)
data_es0ia$windur=z.trans(data_es0ia$Gale_wind_duration_minutes)
data_es0ia$matmap=z.trans(data_es0ia$MAT_MAP_x100)
data_es0ia$elev=z.trans(data_es0ia$Elevation_m)

#Filtering out any NAs in the Hurrecon wind speed column
hurr_sites<-data_es0ia %>% filter(hurrwind!="NA")

#Overall effect size
full.model3 <- rma.mv(-1*yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ia,method = "REML")

summary(full.model3)

#Soil P as single moderator in mixed-effects meta-analysis model
model_P<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP)
summary(model_P)

### Model R2 for soil P mixed-effects model
(sum(full.model3$sigma2) - sum(model_P$sigma2)) / sum(full.model3$sigma2)
#R2 = 0.33

#Reversing standardization of estimate
summary(hurr_sites$Other_soil_P) #mean = 322.1 mg/kg
sd(hurr_sites$Other_soil_P) #SD = 257.89
#Multiply estimate by 2x SD
x = -1.665537*(2*257.89)
x #859.0507
#Add the average
x + 322.1 #-536.9507 is the estimate after reversing the standardization

#Reversing standardization of CI
x_CIub = 1.8121*(2*257.89)
x_CIub #934.6449
#Add the average
x_CIub + 322.1 #1256.745

x_CIlb = 1.5189*(2*257.89)
x_CIlb #934.6449
#Add the average
x_CIlb + 322.1 #1105.518

#Reversed estimate is 

# Calculating R2 for soil P and wind speed mixed-effects model
model_a<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                   data = hurr_sites,method = "REML",
                   mods = ~soilP*hurrwind)
summary(model_a)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_a$sigma2)) / sum(full.model3$sigma2)
#0.39

#Same approach for leaf litterfall resistance

data0alf<-metadat %>% filter(Fraction=="Leaf fall")%>%filter(Cat_TSD_months=="0-0.5")%>%filter(Treatment!="TrimDeb")
data_es0ilf <- escalc(n1i = S_size, n2i = S_size, m1i = Post_Mean, m2i = Pre_Mean, 
                      sd1i = Post_SD, sd2i = Pre_SD, data = data0alf, measure = "ROM")

#Table 2 - Model1a and 2a####
names(data_es0ilf)
#standardizing the moderators
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}
data_es0ilf$soilP=z.trans(data_es0ilf$Other_soil_P)
data_es0ilf$hurrwind=z.trans(data_es0ilf$HURRECON_wind_ms)
data_es0ilf$tsls=z.trans(data_es0ilf$YearsSinceLastStorm)
data_es0ilf$long=z.trans(data_es0ilf$Longitude)
data_es0ilf$stofreq=z.trans(data_es0ilf$StormFrequencyNorm)
data_es0ilf$windur=z.trans(data_es0ilf$Gale_wind_duration_minutes)
data_es0ilf$matmap=z.trans(data_es0ilf$MAT_MAP_x100)
data_es0ilf$elev=z.trans(data_es0ilf$Elevation_m)

#Filtering out any NAs in the Hurrecon wind speed column
hurr_siteslf<-data_es0ilf %>% filter(hurrwind!="NA")

#Overall effect size
full.model3lf <- rma.mv(-1*yi, vi,random = list(~ 1 | Site,~1|DisturbanceName),
                      tdist = TRUE, #here we turn ON the Knapp-Hartung adjustment for CIs
                      data = data_es0ilf,method = "REML")

summary(full.model3lf)

#Soil P as single moderator in mixed-effects meta-analysis model
model_Plf<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_siteslf,method = "REML",
                mods = ~soilP)
summary(model_Plf)

### Model R2 for soil P mixed-effects model
(sum(full.model3lf$sigma2) - sum(model_Plf$sigma2)) / sum(full.model3lf$sigma2)
#R2 = 0.05

# Calculating R2 for soil P and wind speed mixed-effects model
model_alf<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_siteslf,method = "REML",
                mods = ~soilP*hurrwind)
summary(model_alf)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3lf$sigma2) - sum(model_alf$sigma2)) / sum(full.model3lf$sigma2)
#0.39

## Testing alternative models

model_b<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                   data = hurr_sites,method = "REML",
                   mods = ~soilP*long)
summary(model_b)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_b$sigma2)) / sum(full.model3$sigma2)
#0.42

model_c<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*stofreq)
summary(model_c)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_c$sigma2)) / sum(full.model3$sigma2)
#0.54

model_d<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*windur)
summary(model_d)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_d$sigma2)) / sum(full.model3$sigma2)
#0.21

model_e<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP*tsls)
summary(model_e)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_e$sigma2)) / sum(full.model3$sigma2)
#-7.49

model_f<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP + matmap)
summary(model_f)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_f$sigma2)) / sum(full.model3$sigma2)
#-10.76

model_g<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP + Par_Mat)
summary(model_g)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_g$sigma2)) / sum(full.model3$sigma2)
#0.77

model_h<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP + RockP_class)
summary(model_h)
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_h$sigma2)) / sum(full.model3$sigma2)
#0.33

AIC(model_a, model_b, model_c, model_d, model_e, model_f, model_g, model_h)

mods_rt <- list(mixed_a = model_a, mixed_b = model_b)
atab_rt <- aictab(mods_rt)
atab_rt

rma.glmulti<-function(formula, data, ...)
  rma.mv(formula, V=vi, data=hurr_sites, random = list(~ 1 | Site, ~1|DisturbanceName), method="ML", ...)
res<-glmulti(-1*yi~soilP+hurrwind+long, data=hurr_sites, level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=1000)
print(res)
plot(res)

#taking a look at the top models
top <- weightable(res)
top <- top[top$aicc <= min(top$aicc) + 2,]
top
#examining best model
summary(res@objects[[1]])

#multimodel inference
eval(metafor:::.glmulti)
mmi <- as.data.frame(coef(res))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

#Checking model-averaged importance of terms
plot(res, type="s")

# Best model
MRbest<-rma.mv(-1*yi, vi, tdist=TRUE, random = list(~ 1 | Site, ~1|DisturbanceName), 
               mods=~soilP+hurrwind+long, data=hurr_sites)
MRbest

### Model R^2
(sum(full.model3$sigma2) - sum(MRbest$sigma2)) / sum(full.model3$sigma2)
#0.52

## Heterogeneity
mlm.variance.distribution(MRbest)

predict(MRbest, newmods=c(0,0,0), addx=TRUE)

### Using fitstats

model_P<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "REML",
                mods = ~soilP)
summary(model_P)
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_P$sigma2)) / sum(full.model3$sigma2)
#0.33

model_P.1<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "ML",
                mods = ~soilP)
summary(model_P.1)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_P.1$sigma2)) / sum(full.model3$sigma2)
#0.33

model_P.1.5<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                  data = hurr_sites,method = "ML",
                  mods = ~soilP + hurrwind)
summary(model_P.1.5)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_P.1.5$sigma2)) / sum(full.model3$sigma2)
#0.47

model_P.2<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "ML",
                mods = ~soilP + long)
summary(model_P.2)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_P.2$sigma2)) / sum(full.model3$sigma2)
#0.52

model_P.3<-rma.mv(-1*yi,vi,random = list(~ 1 | Site, ~1|DisturbanceName), tdist = TRUE, 
                data = hurr_sites,method = "ML",
                mods = ~soilP + long + hurrwind)
summary(model_P.3)#Both positive predictors, no significant interaction
### Model R^2 for each model
(sum(full.model3$sigma2) - sum(model_P.3$sigma2)) / sum(full.model3$sigma2)
#0.52

### compare fit statistics
fitstats(model_P.1, model_P.1.5,model_P.2,model_P.3)

#Best model, based on AICc is model_P.3, including soilP + long + hurrwind as moderators