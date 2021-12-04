#Old Response Figure Codes

#Figure
preg_tot<-ggplot(metaregplot.1_x2, aes(x=hurr_sites.HURRECON_wind_ms, y=preds_x2.pred))+geom_point(shape=21,aes(col=fac_soilP,size=preds_x2.se),stroke=1.4)
preg_tot<-preg_tot+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_tot<-preg_tot+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#878787","#313695","#8073ac","#4575b4","#74add1","#abd9e9","#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027","#a50026"))
preg_tot<-preg_tot+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_sites.HURRECON_wind_ms,y=preds_x2$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
#geom_smooth(method = 'glm', col="#9f8a89",alpha=0.1,se=FALSE)
#preg_tot
#preg_tot<-preg_tot#+labs(x="Predicted response",y="Observed response")
#rsq_label <- paste('R^2 == 0.73')
preg_tot<-preg_tot+guides(size = guide_legend(override.aes = list(col = "black",shape=21)),color = guide_legend(override.aes = list(size = 8)))
#preg_tot
preg_tot<-preg_tot+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#878787","#313695","#8073ac","#4575b4","#74add1","#abd9e9","#00876c","#3c986d","#60a86d","#84b76e","#a8c671","#cdd476","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027","#a50026"))
Fig5a_new<-preg_tot+theme_pubr()+#scale_color_gradient(low="#544A4A", high="#D8001F")+#breaks=c(400,800,1200,1600,2000))+#scale_color_gradientn(colours = rainbow(5))+
  theme(legend.position = "right", legend.justification = "right",legend.title = element_text (size = 21), legend.text = element_text (size = 20),axis.title=element_text(size=28),axis.text=element_text(size=26))+
  labs(size="",color="Total soil P (mg/kg)")+guides(size=FALSE)+ annotate("text", x = 18, y = 6.5, label = "a Total litterfall", size=10,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5a_new####
Fig5a_new

ggsave(filename = "Fig5a_Final.png",
       plot = Fig5a_new, width = 12, height = 12, units = 'cm',
       scale = 2, dpi = 600)

#Leaf fall
preg_lf<-ggplot(metaregplot.1_x_leaf, aes(x=hurr_siteslf.HURRECON_wind_ms, y=preds_x_leaf.pred))+geom_point(shape=21,aes(col=fac_soilP,size=preds_x_leaf.se),stroke=1.4)
preg_lf<-preg_lf+ scale_size_continuous(range = c(2, 10))#+ scale_fill_fermenter(n.breaks = 12, palette = "Paired")#scale_fill_manual(values = c("#A38D8D","#F0A13C", "#B77878", "#C44474","#DE2BCF", "#CF1A33")) #+scale_fill_viridis_c(guide="legend")  #scale_fill_gradient2(low="black", mid = "white", high="red")+theme_pubr()#scale_size_area(max_size = 10)
preg_lf<-preg_lf+scale_color_manual(values=c("#1a1a1a","#4d4d4d","#8073ac","#4575b4","#abd9e9","#00876c","#84b76e","#a8c671","#f4e07f","#f4c76a","#f3ad5a","#e9774c","#de77ae","#c51b7d","#d73027"))
preg_lf
preg_lf<-preg_lf+ labs(x="HURRECON wind speed (m/s)", y="Predicted response to cyclone")+
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89")+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.lb),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)+
  stat_smooth(method="glm",aes(x=hurr_siteslf.HURRECON_wind_ms,y=preds_x_leaf$ci.ub),formula=y~x,fullrange=T,se=FALSE,size=1,colour="#9f8a89",linetype=3)
#geom_smooth(method = 'glm', col="#9f8a89",alpha=0.1,se=FALSE)
#preg_tot
#preg_tot<-preg_tot#+labs(x="Predicted response",y="Observed response")
#rsq_label <- paste('R^2 == 0.73')
preg_lf<-preg_lf+guides(size = guide_legend(override.aes = list(col = "black",shape=21)),color = guide_legend(override.aes = list(size = 8)))
preg_lf
Fig5b_new<-preg_lf+theme_pubr()+
  theme(legend.position = "top", legend.justification = "center",legend.title = element_text (size = 21), legend.text = element_text (size = 20),axis.title=element_text(size=28),axis.text=element_text(size=26))+
  labs(size="",color="Total soil P \n(mg/kg)")+guides(size=FALSE,color=FALSE)+ annotate("text", x = 16, y = 4.6, label = "b Leaf fall", size=10,hjust=0,colour="black",fontface="bold") #=  bquote('Density Litterfall N and P'~(mg/m^2/day))

#Fig5b_new####
Fig5b_new

Fig5_Response<-Fig5a_new+Fig5b_new+plot_layout(ncol=2)+ plot_layout(guides = 'collect')& theme(legend.justification = "left")
Fig5_Response

ggsave(filename = "Fig5ab_Final2.png",
       plot = Fig5_Response, width = 20, height = 12, units = 'cm',
       scale = 2, dpi = 600)

+scale_color_manual(values=c("#1dabe6","#b35a2d","#c3ced0","#ffa600","#665191","#af060f"))
plot2D<-plot1D+coord_flip()+geom_hline(aes(yintercept=0), lty=2, color = "magenta", cex=1.2, alpha = 0.6)#+geom_hline(aes(yintercept=3.64), lty=2, colour = "#003f5c", cex=0.9