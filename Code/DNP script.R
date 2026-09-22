#load packages
library(Hmisc)
library(tidyverse)
library(lme4)
library(cowplot)
library(lmerTest)
library(agricolae)
library(MASS)
library(emmeans)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(car)
library(sf)
library(ggjoy)
library(MuMIn)
library(plotrix)
library(HH)
library(bestNormalize)
library(svglite)
#Set theme globally
theme_set(theme_cowplot())

#####Assessing survival####

##Lethal dose analysis
datum = read.csv(file = "DNP_survival.csv")

datum2 = melt(datum[,c(2, 7, 8, 9, 10, 11)], id="ID")

datum2.5 <- as.data.frame(str_split_fixed(datum2$ID, "_", 2))
datum2$dose = as.factor(datum2.5$V1)
datum2$variable = as.numeric(datum2$variable)
datum2$variable = as.numeric(datum2$variable-1)

mod.1 = lm(value~dose*variable, data = datum2)
mod.1.comp <- emmeans(mod.1, pairwise~dose | variable)

datum2$dose <- factor(datum2$dose, levels = c("0uM", "0.5uM", "1uM", "2uM", "10uM", "25uM", "50uM", "100uM"))

jpeg(file = "survival curves.jpg", units = "in", width = 6.5, height = 4.5, res = 500)
ggplot(data=datum2, aes(x=variable, y=value, col = dose)) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", size =1, width=0.1) +
  stat_summary(fun = mean, geom = "line", size =1) +
  stat_summary(fun = mean, geom ="point", size = 3, shape =15, show.legend = FALSE)+
  geom_point()+
  scale_x_continuous(name="Day", breaks = seq(0, 4, by = 1))+
  scale_y_continuous(name="Percent survival", breaks = seq(0, 1, by = 0.1))+
  guides(color = guide_legend(title="Concentration DNP"))+
  scale_color_brewer(palette = "Accent")
 dev.off() 
 

 
############## DNP oxygen consumption  ##############################

####10uM 3 day trials
datum3 = read.csv(file = "InhibitorData10uM.csv")
 
 datum3$ID <- factor(datum3$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))

 
mod.2 = lm(abs(Slope_ppm_per_min)~Group * Sex + Weight_mg, data = datum3)
emmeans(mod.2, pairwise ~ Group | Sex)
confint(emmeans(mod.2, pairwise ~ Group | Sex))

#Ridgeline plot
resp.10uM.3days <- ggplot(datum3, aes(x=abs(Slope_ppm_per_min), y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.7, bandwidth = 0.3, rel_min_height =0.01,
   jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "red3", color = 'white')+
  scale_x_continuous(breaks = seq(0,6, by =1), limits = c(0,6)) +
  theme_ridges(center_axis_labels = TRUE)+
  xlab("") +ylab("")+
  ggtitle(label =  bquote('10'~mu*'M'~'DNP 3 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 5, y = 3.5, label = expression(italic(n)*"=36,"~italic(p)*"=0.0359"), size = 3)+
  annotate("text", x = 5, y = 1.5, label = expression(italic(n)*"=27,"~italic(p)*"=0.0011"), size = 3)

jpeg(file = "Resp_10uM_3day.jpg", units = "in", width = 5.5, height = 3, res = 300)
resp.10uM.3days
dev.off()


####10uM 7 day trials
datum4 = read.csv(file = "InhibitorData7Day10uM.csv")

datum4$ID <- factor(datum4$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))


mod.3 = lm(abs(Slope_ppm_per_min)~Group * Sex + Weight_mg, data = datum4)
emmeans(mod.3, pairwise ~ Group | Sex)

#Ridgeline plot
resp.10uM.7days <- ggplot(datum4, aes(x=abs(Slope_ppm_per_min), y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.2, bandwidth = 0.2, rel_min_height =0.01,
    jittered_points = TRUE, position = position_points_jitter(width = 0.00000101, height = 0), point_shape = "|", point_size = 5, point_color = "black",
    quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "red3", color = 'white')+
  #scale_x_continuous(breaks = seq(0,6, by =1), limits = c(0,6)) +
  theme_ridges(center_axis_labels = TRUE)+
  xlab("") +ylab("")+
  ggtitle(label =  bquote('10'~mu*'M'~'DNP 7 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 3, y = 3.5, label = expression(italic(n)*"=18,"~italic(p)*"=0.1416"), size = 3)+
  annotate("text", x = 3, y = 1.5, label = expression(italic(n)*"=18,"~italic(p)*"=0.9646"), size = 3)

jpeg(file = "Resp_10uM_7day.jpg", units = "in", width = 5.5, height = 3, res = 300)
resp.10uM.7days
dev.off()


###2uM analysis
datum5 = read.csv(file = "InhibitorData7Day2uM.csv")

datum5$ID <- factor(datum5$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))

mod.4 = lm(abs(Slope_ppm_per_min)~Group * Sex + Weight_mg, data = datum5)
emmeans(mod.4, pairwise ~ Group | Sex)
confint(emmeans(mod.4, pairwise ~ Group | Sex))

#Ridgeline plot

 resp.2uM.7day<- ggplot(datum5, aes(x=abs(Slope_ppm_per_min), y=ID, fill = ..x..)) +
    geom_density_ridges_gradient(scale=1.3, bandwidth = 0.2, rel_min_height =0.01,
     jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                                 quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "red3", color = 'white')+
   #scale_x_continuous(breaks = seq(0,6, by =1), limits = c(0,6)) +
    theme_ridges(center_axis_labels = TRUE)+
    xlab(bquote('Respiration rate (mmol O2'~min^-1*')')) +ylab("")+
    ggtitle(label =  bquote('2'~mu*'M'~'DNP 7 days - red stock copepods'))+
    theme(legend.position = "non")+
    annotate("text", x = 2, y = 3.5, label = expression(italic(n)*"=22,"~italic(p)*"=0.0346"), size = 3)+
    annotate("text", x = 2, y = 1.5, label = expression(italic(n)*"=28,"~italic(p)*"=0.0357"), size = 3)

  jpeg(file = "Resp_2uM_7day.jpg", units = "in", width = 5.5, height = 3, res = 300)
  resp.2uM.7day
  dev.off()



#####Yeast 2 um 7 day #####
datum5yeast = read.csv(file = "Yeast.data.csv")
str(datum5yeast)
#Remove erroneous data point
datum5yeast = datum5yeast[-c(18),]

#Make respiration data absolute
datum5yeast$resp.final = abs(datum5yeast$resp.final)

#order ID factor
datum5yeast$ID <- factor(datum5yeast$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))

#Run model and check results
mod.5 = lm(abs(Slope_mmol_per_min)~Group * Sex +Weight_mg, data = datum5yeast)
emmeans(mod.5, pairwise ~ Group | Sex)
confint(emmeans(mod.5, pairwise ~ Group | Sex))

#Ridgeline plot

resp.2uM.yeast <- ggplot(datum5yeast, aes(x=abs(Slope_mmol_per_min), y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.3, bandwidth = 0.1, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "red3", color = 'white')+
  theme_ridges(center_axis_labels = TRUE)+
  ggtitle(label = bquote('2'~mu*'M'~'DNP 7 days - color-restored copepods'))+
  scale_x_continuous(breaks = seq(0,3, by =0.5), limits = c(0,3)) + 
  xlab(bquote('Respiration rate (mmol O2'~min^-1*')')) +ylab("")+
  theme(legend.position = "non")+
  annotate("text", x = 2, y = 3.5, label = expression(italic(n)*"=54,"~italic(p)*"=0.0353"), size = 3)+
  annotate("text", x = 2, y = 1.5, label = expression(italic(n)*"=36,"~italic(p)*"=0.1298"), size = 3)

jpeg(file = "Resp_Yeast_7day_2uM.jpg", units = "in", width = 7, height = 5, res = 500)
resp.2uM.yeast
dev.off()



#Combine algae respiration figures into one three-panel figure
Fig.resp.ridgelines <- ggarrange(resp.10uM.3days, resp.10uM.7days, resp.2uM.7day,  
                                             labels = c("A", "B", "C"), 
                                             ncol=1, nrow=3, common.legend = FALSE, legend = "none")

jpeg(file="Algae resp ridgelines.jpg", units="in", width=6, height=10, res=500)
Fig.resp.ridgelines
dev.off()

#export as svg file
ggsave(file="Figure 2.svg", plot=Fig.resp.ridgelines, width=6, height=10)

############ HPLC data analysis ############

carotenoid.datum = read.csv(file = "carotenoid.data.csv")

carotenoid.datum$log.resp.min = log(abs(carotenoid.datum$resp.min))

#Split data set into separate data frames
#Plot data and run models for each group comparing across males and females
#algae fed copepods
carotenoid.datum.algae = subset(carotenoid.datum, diet == "algae" )
carotenoid.datum.algae = droplevels(carotenoid.datum.algae)

#Compare all identifiable carotenoids in algae dataset
#Subset desired samples from data frame
datum.algae_sub = carotenoid.datum.algae[,c(1,13,14,15)]
##Then rearrange your data frame
dd.algae = melt(datum.algae_sub, id=c("HPLC.ID"))

mod.carot.comparison <- lm(value~variable, data = dd.algae)
summary(mod.carot.comparison)
emmeans(mod.carot.comparison, pairwise ~ variable)

#Repeat for just beta-carot vs hydroxy
datum.algae_sub.2 = carotenoid.datum.algae[,c(1,14,15)]
##Then rearrange your data frame
dd.algae.2 = melt(datum.algae_sub.2, id=c("HPLC.ID"))

mod.carot.comparison.2 <- lm(value~variable, data = dd.algae.2)
summary(mod.carot.comparison.2)
confint(mod.carot.comparison.2)


#Ridgeline plot for all carots
all.carots.algae <- ggplot(dd.algae, aes(x=value, y=variable, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.3, bandwidth = 0.1, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,2.5, by =0.5), limits = c(0,3)) + 
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="yellow", high="red")+
  xlab(bquote('')) +ylab("")+
  ggtitle(label = "Red stock copepods")+
  scale_y_discrete(labels = c('Astaxanthin','Beta-carotene', 'Hydroxyechinenone'))+
  theme(legend.position = "non")

jpeg(file = "carotenoid comparison ridgeline.jpg", units = "in", width = 6, height = 4, res = 500)
all.carots.algae
dev.off()

#Ridgeline plot for just beta carot and hydroxy
dietary.carots.algae <- ggplot(dd.algae.2, aes(x=value, y=variable, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.3, bandwidth = 0.01, rel_min_height =0.01,
                               jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.2, by =0.05), limits = c(0,0.2)) + 
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="yellow", high="red")+
  xlab(bquote('')) +ylab("")+
  ggtitle(label = "Red stock copepods")+
  scale_y_discrete(labels = c('Beta-carotene', 'Hydroxyechinenone'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.15, y = 1.5, label = expression(italic(n)*"=162,"~italic(p)*"=0.0184"), size = 3)





####Yeast fed copepods 
carotenoid.datum.yeast = subset(carotenoid.datum, diet == "yeast" )
carotenoid.datum.yeast = droplevels(carotenoid.datum.yeast)
carotenoid.datum.yeast$ID <- factor(carotenoid.datum.yeast$ID, levels = c("Colorless Female", "Colorless Male","Control Female", "DNP Female", "Control Male", "DNP Male"))

#Compare all identifiable carotenoids in color restored copepods only (subset out colorless)
#Subset desired samples from data frame
datum.yeast_sub = carotenoid.datum.yeast[c(1:90),c(1,13,14,15)]
##Then rearrange your data frame
dd.yeast = melt(datum.yeast_sub, id=c("HPLC.ID"))

mod.carot.comparison.yeast <- lm(value~variable, data = dd.yeast)
summary(mod.carot.comparison.yeast)
emmeans(mod.carot.comparison.yeast, pairwise ~ variable)

#Repeat for just beta-carot and hydroxy
datum.yeast_sub.2 = carotenoid.datum.yeast[c(1:90),c(1,14,15)]
##Then rearrange your data frame
dd.yeast.2 = melt(datum.yeast_sub.2, id=c("HPLC.ID"))

mod.carot.comparison.yeast.2 <- lm(value~variable, data = dd.yeast.2)
summary(mod.carot.comparison.yeast.2)
confint(mod.carot.comparison.yeast.2)

#Ridgeline plot of all carotenoids
all.carots.yeast <- ggplot(dd.yeast, aes(x=value, y=variable, fill = ..x..)) +
  geom_density_ridges_gradient(scale=2, bandwidth = 0.1, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black", 
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,2.5, by =0.5), limits = c(0,3)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="yellow", high="red")+
  xlab(bquote('Concentration ('~mu*'g'~'mg'^-1*') across all samples')) +ylab("")+
  ggtitle(label = "Color-restored copepods")+
  scale_y_discrete(labels = c('Astaxanthin','Beta-carotene', 'Hydroxyechinenone'))+
  theme(legend.position = "non")

jpeg(file = "carotenoid comparison ridgeline yeast.jpg", units = "in", width = 6, height = 4, res = 500)
all.carots.yeast
dev.off()

#Ridgeline plot for just beta carot and hydroxy
dietary.carots.yeast <- ggplot(dd.yeast.2, aes(x=value, y=variable, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.3, bandwidth = 0.01, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.2, by =0.05), limits = c(0,0.2)) + 
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="yellow", high="red")+
  xlab(bquote('Concentration ('~mu*'g'~'mg'^-1*') across all samples')) +ylab("")+
  ggtitle(label = "Color-restored copepods")+
  scale_y_discrete(labels = c('Beta-carotene', 'Hydroxyechinenone'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.15, y = 1.5, label = expression(italic(n)*"=180,"~italic(p)*"<0.001"), size = 3)

#Combine yeast and algae carotenoids data into one figure
Fig.all.carots <- ggarrange(all.carots.algae,  dietary.carots.algae, all.carots.yeast, dietary.carots.yeast, 
                                         labels = c("A", "B", "C", "D"), 
                                         ncol=2, nrow=2, common.legend = FALSE, legend = "none")

jpeg(file="All carots.jpg", units="in", width=11.5, height=7, res=500)
Fig.all.carots
dev.off()



#Comparing astaxanthin concentration (subset out colorless copepods)
mod.yeast = lm(asta.conc~treatment * sex, data = carotenoid.datum.yeast[1:90,])
summary(mod.yeast)
emmeans(mod.yeast, pairwise ~ treatment | sex)
confint(emmeans(mod.yeast, pairwise ~ treatment | sex))

#Comparing astaxanthin concentration in colorless males and females
mod.yeast.b = lm(asta.conc~ID, data = carotenoid.datum.yeast[91:103,])
summary(mod.yeast.b)

#Compare astaxanthin concentration in colorless to control and DNP copepods 
mod.yeast.c = lm(asta.conc~ID * sex, data = carotenoid.datum.yeast)
emmeans(mod.yeast.c, pairwise ~ ID | sex)


#Ridgeline plot
asta.2uM.yeast <- ggplot(carotenoid.datum.yeast, aes(x=asta.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.7, bandwidth = 0.03, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,1, by =0.2), limits = c(0,1)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*'dry body weight)')) +ylab("")+
  ggtitle(label = bquote('2'~mu*'M'~'DNP 7 days - color-restored copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.88, y = 3.5, label = expression(italic(n)*"=36,"~italic(p)*"=0.0312"), size = 3)+
  annotate("text", x = 0.88, y = 5.5, label = expression(italic(n)*"=54,"~italic(p)*"=0.5836"), size = 3)+
  annotate("text", x = 0.35, y = 1.5, label = expression(italic(n)*"=12"), size = 3)

jpeg(file = "Asta_Yeast_7day_2uM.jpg", units = "in", width = 7, height = 5, res = 500)
asta.2uM.yeast
dev.off()



#Combine 2uM yeast respiration figure and asta figure into one
Fig.2uM.yeast.resp.and.asta <- ggarrange(resp.2uM.yeast, asta.2uM.yeast,  
                                 labels = c("A", "B"), 
                                 ncol=1, nrow=2, common.legend = FALSE, legend = "none")

jpeg(file="Resp and asta figure 2uM yeast.jpg", units="in", width=6, height=9, res=500)
Fig.2uM.yeast.resp.and.asta
dev.off()

#export as svg file
ggsave(file="Figure 4.svg", plot=Fig.2uM.yeast.resp.and.asta, width=6, height=9)



#Comparing dietary carotenoid concentration
mod.yeast.dietary = lm(dietary.conc~treatment * sex, data = carotenoid.datum.yeast)
summary(mod.yeast.dietary)
emmeans(mod.yeast.dietary, pairwise ~ treatment | sex)
#Ridgeline plot

dietary.2uM.yeast<- ggplot(carotenoid.datum.yeast[1:90,], aes(x=dietary.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.1, bandwidth = 0.003, rel_min_height =0.01,
   jittered_points = TRUE, position = position_points_jitter(width = 0.00000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.1, by =0.02), limits = c(0,0.1)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*'dry body weight)')) +ylab("")+
  ggtitle(label = bquote('2'~mu*'M'~'DNP 7 days - color-restored copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.08, y = 1.5, label = expression(italic(n)*"=36,"~italic(p)*"=0.5270"), size = 3)+
  annotate("text", x = 0.08, y = 3.5, label = expression(italic(n)*"=54,"~italic(p)*"=0.7660"), size = 3)

jpeg(file = "Dietary_Yeast_7day_2uM.jpg", units = "in", width = 7, height = 5, res = 500)
dietary.2uM.yeast
dev.off()





####10uM copepods for 3 days
carotenoid.datum.10uM = subset(carotenoid.datum.algae, conc.dnp == "10uM")
carotenoid.datum.10uM = droplevels(carotenoid.datum.10uM)
carotenoid.datum.10uM$ID <- factor(carotenoid.datum.10uM$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))

#Comparing astaxanthin concentration
mod.10uM = lm(asta.conc~treatment * sex , data = carotenoid.datum.10uM)
summary(mod.10uM)
emmeans(mod.10uM, pairwise ~ treatment | sex)
#Ridgeline plot

asta.10uM.3day <- ggplot(carotenoid.datum.10uM, aes(x=asta.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.1, bandwidth = 0.15, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,3, by =0.5), limits = c(0,3)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('')) +ylab("")+
  ggtitle(label = bquote('10'~mu*'M'~'DNP 3 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 2.5, y = 1.5, label = expression(italic(n)*"=22,"~italic(p)*"=0.5159"), size = 3)+
  annotate("text", x = 2.5, y = 3.5, label = expression(italic(n)*"=24,"~italic(p)*"=0.9071"), size = 3)


jpeg(file = "Asta_10uM_3day.jpg", units = "in", width = 7, height = 5, res = 500)
asta.10uM.3day
dev.off()

#Comparing dietary carotenoid concentration
mod.10uM.dietary = lm(dietary.conc~treatment * sex, data = carotenoid.datum.10uM)
summary(mod.10uM.dietary)
emmeans(mod.10uM.dietary, pairwise ~ treatment | sex)
#Ridgeline plot

dietary.10uM.3day <- ggplot(carotenoid.datum.10uM, aes(x=dietary.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.1, bandwidth = 0.009, rel_min_height =0.01,
   jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.2, by =0.05), limits = c(0,0.2)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('')) +ylab("")+
  ggtitle(label = bquote('10'~mu*'M'~'DNP 3 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.15, y = 1.5, label = expression(italic(n)*"=22,"~italic(p)*"=0.9594"), size = 3)+
  annotate("text", x = 0.15, y = 3.5, label = expression(italic(n)*"=24,"~italic(p)*"=0.5441"), size = 3)

jpeg(file = "Dietary_10uM_3day.jpg", units = "in", width = 7, height = 5, res = 500)
dietary.10uM.3day
dev.off()






#2uM copepods for 7 days
carotenoid.datum.2uM = subset(carotenoid.datum.algae, conc.dnp == "2uM")
carotenoid.datum.2uM = droplevels(carotenoid.datum.2uM)
carotenoid.datum.2uM$ID <- factor(carotenoid.datum.2uM$ID, levels = c("Control Female", "DNP Female", "Control Male", "DNP Male"))

#Comparing astaxanthin concentration
mod.2uM = lm(asta.conc~treatment * sex, data = carotenoid.datum.2uM)
summary(mod.2uM)
emmeans(mod.2uM, pairwise ~ treatment | sex)
confint(emmeans(mod.2uM, pairwise ~ treatment | sex))
#Ridgeline plot

asta.2uM.7days <- ggplot(carotenoid.datum.2uM, aes(x=asta.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.3, bandwidth = 0.15, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,3, by =0.5), limits = c(0,3)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*'dry body weight)')) +ylab("")+
  ggtitle(label = bquote('2'~mu*'M'~'DNP 7 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 2.5, y = 1.5, label = expression(italic(n)*"=20,"~italic(p)*"=0.3948"), size = 3)+
  annotate("text", x = 2.5, y = 3.5, label = expression(italic(n)*"=15,"~italic(p)*"=0.0042"), size = 3)

jpeg(file = "Asta_2uM_7day.jpg", units = "in", width = 7, height = 5, res = 500)
asta.2uM.7days
dev.off()

#Comparing dietary carotenoid concentration
mod.2uM.dietary = lm(dietary.conc~treatment * sex, data = carotenoid.datum.2uM)
summary(mod.2uM.dietary)
emmeans(mod.2uM.dietary, pairwise ~ treatment | sex)
#Ridgeline plot

dietary.2uM.7days <- ggplot(carotenoid.datum.2uM, aes(x=dietary.conc, y=ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.1, bandwidth = 0.009, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "blue", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.2, by =0.05), limits = c(0,0.2)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="red3")+
  xlab(bquote('Dietary carotenoids ('~mu*'g'~'mg'^-1*'dry body weight)')) +ylab("")+
  ggtitle(label = bquote('2'~mu*'M'~'DNP 7 days - red stock copepods'))+
  theme(legend.position = "non")+
  annotate("text", x = 0.15, y = 1.5, label = expression(italic(n)*"=20,"~italic(p)*"=0.8594"), size = 3)+
  annotate("text", x = 0.15, y = 3.5, label = expression(italic(n)*"=15,"~italic(p)*"=0.3012"), size = 3)

jpeg(file = "Dietary_2uM_7day.jpg", units = "in", width = 7, height = 5, res = 500)
dietary.2uM.7days
dev.off()




#Combine astaxanthin figures for 10uM 3day and 2uM 7day into one
Fig.asta.10uM.and.2uM <- ggarrange(asta.10uM.3day, asta.2uM.7days,  
                                         labels = c("A", "B"), 
                                         ncol=1, nrow=2, common.legend = FALSE, legend = "none")

jpeg(file="Asta 10uM and 2uM.jpg", units="in", width=6, height=8, res=500)
Fig.asta.10uM.and.2uM
dev.off()

#export as svg file
ggsave(file="Figure 3.svg", plot=Fig.asta.10uM.and.2uM, width=6, height=8)

#Combine dietary carotenoids figures for 10uM 3day and 2uM 7day into one
Fig.dietary.10uM.and.2uM <- ggarrange(dietary.10uM.3day, dietary.2uM.7days,  
                                         labels = c("A", "B"), 
                                         ncol=1, nrow=2, common.legend = FALSE, legend = "none")

jpeg(file="Dietary 10uM and 2uM.jpg", units="in", width=6, height=8, res=500)
Fig.dietary.10uM.and.2uM
dev.off()




####Comparing astaxanthin vs respiration controlling for effect of diet and sex####
#DNP copepods
mod.asta.vs.resp.DNP = lmer(asta.conc~abs(resp.min.mg) + sex + (1|diet), 
                            data = subset(carotenoid.datum, treatment == "DNP"), REML = FALSE)
summary(mod.asta.vs.resp.DNP)
confint(mod.asta.vs.resp.DNP)
r.squaredGLMM(mod.asta.vs.resp.DNP)

mod.asta.vs.resp.DNP2 = lm(asta.conc~abs(resp.min.mg) + sex + diet, 
                            data = subset(carotenoid.datum, treatment == "DNP"))
AIC(mod.asta.vs.resp.DNP2, mod.asta.vs.resp.DNP)
#plot asta vs resp
asta.vs.resp.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = abs(resp.min.mg), y=asta.conc)) +
  geom_point(aes(y = asta.conc),size = 5, alpha =0.9)+
  scale_x_continuous(breaks = seq(0, 80, by =20))+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "DNP treated")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  annotate("text", 70, 0.5, label =  bquote('p=0.029, R'^2*'=0.754'), size=4)

#Control copepods
mod.asta.vs.resp.control = lmer(asta.conc~abs(resp.min.mg)+ sex+(1|diet), 
                            data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.asta.vs.resp.control)
r.squaredGLMM(mod.asta.vs.resp.control)
#plot asta vs resp
asta.vs.resp.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = abs(resp.min.mg), y=asta.conc)) +
  geom_point(aes(y = asta.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mm'^-1*')'))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  annotate("text", 70, 0.5, label =  bquote('p=0.565, R'^2*'=0.768'), size=4)


###Comparing dietary carotenoids vs respiration controlling for effect of diet
xtrans <- bestNormalize(carotenoid.datum$dietary.conc)
carotenoid.datum$tf.dietary = xtrans$x.t

#DNP copepods
mod.dietary.vs.resp.DNP = lmer(tf.dietary~abs(resp.min.mg) + sex + (1|diet), 
                            data = subset(carotenoid.datum, treatment == "DNP"))
summary(mod.dietary.vs.resp.DNP)
confint(mod.dietary.vs.resp.DNP)
r.squaredGLMM(mod.dietary.vs.resp.DNP)
#plot asta vs resp
dietary.vs.resp.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = abs(resp.min.mg), y=tf.dietary)) +
  geom_point(aes(y = tf.dietary),size = 5, alpha =0.9)+
  scale_x_continuous(breaks = seq(0, 80, by =20))+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Dietary carotenoids concentration (transformed)'))  +
  ggtitle(label = "DNP treated")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 75, -1, label =  bquote('p=0.016, R'^2*'=0.134'), size=4)

#Control copepods
mod.dietary.vs.resp.control = lmer(tf.dietary~abs(resp.min.mg)+ sex+(1|diet), 
                                data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.dietary.vs.resp.control)
confint(mod.dietary.vs.resp.control)
r.squaredGLMM(mod.dietary.vs.resp.control)
#plot asta vs resp
dietary.vs.resp.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = abs(resp.min.mg), y=tf.dietary)) +
  geom_point(aes(y = tf.dietary),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Dietary carotenoids concentration (transformed)'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 70, -1, label =  bquote('p=0.608, R'^2*'=0.013'), size=4)


#Combine ALL scatterplots (astaxanthin and dietary) into one four-panel figure
Fig.carotenoids.vs.resp.overall <- ggarrange(asta.vs.resp.control, asta.vs.resp.dnp, dietary.vs.resp.control, dietary.vs.resp.dnp,  
                                             labels = c("A", "B", "C", "D"), 
                                         ncol=2, nrow=2, common.legend = TRUE, legend = "top")

jpeg(file="carotenoids.vs.resp.scatters.jpg", units="in", width=11, height=11, res=200)
Fig.carotenoids.vs.resp.overall
dev.off()

#export as svg file
ggsave(file="Figure 5.svg", plot=Fig.carotenoids.vs.resp.overall, width=11, height=11)

##HYDROXYECHINENONE
#DNP copepods
mod.hydroxy.vs.resp.DNP = lm(hydroxy.conc~abs(resp.min.mg) + sex + diet, 
                               data = subset(carotenoid.datum, treatment == "DNP"))
summary(mod.hydroxy.vs.resp.DNP)
confint(mod.hydroxy.vs.resp.DNP)
r.squaredGLMM(mod.hydroxy.vs.resp.DNP)
#plot asta vs resp
hydroxy.vs.resp.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = abs(resp.min.mg), y=hydroxy.conc)) +
  geom_point(aes(y = hydroxy.conc),size = 5, alpha =0.9)+
  scale_x_continuous(breaks = seq(0, 80, by =20))+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "DNP treated")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 75, 0.075, label =  bquote('p<0.001, R'^2*'=0.178'), size=4)

#Control copepods
mod.hydroxy.vs.resp.control = lm(hydroxy.conc~abs(resp.min.mg)+ sex+diet, 
                                   data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.hydroxy.vs.resp.control)
confint(mod.hydroxy.vs.resp.control)
r.squaredGLMM(mod.hydroxy.vs.resp.control)
#plot asta vs resp
hydroxy.vs.resp.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = abs(resp.min.mg), y=hydroxy.conc)) +
  geom_point(aes(y = hydroxy.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Dietary carotenoids concentration (transformed)'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 70, 0.06, label =  bquote('p=0.063, R'^2*'=0.088'), size=4)

##HYDROXYECHINENONE
#DNP copepods
mod.hydroxy.vs.resp.DNP = lm(hydroxy.conc~abs(resp.min.mg) + sex + diet, 
                             data = subset(carotenoid.datum, treatment == "DNP"))
summary(mod.hydroxy.vs.resp.DNP)
confint(mod.hydroxy.vs.resp.DNP)
r.squaredGLMM(mod.hydroxy.vs.resp.DNP)
#plot asta vs resp
hydroxy.vs.resp.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = abs(resp.min.mg), y=hydroxy.conc)) +
  geom_point(aes(y = hydroxy.conc),size = 5, alpha =0.9)+
  scale_x_continuous(breaks = seq(0, 80, by =20))+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "DNP treated")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 75, 0.075, label =  bquote('p<0.001, R'^2*'=0.178'), size=4)

#Control copepods
mod.hydroxy.vs.resp.control = lm(hydroxy.conc~abs(resp.min.mg)+ sex+diet, 
                                 data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.hydroxy.vs.resp.control)
confint(mod.hydroxy.vs.resp.control)
r.squaredGLMM(mod.hydroxy.vs.resp.control)
#plot asta vs resp
hydroxy.vs.resp.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = abs(resp.min.mg), y=hydroxy.conc)) +
  geom_point(aes(y = hydroxy.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Hydroxyechinenone ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 70, 0.06, label =  bquote('p=0.063, R'^2*'=0.088'), size=4)

##HYDROXYECHINENONE
#DNP copepods
mod.beta.vs.resp.DNP = lm(b.carot.conc~abs(resp.min.mg) + sex + diet, 
                             data = subset(carotenoid.datum, treatment == "DNP"))
summary(mod.beta.vs.resp.DNP)
confint(mod.beta.vs.resp.DNP)
r.squaredGLMM(mod.beta.vs.resp.DNP)
#plot asta vs resp
beta.vs.resp.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = abs(resp.min.mg), y=b.carot.conc)) +
  geom_point(aes(y = b.carot.conc),size = 5, alpha =0.9)+
  scale_x_continuous(breaks = seq(0, 80, by =20))+
  geom_smooth(method='lm')+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')'))+
  ylab(bquote('Beta-carotene ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "DNP treated")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 75, 0.06, label =  bquote('p=0.070, R'^2*'=0.243'), size=4)

#Control copepods
mod.beta.vs.resp.control = lm(b.carot.conc~abs(resp.min.mg)+ sex+diet, 
                                 data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.beta.vs.resp.control)
confint(mod.beta.vs.resp.control)
r.squaredGLMM(mod.beta.vs.resp.control)
#plot asta vs resp
beta.vs.resp.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = abs(resp.min.mg), y=b.carot.conc)) +
  geom_point(aes(y = b.carot.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Beta-carotene ('~mu*'g'~'mg'^-1*')'))+
  ylab(bquote('Dietary carotenoids concentration (transformed)'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 14))+
  annotate("text", 70, 0.06, label =  bquote('p=0.501, R'^2*'=0.112'), size=4)


#Combine ALL scatterplots (astaxanthin and dietary) into one four-panel figure
Fig.dietary.carotenoids.vs.resp.overall <- ggarrange(hydroxy.vs.resp.control, hydroxy.vs.resp.dnp, beta.vs.resp.control, beta.vs.resp.dnp,  
                                             labels = c("A", "B", "C", "D"), 
                                             ncol=2, nrow=2, common.legend = TRUE, legend = "top")

jpeg(file="dietary carotenoids.vs.resp.scatters.jpg", units="in", width=11, height=11, res=200)
Fig.dietary.carotenoids.vs.resp.overall
dev.off()




##Is there a relationship between astaxanthin and dietary carotenoids? 
#DNP copepods
mod.dietary.vs.asta = lmer(asta.conc~tf.dietary+ sex+(1|diet), 
                                   data = subset(carotenoid.datum, treatment == "DNP"))
summary(mod.dietary.vs.asta)
r.squaredGLMM(mod.dietary.vs.asta)

#plot asta vs dietary
dietary.vs.asta.dnp <- subset(carotenoid.datum, treatment == "DNP") %>%
  ggplot(aes(x = tf.dietary, y=asta.conc)) +
  geom_point(aes(y = asta.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Dietary carotenoids concentration (transformed'))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "DNP")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))+
  annotate("text", 2, 0, label =  bquote('p=0.894, R'^2*'=0.777'), size=4)

#Control copepods
mod.dietary.vs.asta.control = lmer(asta.conc~tf.dietary+ sex+(1|diet), 
                           data = subset(carotenoid.datum, treatment == "Control"))
summary(mod.dietary.vs.asta.control)
confint(mod.dietary.vs.asta.control)
r.squaredGLMM(mod.dietary.vs.asta.control)

#plot asta vs dietary
dietary.vs.asta.control <- subset(carotenoid.datum, treatment == "Control") %>%
  ggplot(aes(x = tf.dietary, y=asta.conc)) +
  geom_point(aes(y = asta.conc),size = 5, alpha =0.9)+
  geom_smooth(method='lm')+
  xlab(bquote('Dietary carotenoids concentration (transformed)'))+
  ylab(bquote('Astaxanthin ('~mu*'g'~'mg'^-1*')'))  +
  ggtitle(label = "Control")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))+
  annotate("text", 1.5, 0, label =  bquote('p=0.001, R'^2*'=0.774'), size=4)

#Combine scatterplots into one figure
Fig.dietary.vs.asta <- ggarrange(dietary.vs.asta.control, dietary.vs.asta.dnp,  labels = c("A", "B"), 
                                         ncol=2, nrow=1, common.legend = TRUE, legend = "top")

jpeg(file="dietary.vs.asta.scatters.jpg", units="in", width=11, height=6, res=500)
Fig.dietary.vs.asta
dev.off()

#export as svg file
ggsave(file="Figure 6.svg", plot=Fig.dietary.vs.asta, width=11, height=6)


#Is there a significant difference in mass between DNP-treated and control copepods?
mod.mass = lmer(mass~treatment + sex + (1|diet), data = carotenoid.datum[1:171,])
summary(mod.mass)
emmeans(mod.mass, pairwise~ treatment)

#Try with fixed effect of diet to see if difference between diets
mod.mass2 = lm(mass~treatment + sex + diet, data = carotenoid.datum[1:171,])
summary(mod.mass2)


#Ridgeline plot (subset out colorless data, irrelevant)
jpeg(file = "Mass_vs_treatment.jpg", units = "in", width = 7, height = 5, res = 500)
ggplot(carotenoid.datum[1:171,], aes(x=mass, y=treatment, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.1, bandwidth = 0.003, rel_min_height =0.01,
  jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,0.1, by =0.025), limits = c(0,0.1)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="grey", high="darkgreen")+
  xlab(bquote('Mass (mg dry body weight)')) +ylab("")+
  theme(legend.position = "non")+
  annotate("text", x = 0.085, y = 1.5, label = expression(italic(n)*"=171,"~italic(p)*"=0.945"), size = 3)
dev.off()


#Is there a significant difference in respiration rate between algae raised and yeast raised copepods? 
mod.resp.diet = lm(abs(resp.min)~diet + mass + sex, data = carotenoid.datum)
summary(mod.resp.diet)
emmeans(mod.resp.diet, pairwise~ treatment)


#Ridgeline plot (subset out colorless data, irrelevant)
jpeg(file = "resp.vs.diet.jpg", units = "in", width = 7, height = 5, res = 500)
ggplot(carotenoid.datum, aes(x=abs(resp.min.mg), y=diet, fill = ..x..)) +
  geom_density_ridges_gradient(scale=1.7, bandwidth = 5, rel_min_height =0.01,
                               jittered_points = TRUE, position = position_points_jitter(width = 0.000001, height = 0), point_shape = "|", point_size = 5, point_color = "black",
                               quantile_lines = TRUE, quantile_fun = mean, size =1, vline_color = "black", color = 'white')+
  scale_x_continuous(breaks = seq(0,100, by =20), limits = c(0,100)) +
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_gradient(low="lightblue", high="darkblue")+
  xlab(bquote('Respiration rate (mmol O2'~min^-1~'mg'^-1*')')) +ylab("")+
  theme(legend.position = "non")+
  annotate("text", x = 80, y = 1.8, label = expression(italic(n)*"=183,"~italic(p)*"<0.001"), size = 3.5)
dev.off()



####Chromatogram figure####
datumhplc = read.csv(file = "hplc.chromatograms.csv")

#Trim first 13 minutes off to remove solvent front peak and reduce white space of graph
#trim off last 5 minutes to remove equilibration period and reduce white space of graph
datumhplctrim = subset(datumhplc, time >=13 & time <=26)

#Subset desired samples from data frame
dd_sub = datumhplctrim[,c(1,3,4,7,8)]
dd_subtetra = datumhplctrim[,c(1,2,3,4)]

#You will notice numbers after the Y variable call for the standard mix. This is just an adjustment to 
#all of the values of the intensity column for that standard to correct for the baseline drift
#of the HPLC system

#Adjust y values for baseline drift during hplc
dd_sub$X3HE = dd_sub$X3HE+200
dd_subtetra$X3HE = dd_subtetra$X3HE+200

##Then rearrange your data frame
dd = melt(dd_sub, id=c("time"))
ddtetra = melt(dd_subtetra, id=c("time"))

#Make plot comparing DNP and control copepod chromatograms representatives
chromatogram.yeast <- ggplot(dd) + geom_line(aes(x=time, y=value, colour=variable)) +
  scale_colour_manual(name="",values=c("grey","goldenrod3", "red", "blue"),
                      labels=c("Standard Mix (Asta, zeax, canthax, beta-carot)","Hydroxyechinenone standard", "Control", "DNP"))+
  scale_y_continuous(breaks=seq(0,8000,1000))+
  scale_x_continuous(breaks=seq(15,27,3))+
  ylab("Intensity (mv)")+ xlab("Time (min)")+
  annotate("text", x = 15.7, y = 4800, label = "1",size = 4, color ="black")+
  annotate("text", x = 17.4, y = 5300, label = "2*", size = 4, color ="black")+
  annotate("text", x = 19, y = 4800, label = "3*", size = 4, color ="black")+
  annotate("text", x = 20.5, y = 8000, label = "4", size = 4, color ="black")+
  annotate("text", x = 23.9, y = 6800, label = "5",  size = 4, color ="black")+
  annotate("text", x = 24, y = 8000, label = "*not detected in any copepod samples",  size = 3, color ="black")+
  guides(colour = guide_legend(override.aes = list(size=2), nrow = 2, ncol = 2))+
  theme(legend.position = "top", legend.direction = "horizontal")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#Print plot
jpeg(filename = "representative chromatogram.jpg", units = "in", width = 7, height = 5, res=500)
chromatogram.yeast
dev.off()
