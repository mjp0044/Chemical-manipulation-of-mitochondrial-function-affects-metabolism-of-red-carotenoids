###LH Data Analysis###
library(ggpubr) #For ggplot figures
library(cowplot) #clean up ggplot figures and set theme
theme_set(theme_cowplot())
library(lme4) #For glm models
library(MASS) #For glm.nb models
library(lmerTest) #To attach p-value reports to model summaries
library(blmeco) #For data handling with bayesian models and for checking dispersion of glmer models
library(Hmisc) #For data handling and manipulation for plotting
library(tidyverse) #For data handling and manipulation of data frames
library(MCMCglmm) #For bayesian MCMC models
library(HH) #contains the interval () function that calculates the 95% ci from glm 

#Read in data set and check variables
datum=read.csv(file="n1data.csv",na.strings = "na")
str(datum)

#Aggregate data measured as a single value per mother in a new data frame called "first"
first <- aggregate(firstfecun ~ momID, datum, mean)
total <- aggregate(totalfecun ~ momID, datum, mean)
average <- aggregate(avfecun ~ momID, datum, mean)
life <- aggregate(life ~ momID, datum, mean)
totclutch <- aggregate(totclutch~momID, datum, mean)
peak <- aggregate(peakfecun~momID, datum, mean)
surv <- aggregate(survavg~momID, datum, mean)
cope <- aggregate(copeavg~momID, datum, mean)
size <- aggregate(sizeavg~momID, datum, mean)
last <- aggregate(lastfecun~momID, datum, mean)
eggdevavg <- aggregate(eggdevavg~momID, datum, mean)

first$eggdevavg = eggdevavg$eggdevavg
first$last = last$lastfecun
first$peakfecun = peak$peakfecun
first$totalfecun = total$totalfecun
first$avfecun = average$avfecun
first$life = life$life
first$totclutch = totclutch$totclutch
first$survavg = surv$survavg
first$copeavg= cope$copeavg
first$sizeavg= size$sizeavg

write.csv(first, file = "n1data.avg.csv")

#Aggregate data measured as a single value per clutch in a new data frame called "second"
second <- aggregate(naup~clutchID, datum, mean)
eggdev <- aggregate(eggdev~clutchID, datum, mean)
momage <- aggregate(momage~clutchID, datum, mean)
clutchnum <- aggregate(clutchnum~clutchID, datum, mean)
ID <- aggregate(ID~clutchID, datum, mean)

second$eggdev = eggdev$eggdev
second$momage= momage$momage
second$clutchnum= clutchnum$clutchnum
second$ID= ID$ID
second$momID = substr(second$clutchID, 1, 4)

#Aggregate data based on clutch accounting for missing survival and development values
datum.missing = datum[-c(20,21,289,293,302,303,304,305,306,307,308,309,314,316),]
third <- aggregate(naup~clutchID, datum.missing, mean)
surv <- aggregate(totalsurv~clutchID, datum.missing, mean)
cope <- aggregate(copesurv~clutchID, datum.missing, mean)
survrat <- aggregate(survival~clutchID, datum.missing, mean)
coperat <- aggregate(development~clutchID, datum.missing, mean)
momage <- aggregate(momage~clutchID, datum.missing, mean)
clutchnum <- aggregate(clutchnum~clutchID, datum.missing, mean)
eggdev <- aggregate(eggdev~clutchID, datum.missing, mean)
life <- aggregate(life~clutchID, datum.missing, mean)

third$eggdev = eggdev$eggdev
third$survival = surv$totalsurv
third$copedev= cope$copesurv
third$survivalratio = survrat$survival
third$copedevratio= coperat$development
third$momage= momage$momage
third$clutchnum = clutchnum$clutchnum
third$momID = substr(third$clutchID, 1, 4)
third$life = life$life

third = droplevels(third)
str(third)

#Scale size data to curb warnings for large difference in variable scale in some models
datum$areaum2scale = scale(datum$areaum2)
#Create data frame with missing values in offspring size variable removed
datum.no.missing.size <- datum[!is.na(datum$areaum2),]
datum.no.missing.size = droplevels(datum.no.missing.size)
  
####Is the size of the first clutch a good proxy for fitness in Tigriopus, as used in other studies?####


####Analysis with lifespan as a covariate
#What is the correlation between the size of the first clutch and total fecundity?
res1 = glm.nb(totalfecun ~ firstfecun, data=first)
summary(res1)
exp(0.03973)
exp(confint(res1))

#Does first clutch size predict average fecundity? 
res1b = lm(avfecun ~ firstfecun, data=first)
summary(res1b)
confint(res1b)

#Does first clutch size predict peak clutch size? 
res1c = glm.nb(peakfecun~firstfecun, data=first)
summary(res1c)
exp(coef(res1c))
exp(confint(res1c))

#Does first clutch size predict life span? 
res1d = glm.nb(life~firstfecun, data=first)
summary(res1d)
exp(coef(res1d))
exp(confint(res1d))

#Does first clutch size predict total reproductive bouts?
res1e = glm.nb(totclutch ~ firstfecun, data=first)
summary(res1e) 
exp(coef(res1e))
exp(confint(res1e))

#Does first clutch size predict offspring survival, development rate and average body size?
#Survival average
res1f = glm(survavg ~ firstfecun, data=first, family = binomial)
summary(res1f)
exp(coef(res1f))
exp(confint(res1f))

#Development ratio
res1g = glm(copeavg ~ firstfecun , family = binomial, data=first)
summary(res1g)
exp(coef(res1g))
exp(confint(res1g))

#Offspring size
res1h = lm(sizeavg~firstfecun, data=first)
summary(res1h)
confint(res1h)

#Does first clutch size correlate with egg gestation period?
res1i = lm(eggdevavg~firstfecun, data=first)
summary(res1i)
confint(res1i)

## Plotting ##

#Total fecun P1
#Extract the predicted values and standard error  
pred = predict(res1,type = "response", se.fit = TRUE)
pred2=interval(res1, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds1= pred2[,1]
upr1 = pred2[,3]
lwr1 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p1 <- first %>% 
  mutate(pred = preds1,
         upr= upr1,
         lwr = lwr1) %>% 
  ggplot(aes(x=firstfecun, y=totalfecun))  + geom_point(aes(y=totalfecun), size=2)+
  xlab("First clutch size") + ylab("Total Offspring")+
  geom_line(aes(y=preds1), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr1,ymax=upr1), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=8.5), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.1), "cm"))+
    annotate("text", x = 7, y = 240, label = expression(italic(beta)*"=4.05,"~italic(p)*"=0.020"), 
         size = 3, color ="blue")


#Average fecun P2
p2 <- ggplot(first, aes(x=firstfecun, y=avfecun))  + geom_point(aes(y=avfecun), size=2)+
  xlab("First clutch size") + ylab("Avg clutch size")+
  geom_smooth(method='lm', color = "black", lwd = 1.5)+
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 7, y = 18, label = expression(italic(beta)*"=0.22,"~italic(p)*"=0.039"), 
           size = 3, color ="blue")


#Peak fecundity P3
#Extract the predicted values and standard error  
pred = predict(res1c,type = "response", se.fit = TRUE)
pred2=interval(res1c, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds3= pred2[,1]
upr3 = pred2[,3]
lwr3 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p3 <- first %>% 
  mutate(pred = preds3,
         upr= upr3,
         lwr = lwr3) %>% 
  ggplot(aes(x=firstfecun, y=peakfecun))  + geom_point(aes(y=peakfecun), size=2)+
  xlab("First clutch size") + ylab("Largest clutch size")+
  geom_line(aes(y=preds3), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr3,ymax=upr3), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0,0.2,0.2,0.5), "cm"))+
  annotate("text", x = 5.5, y = 38, label = expression(italic(beta)*"=1.67,"~italic(p)*"=0.072"), 
           size = 3.5, color ="blue")


#Lifespan P4
pred = predict(res1d,type = "response", se.fit = TRUE)
pred2=interval(res1d, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds4= pred2[,1]
upr4 = pred2[,3]
lwr4 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p4 <- first %>% 
  mutate(pred = preds4,
         upr= upr4,
         lwr = lwr4) %>% 
  ggplot(aes(x=firstfecun, y=life))  + geom_point(aes(y=life), size=2)+
  xlab("First clutch size") + ylab("Adult lifespan (days)")+
  geom_line(aes(y=preds4), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr4,ymax=upr4), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 7, y = 85, label = expression(italic(beta)*"=2.79,"~italic(p)*"=0.020"), 
           size = 3, color ="blue")


#Total number of clutches P5
pred = predict(res1e,type = "response", se.fit = TRUE)
pred2=interval(res1e, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds5= pred2[,1]
upr5 = pred2[,3]
lwr5 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p5 <- first %>% 
  mutate(pred = preds5,
         upr= upr5,
         lwr = lwr5) %>% 
  ggplot(aes(x=firstfecun, y=totclutch))  + geom_point(aes(y=totclutch), size=2)+
  xlab("First clutch size") + ylab("Total Clutches")+
  geom_line(aes(y=preds5), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr5,ymax=upr5), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0,0,0.2,0.5), "cm"))+
  annotate("text", x = 5.5, y = 15.3, label = expression(italic(beta)*"=2.20,"~italic(p)*"=0.092"), 
           size = 3.5, color ="blue")


#Survival average P6
pred = predict(res1f,type = "response", se.fit = TRUE)
pred2=interval(res1f, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds6= pred2[,1]
upr6 = pred2[,3]
lwr6 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p6 <- first %>% 
  mutate(pred = preds6,
         upr= upr6,
         lwr = lwr6) %>% 
  ggplot(aes(x=firstfecun, y=survavg))  + geom_point(aes(y=survavg), size=2)+
  xlab("First clutch size") + ylab("Avg survival ratio")+
  geom_line(aes(y=preds6), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr6,ymax=upr6), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_blank(),
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 5.5, y = 1.1, label = expression(italic(beta)*"=1.71,"~italic(p)*"=0.774"), 
           size = 3.5, color ="blue")


#Development average P7
pred = predict(res1g,type = "response", se.fit = TRUE)
pred2=interval(res1g, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
preds7= pred2[,1]
upr7 = pred2[,3]
lwr7 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p7 <- first %>% 
  mutate(pred = preds7,
         upr= upr7,
         lwr = lwr7) %>% 
  ggplot(aes(x=firstfecun, y=copeavg))  + geom_point(aes(y=copeavg), size=2)+
  xlab("First clutch size") + ylab("Avg development ratio")+
  geom_line(aes(y=preds7), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwr7,ymax=upr7), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_blank(),
        plot.margin = unit(c(0.5,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 5.5, y = 1.1, label = expression(italic(beta)*"=1.84,"~italic(p)*"=0.757"), 
           size = 3.5, color ="blue")


#Average offspring body size P8
p8 <- ggplot(first, aes(x=firstfecun, y=sizeavg))  + geom_point(aes(y=sizeavg), size=2)+
  xlab("First clutch size") + ylab(expression(paste("Avg offspring body size (", uM^2,")")))+
  geom_smooth(method='lm', color = "black", lwd = 1.5)+
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=8.5),
        axis.title.y = element_text(size=10.5), axis.title.x = element_blank(),
        plot.margin = unit(c(0.5,0,0.2,0.2), "cm"))+
  annotate("text", x = 6, y = 15000, label = expression(italic(beta)*"=13.25,"~italic(p)*"=0.738"), 
           size = 3.5, color ="blue")

#Average egg gestation time P9
p9 <- ggplot(first, aes(x=firstfecun, y=eggdevavg))  + geom_point(aes(y=eggdevavg), size=2)+
  xlab("First clutch size") + ylab("Avg egg gestation (days)")+ylim(2.4,4)+
  geom_smooth(method='lm', color = "black", lwd = 1.5)+
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=10.5), axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0.0,0.5), "cm"))+
  annotate("text", x = 5.5, y = 4, label = expression(italic(beta)*"=0.03,"~italic(p)*"=0.166"), 
           size = 3.5, color ="blue")

#Significant results
figure1 <- ggarrange(p4, p1, p2, labels = c("A", "B", "C"),  ncol=1, nrow=3, hjust = -0.2, vjust = 1.1)

jpeg(file="Fig1.jpg", units="in", width=3, height=7, res=300)
figure1
dev.off()

#Non-significant results
figureS3 <- ggarrange(p8, p7, p9, p6, p5, p3, labels = c("A", "B", "C", "D", "E", "F"),  ncol=2, nrow=3,
                     hjust = -0.2, vjust = 1.1)

jpeg(file= "FigS3.jpg", units = "in", width = 6.5, height =7.5, res=300)
figureS3
dev.off()



####Run simulation to test likelihood that relationship between first clutch and total output is signficant####
####
#LRT with Null Simulation outline

## (1)  Create empty column in mother averaged data set to hold our simulated total offspring variable
## (2)  Define simulation function 
## (2b) Simulate new first clutch size from poisson distribution using mean scalar from original data
## (2c) Independently generate remaining reproductive output that adds to total offspring
## (2d) run the null simulated model testing effect of simulated first clutch on simulated total offspring
## (2e) Extract effect size and store in temporary list
## (3)  Run original model to get observed, original effect size
## (3b) Run simulation 1000 times to generate distribution of effect sizes under the null hypothesis
## (4)  Plot distribution of effect sizes from simulation with line marking original estimate 
## (4b) Calculate the associated p-value from the LRT of the two models (null simulation vs observed data model)
####

#Set the sim seed
set.seed(999)

#### Begin Define Simualtion Function for Null Hypothesis

## (1) Set up column in data set to hold simulated total offspring values under the Null 
first$Fsim=NA
## (2) Define simulation function
firstclutchsim = function() {
  
  #Retrieve mean first clutch size from dataset
  FCF=mean(first$firstfecun)
  #Retrieve mean total offspring from data set
  F=mean(first$totalfecun)
  #Take the difference of the observed data
  Diff=F-FCF
  
  # (2b) Randomly generate a new set of first clutch sizes using the original mean
  FCFsim=rpois(23, FCF)
  #Generate a set of differences between total offspring-first clutch from mean diff 
  # (2c) Independent of values of simualated first clutch (f1) set, isolating direct effect of f1
  Diffsim=rpois(23, Diff)
  
  #Create new simulated response variable under null (only direct effect of first clutch being added)
  first$Fsim=FCFsim+Diffsim
  #Add simulated first clutch sizes to data frame
  first$simf1 = FCFsim
  
  ##(2d) Run the model with the simulated null response variable against the simulated first clutch size
  simmodel=glm.nb(Fsim~simf1, data=first)
  
  ## (2e) store effect size from model in a temporary variable called tmp_beta
  tmp_beta=summary(simmodel)$coefficients["simf1","Estimate"]
  tmp_beta
}

#################################### End Define Simulation Function 
## (3) Get observed beta from best AIC model     
#model with least number of variables
with=glm.nb(totalfecun~firstfecun, data=first)
summary(with)
#Get original effect size estimate from glm model
obs_betaglm=summary(with)$coefficients["firstfecun","Estimate"]

#Test once to see if simulation works, ignore sporadic warnings for theta going to infinity. 
#Happens approximately every three iterations. Product of fitting negative binomial model to data that's not that overdispersed.
firstclutchsim()

#### (3b) Run 1000 times to get distribution of betas
distobetas=replicate(1000,firstclutchsim())
obs_betaglm

################################ Compare simulated betas vs observed beta 
##(4) plot distribution of effect sizes and add a line marking observed effect size from original model
#If blue line falls outside simulated null model distribution, 
#there is an added effect between first clutch and total offspring that cannot be explained due to simple addition
#histogram of simulated effect sizes from null model
tiff(file = "Figure S2 First clutch Hist LRT.tiff", units = "in", width = 5, height = 5.5, res = 300)
hist(distobetas, breaks=20, xlim = c(-0.06,0.06), main = "", xlab = "Distribution of null effect sizes")
text(0.05, 100, "Observed 
effect", col = "blue", cex=0.7)
#add line marking original effect size
abline(v = obs_betaglm, col="blue", lwd=3)
dev.off()

## (4b) Calculate the p-value from LRT of effect sizesf
#Likelihood that observed effect is significantly different from simulated effect. 
# Calculate the number of betas under the null model greater than or equal to the observed
moreEx <- which(distobetas >= obs_betaglm)
# Calculate proportion of cases that were more extreme to get p-value
length(moreEx) / 1000  
#p-value extremely small means that observed effect never falls within simulated effect distribution

#### What is the effect of aging and increasing reproductive bouts on offspring? ####

#single random effect prior, weak prior chosen to maximize influence of data instead of assumptions
prior=list(R=list(V=diag(2), nu=1), G=list(G1=list(V=diag(2), nu=1)))

#1 What is the effect of mom age on clutch size?
res2 = glmer.nb(naup~ momage + (1|momID), data=second)
summary(res2)
exp(-0.018699)
exp(confint.merMod(res2))
#Check overdispersion
dispersion_glmer(res2)

#Estimate covariances among and between mothers
m2 <- MCMCglmm(cbind(naup, momage) ~ trait-1,
               random=~us(trait):momID,
               rcov = ~ us(trait):units, 
               prior = prior,
               family =  c("poisson", "poisson"), nitt = 100000, burnin = 20000,
               thin = 25, data = second)
summary(m2)

#2 What is the effect of clutch number on clutch size?
res2b = glmer.nb(naup~clutchnum + (1|momID), data=second, family=poisson)
summary(res2b)
exp(-0.05318)
exp(confint.merMod(res2b))
dispersion_glmer(res2b)

#Estimate covariances among and between mothers
m2b <- MCMCglmm(cbind(naup, clutchnum) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("poisson", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = second)
summary(m2b)

#3 What is the effect of age on offspring size?
res2c = lmer(areaum2 ~ momage + (1|momID),  data=datum.no.missing.size)
summary(res2c)
confint(res2c)

#Estimate covariances among and between mothers
m2c <- MCMCglmm(cbind(areaum2, momage) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = datum.no.missing.size)
summary(m2c)

#4What is the effect of clutch number on offspring size?
res2d = lmer(areaum2~clutchnum + (1|momID), data = datum.no.missing.size)
summary(res2d)
confint(res2d)

#Estimate covariances among and between mothers
m2d <- MCMCglmm(cbind(areaum2, clutchnum) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = datum.no.missing.size)
summary(m2d)

#5What is the effect of age on offspring survival?
res2e = glmer(survivalratio~momage + (1|momID), data = third, family = binomial)
summary(res2e)
exp(-0.02032)
exp(confint.merMod(res2e))

#Estimate covariances among and between mothers
m2e <- MCMCglmm(cbind(survivalratio, momage) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m2e)

#6What is the effect of clutch number on offspring survival?
res2f = glmer(survivalratio~clutchnum + (1|momID), data = third, family = binomial)
summary(res2f)
exp(-0.03047)
exp(confint.merMod(res2f))

#Estimate covariances among and between mothers
m2f <- MCMCglmm(cbind(survivalratio, clutchnum) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m2f)

#7What is the effect of age on offspring development?
res2g = glmer(copedevratio~momage + (1|momID), data = third, family = binomial)
summary(res2g)
exp(-0.02269)
exp(confint.merMod(res2g))

#Estimate covariances among and between mothers
m2g <- MCMCglmm(cbind(copedevratio, momage) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m2g)

#8What is the effect of clutch number on offspring development?
res2h= glmer(copedevratio~clutchnum + (1|momID), data = third, family = binomial)
summary(res2h)
exp(-0.05284)
exp(confint.merMod(res2h))

#Estimate covariances among and between mothers
m2h <- MCMCglmm(cbind(copedevratio, clutchnum) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m2h)

#9What is the effect of age on egg development time?
res2i = glmer(eggdev~momage + (1|momID) + (0+momage|ID), data = second, family = poisson)
summary(res2i)
exp(0.007807)
exp(confint.merMod(res2i))
dispersion_glmer(res2i)

#Estimate covariances among and between mothers
m2i <- MCMCglmm(cbind(eggdev, momage) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("poisson", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = second)
summary(m2i)

#10What is the effect of clutch number on egg development time?
res2j = glmer(eggdev~clutchnum+ (1|momID), data = second, family = poisson)
summary(res2j)
exp(0.01464)
exp(confint.merMod(res2j))
dispersion_glmer(res2j)

#Estimate covariances among and between mothers
m2j <- MCMCglmm(cbind(eggdev, clutchnum) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("poisson", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = second)
summary(m2j)

#### Plot effects of aging and increasing reproductive bouts####  
jpeg(file = "Fig2.jpeg", units = "in", height = 8.5, width = 7, res = 300)
#1Mom age vs clutch size (sig)
par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(4,5,0.1,0.1), cex.lab=1.6)
plot(naup~momage,data=second, col="grey", pch=15.5, xlab="",ylab="Clutch size (# Indv.)",ylim=c(0,40),xlim=c(0,60), bty="l")
### Go through each mom ID:
for(i in unique(as.character(second$momID))){
  ### Subset data by mom
  sub_naup <- subset(second, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$momage), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$momage,
        y = predict(res2, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- second[order(second$momage), ]
lines(x = ordersecond$momage, y = predict(glm(naup~momage, data=second, family=poisson), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(49, 31, expression(italic(beta)*"=-1.85,"~italic(p)*"<0.001"), col="blue", cex=1.2)
mtext("A", side=2, adj=4, padj = -9, las=2)

#2Clutch number vs clutch size (sig)
plot(naup~clutchnum,data=second, col="grey", pch=15.5, xlab="",ylab="",ylim=c(0,40),xlim=c(0,16), bty="l")
### Go through each mom ID:
for(i in unique(as.character(second$momID))){
  ### Subset data by mom
  sub_naup <- subset(second, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$clutchnum), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$clutchnum,
        y = predict(res2b, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- second[order(second$clutchnum), ]
lines(x = ordersecond$clutchnum, y = predict(glm(naup~clutchnum, data=second, family=poisson), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(13, 31, expression(italic(beta)*"=-5.18,"~italic(p)*"=0.011"), col="blue", cex=1.2)
mtext("B", side=2, adj=4, padj = -9, las=2)

#3mom age vs offspring body size (sig)
plot(areaum2~momage,data=datum, col="grey", pch=15.5, xlab="",ylab= expression(paste("Offspring body size (uM" ^"  2", ")")),ylim=c(5000,22000),xlim=c(0,60), bty="l")
### Go through each mom ID:
for(i in unique(as.character(datum$momID))){
  ### Subset data by mom
  sub_naup <- subset(datum, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$momage), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$momage,
        y = predict(res2c, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- datum[order(datum$momage), ]
lines(x = ordersecond$momage, y = predict(lm(areaum2~momage, data=datum, na.action = na.exclude), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(49, 21000, expression(italic(beta)*"=-36.64,"~italic(p)*"=0.004"), col="blue", cex=1.2)
mtext("C", side=2, adj=4, padj = -9, las=2)

#4clutch number vs offspring body size (sig)
plot(areaum2~clutchnum,data=datum, col="grey", pch=15.5, xlab="",ylab="",ylim=c(5000,22000),xlim=c(0,16), bty="l")
### Go through each mom ID:
for(i in unique(as.character(datum$momID))){
  ### Subset data by mom
  sub_naup <- subset(datum, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$clutchnum), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$clutchnum,
        y = predict(res2d, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- datum[order(datum$clutchnum), ]
lines(x = ordersecond$clutchnum, y = predict(lm(areaum2~clutchnum, data=datum, na.action = na.exclude), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(13, 21000, expression(italic(beta)*"=-131.69,"~italic(p)*"=0.007"), col="blue", cex=1.2)
mtext("D", side=2, adj=4, padj = -9, las=2)

#9Mom age vs egg development time
res2iplot = glmer(eggdev~momage + (1|momID) + (0+momage|ID), data = second, family = poisson)
plot(eggdev~momage,data=second, col="grey", pch=15, xlab="Mother age (days)",ylab="Egg gestation time (days)",ylim=c(0,10),xlim=c(0,60), bty="l")
### Go through each mom ID:
for(i in unique(as.character(second$momID))){
  ### Subset data by mom
  sub_naup <- subset(second, momID == i)
  #### Order subsetted data by increasing mom age
  sub_naup <- sub_naup[order(sub_naup$momage), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$momage,
        y = predict(res2iplot, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- second[order(second$momage), ]
lines(x = ordersecond$momage, y = predict(glm(eggdev~momage, data=second, family = poisson), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(49, 8, expression(italic(beta)*"=0.78,"~italic(p)*"=0.007"), col="blue", cex=1.2)
mtext("E", side=2, adj=4, padj = -9, las=2)

#10Clutch number vs egg development time
res2jplot = glmer(eggdev~clutchnum + (1|momID) + (0+clutchnum|ID), data = second, family = poisson)
plot(eggdev~clutchnum,data=second, col="grey", pch=15, xlab="Clutch number",ylab="",ylim=c(0,10),xlim=c(0,16), bty="l")
### Go through each mom ID:
for(i in unique(as.character(second$momID))){
  ### Subset data by mom
  sub_naup <- subset(second, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$clutchnum), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$clutchnum,
        y = predict(res2jplot, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- second[order(second$clutchnum), ]
lines(x = ordersecond$clutchnum, y = predict(glm(eggdev~clutchnum, data=second, family = poisson), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(13, 8, expression(italic(beta)*"=1.47,"~italic(p)*"=0.244"), col="blue", cex=1.2)
mtext("F", side=2, adj=4, padj = -9, las=2)
par(xpd=NA)
legend(-8, -1.3, legend = c("Individual females", "Mean response"), cex=1.4, bty="n",
       pch = c(NA, NA), lty = c(1, 2),
       col = c("darkgrey","black"), text.col = c("darkgrey", "black"))
dev.off()

jpeg(file="FigS4.jpg", units = "in", height = 6, width =7, res=300)
par(mfrow=c(2,2), oma=c(1,1,1,1), mar=c(4,4.5,0.1,0.1), cex.lab=1)
#7Mom age vs copepod development ratio
plot(copedevratio~momage,data=third, col="grey", pch=16, xlab="",ylab="Development ratio",ylim=c(0,1),xlim=c(0,60), bty="l")
### Go through each mom ID:
for(i in unique(as.character(third$momID))){
  ### Subset data by mom
  sub_naup <- subset(third, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$momage), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$momage,
        y = predict(res2g, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- third[order(third$momage), ]
lines(x = ordersecond$momage, y = predict(glm(copedevratio~momage, data=third, family = binomial), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(47, 0.1, expression(italic(beta)*"=-2.24,"~italic(p)*"=0.062"), col="blue", cex=1)
mtext("A",  side=2, adj=4, padj = -9, las=2)

#8Clutch number vs copepod development ratio
plot(copedevratio~clutchnum,data=third, col="grey", pch=16, xlab="",ylab="",ylim=c(0,1),xlim=c(0,16), bty="l")
### Go through each mom ID:
for(i in unique(as.character(third$momID))){
  ### Subset data by mom
  sub_naup <- subset(third, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$clutchnum), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$clutchnum,
        y = predict(res2h, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- third[order(third$clutchnum), ]
lines(x = ordersecond$clutchnum, y = predict(glm(copedevratio~clutchnum, data=third, family = binomial), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(13, 0.1, expression(italic(beta)*"=-5.15,"~italic(p)*"=0.331"), col="blue", cex=1)
mtext("B",  side=2, adj=4, padj = -9, las=2)

#5Mother age vs offspring survival ratio
plot(survivalratio~momage,data=third, col="grey", pch=16, xlab="Mother age (days)",ylab="Survival proportion",ylim=c(0,1),xlim=c(0,60), bty="l")
### Go through each mom ID:
for(i in unique(as.character(third$momID))){
  ### Subset data by mom
  sub_naup <- subset(third, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$momage), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$momage,
        y = predict(res2e, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- third[order(third$momage), ]
lines(x = ordersecond$momage, y = predict(glm(survivalratio~momage, data=third, family = binomial), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(47, 0.1, expression(italic(beta)*"=-2.01,"~italic(p)*"=0.091"), col="blue", cex=1)
mtext("C",  side=2, adj=4, padj = -9, las=2)

#6Clutch number vs offspring survival ratio
plot(survivalratio~clutchnum,data=third, col="grey", pch=16, xlab="Clutch number",ylab="",ylim=c(0,1),xlim=c(0,16), bty="l")
### Go through each mom ID:
for(i in unique(as.character(third$momID))){
  ### Subset data by mom
  sub_naup <- subset(third, momID == i)
  #### Order subsetted data by increasing clutch number
  sub_naup <- sub_naup[order(sub_naup$clutchnum), ]
  ### Plot model prediction with type set to response
  lines(x = sub_naup$clutchnum,
        y = predict(res2f, type="response",
                    newdata = sub_naup),
        col = "grey", lwd = 1.5, lty=1)
}
#### Add mean regression line
ordersecond <- third[order(third$clutchnum), ]
lines(x = ordersecond$clutchnum, y = predict(glm(survivalratio~clutchnum, data=third, family = binomial), type="response", newdata = ordersecond), lwd = 5, lty=2)
text(13, 0.1, expression(italic(beta)*"=-3.01,"~italic(p)*"=0.586"), col="blue", cex=1)
mtext("D",  side=2, adj=4, padj = -9, las=2)
par(xpd=NA)
legend(-9, -0.2, legend = c("Individual females", "Mean response"), cex=1.1, bty="n",
       pch = c(NA, NA), lty = c(1, 2),
       col = c("darkgrey","black"), text.col = c("darkgrey", "black"))
dev.off()

####Look for tradeoffs in reproductive variables using ####

#1Trade-off between clutch size and offspring size (control for mom age) (sig)
res3 = glmer.nb(naup ~ areaum2scale + momage+ (1|momID), data=datum.no.missing.size)
summary(res3)
exp(-0.102265)
exp(-0.102265-1.96*0.040626)
exp(-0.102265+1.96*0.040626)

#Estimate covariances among and between mothers
m3 <- MCMCglmm(cbind(naup, areaum2) ~ trait-1,
               random=~us(trait):momID,
               rcov = ~ us(trait):units, 
               prior = prior,
               family =  c("poisson", "gaussian"), nitt = 100000, burnin = 20000,
               thin = 25, data = datum.no.missing.size)
summary(m3)

#2 Relationship between offspring size and survival (control for mom age)? (non-sig)
#Create data frame with missing values in offspring size and survival variables removed
datum.no.missing.size.surv <- datum.no.missing.size[!is.na(datum.no.missing.size$survival),]

res3b = glmer(survival ~ areaum2scale + momage + (1|momID), data=datum.no.missing.size.surv, family = binomial)
summary(res3b)
exp(0.11280)
exp(0.11280-1.96*0.20103)
exp(0.11280+1.96*0.20103)


#Estimate covariances among and between mothers
m3b <- MCMCglmm(cbind(survival, areaum2scale) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "gaussian"), nitt = 100000, burnin = 20000,
                thin = 25, data = datum.no.missing.size.surv)
summary(m3b)

#3 Relationship betweeen clutch size and survival (non-sig)
res3c = glmer(survivalratio ~ naup +momage +(1|momID), data=third, family = binomial)
summary(res3c)
exp(0.04522)
exp(0.04522-1.96*0.02374)
exp(0.04522+1.96*0.02374)

#Estimate covariances among and between mothers
m3c <- MCMCglmm(cbind(survivalratio, naup) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m3c)

#4 Relationship between offspring size and development ratio (control for mom age) (non-sig)
#use scaled size data again
#Create data frame with missing values in offspring size and development variables removed
datum.no.missing.size.dev <- datum.no.missing.size[!is.na(datum.no.missing.size$development),]

res3d = glmer(development ~ areaum2scale + momage+ (1|momID), data=datum.no.missing.size.dev, family = binomial)
summary(res3d)
exp(0.02846)
exp(0.02846-1.96*0.18203)
exp(0.02846+1.96*0.18203)

#Estimate covariances among and between mothers
m3d <- MCMCglmm(cbind(development, areaum2scale) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "gaussian"), nitt = 100000, burnin = 20000,
                thin = 25, data = datum.no.missing.size.dev)
summary(m3d)

#5 Relationship between clutch size and development ratio (control for mom age) (non-sig)
res3e = glmer(copedevratio ~ naup +momage+ (1|momID), data=third, family = binomial, na.action = na.exclude)
summary(res3e)
exp(0.03392)
exp(0.03392-1.96*0.02271)
exp(0.03392+1.96*0.02271)

#Estimate covariances among and between mothers
m3e <- MCMCglmm(cbind(copedevratio, naup) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m3e)

#6relationship between egg dev and size (control for mom age) (sig)
res3f = lmer(areaum2 ~ eggdev +momage+(1|momID), data=datum.no.missing.size)
summary(res3f)
confint(res3f)

#Estimate covariances among and between mothers
m3f <- MCMCglmm(cbind(areaum2, eggdev) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = datum.no.missing.size)
summary(m3f)

#7relationship between egg dev and clutch size (control for mom age) (sig)
res3g = glmer.nb(naup ~ eggdev +momage+(1|momID), data=second, family=poisson)
summary(res3g)
exp(-0.111763)
exp(-0.111763-1.96*0.037794)
exp(-0.111763+1.96*0.037794)
dispersion_glmer(res3g)

#Estimate covariances among and between mothers
m3g <- MCMCglmm(cbind(naup, eggdev) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("poisson", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = second)
summary(m3g)

#8relationship between egg dev and survival ratio (sig)
res3h = glmer(survivalratio ~ eggdev +momage+ +(1|momID), data=third, family=binomial)
summary(res3h)
exp(-0.30962)
exp(-0.30962-1.96*0.13858)
exp(-0.30962+1.96*0.13858)

#Estimate covariances among and between mothers
m3h <- MCMCglmm(cbind(survivalratio, eggdev) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m3h)

#9relationship between egg dev and development ratio (control for mom age) (nonsig)
res3i = glmer(copedevratio ~ eggdev +momage+(1|momID), data=third, family=binomial)
summary(res3i)
exp(-0.21632)
exp(-0.21632-1.96*0.13003)
exp(-0.21632+1.96*0.13003)

#Estimate covariances among and between mothers
m3i <- MCMCglmm(cbind(copedevratio, eggdev) ~ trait-1,
                random=~us(trait):momID,
                rcov = ~ us(trait):units, 
                prior = prior,
                family =  c("gaussian", "poisson"), nitt = 100000, burnin = 20000,
                thin = 25, data = third)
summary(m3i)

####Plotting for res 3###################### 

#1Trade-off between clutch size and offspring size (control for mom age) (sig)
#Figure 3
#Run simple model to get predicted response and confidence intervals
res3plot = glm.nb(naup ~ areaum2, data=datum.no.missing.size)
summary(res3plot)
#Extract the predicted values and standard error
pred = predict(res3plot,type = "response", se.fit = TRUE)
pred2=interval(res3plot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3= pred2[,1]
uprres3 = pred2[,3]
lwrres3 = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p10 <- datum.no.missing.size %>% 
  mutate(pred = predsres3,
         upr= uprres3,
         lwr = lwrres3) %>% 
  ggplot(aes(x=areaum2, y=naup))  + geom_point(aes(y=naup), size=2)+
  xlab(expression(paste("Offspring body size (", uM^2,")"))) + ylab("Clutch size (# Indv.)")+
  geom_line(aes(y=predsres3), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3,ymax=uprres3), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=17), axis.title.x = element_text(size=17),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 7500, y = 33, label = expression(italic(beta)*"=-9.73,"~italic(p)*"=0.012"), 
           size = 5, color ="blue")

jpeg(file="Fig3.jpg", units="in", width=8, height=6, res=300)
p10
dev.off()


#2 Relationship between offspring size and survival (non-sig)
#Run simple model to get predicted response and confidence intervals
res3bplot = glm(survival ~ areaum2, data=datum.no.missing.size.surv, family= binomial)
summary(res3bplot)
#Extract the predicted values and standard error
pred = predict(res3bplot,type = "response", se.fit = TRUE)
pred2=interval(res3bplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3b= pred2[,1]
uprres3b = pred2[,3]
lwrres3b = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p11 <- datum.no.missing.size.surv %>% 
  mutate(pred = predsres3b,
         upr= uprres3b,
         lwr = lwrres3b) %>% 
  ggplot(aes(x=areaum2, y=survival))  + geom_point(aes(y=survival), size=2)+
  xlab("") + ylab("Clutch survival ratio")+
  geom_line(aes(y=predsres3b), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3b,ymax=uprres3b), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 19000, y = 0.13, label = expression(italic(beta)*"=11.94,"~italic(p)*"=0.575"), 
           size = 4, color ="blue")


#3 Relationship betweeen clutch size and survival (non-sig)
#Run simple model to get predicted response and confidence intervals
res3cplot = glm(survivalratio ~ naup, data=third, family = binomial)
summary(res3cplot)
#Extract the predicted values and standard error
pred = predict(res3cplot,type = "response", se.fit = TRUE)
pred2=interval(res3cplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3c= pred2[,1]
uprres3c = pred2[,3]
lwrres3c = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p12 <- third %>% 
  mutate(pred = predsres3c,
         upr= uprres3c,
         lwr = lwrres3c) %>% 
  ggplot(aes(x=naup, y=survivalratio))  + geom_point(aes(y=survivalratio), size=2)+
  xlab("") + ylab("")+
  geom_line(aes(y=predsres3c), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3c,ymax=uprres3c), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 32, y = 0.15, label = expression(italic(beta)*"=4.63,"~italic(p)*"=0.057"), 
           size = 4, color ="blue")


#4 Relationship between offspring size and development ratio (control for mom age) (non-sig)
#use scaled size data again
#Run simple model to get predicted response and confidence intervals
res3dplot = glm(development ~ areaum2, data=datum.no.missing.size.dev, family =binomial)
summary(res3dplot)
#Extract the predicted values and standard error
pred = predict(res3dplot,type = "response", se.fit = TRUE)
pred2=interval(res3dplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3d= pred2[,1]
uprres3d = pred2[,3]
lwrres3d = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p13 <- datum.no.missing.size.dev %>% 
  mutate(pred = predsres3d,
         upr= uprres3d,
         lwr = lwrres3d) %>% 
  ggplot(aes(x=areaum2, y=development))  + geom_point(aes(y=development), size=2)+
  xlab(expression(paste("Offspring body size (", uM^2,")"))) + ylab("Clutch development ratio")+
  geom_line(aes(y=predsres3d), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3d,ymax=uprres3d), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 19000, y = 0.13, label = expression(italic(beta)*"=2.89,"~italic(p)*"=0.876"), 
           size = 4, color ="blue")


#5 Relationship between clutch size and development ratio (control for mom age) (non-sig)
#Run simple model to get predicted response and confidence intervals
res3eplot = glm(copedevratio ~ naup, data=third, family = binomial)
summary(res3eplot)
#Extract the predicted values and standard error
pred = predict(res3eplot,type = "response", se.fit = TRUE)
pred2=interval(res3eplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3e= pred2[,1]
uprres3e = pred2[,3]
lwrres3e = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p14 <- third %>% 
  mutate(pred = predsres3e,
         upr= uprres3e,
         lwr = lwrres3e) %>% 
  ggplot(aes(x=naup, y=copedevratio))  + geom_point(aes(y=copedevratio), size=2)+
  xlab("Clutch size") + ylab("")+
  geom_line(aes(y=predsres3e), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3e,ymax=uprres3e), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12), axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 32, y = 0.15, label = expression(italic(beta)*"=3.45,"~italic(p)*"=0.135"), 
           size = 4, color ="blue")


#combine 11 through 14 together
figureS5 <- ggarrange(p11, p12, p13, p14, labels = c("A", "B", "C", "D"),  ncol=2, nrow=2)

jpeg(file="FigS5.jpg", units="in", width=9, height=6.5, res=300)
figureS5
dev.off()


#6relationship between egg dev and size (control for mom age) (sig)
p15 <- ggplot(datum, aes(x=eggdev, y=areaum2))  + geom_point(aes(y=areaum2), size=2)+
  ylim(4000,22000)+
  xlab("") + ylab(expression(paste("Offspring body size (", uM^2,")")))+
  geom_smooth(method='lm', color = "black", lwd = 2)+
  theme(axis.text.x = element_text(size=9))+
  theme(axis.text.y = element_text(size=9))+
  theme(axis.title.y = element_text(size=10))+
  theme(axis.title.x = element_text(size=11))+
  annotate("text", x = 13, y = 20000, label = expression(italic(beta)*"=316.59,"~italic(p)*"<0.001"), 
           size = 3, color ="blue")

#7relationship between egg dev and clutch size (control for mom age) (sig)
#Run simple model to get predicted response and confidence intervals
res3gplot = glm(naup ~ eggdev, data=third, family = poisson)
summary(res3gplot)
#Extract the predicted values and standard error
pred = predict(res3gplot,type = "response", se.fit = TRUE)
pred2=interval(res3gplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3g= pred2[,1]
uprres3g = pred2[,3]
lwrres3g = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p16 <- third %>% 
  mutate(pred = predsres3g,
         upr= uprres3g,
         lwr = lwrres3g) %>% 
  ggplot(aes(x=eggdev, y=naup))  + geom_point(aes(y=naup), size=2)+
  xlab("") + ylab("Clutch Size")+
  geom_line(aes(y=predsres3g), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3g,ymax=uprres3g), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_text(size=11),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 13, y = 32, label = expression(italic(beta)*"=-10.57,"~italic(p)*"=0.003"), 
           size = 3, color ="blue")


#8relationship between egg dev and survival ratio (sig)
#Run simple model to get predicted response and confidence intervals
res3hplot = glm(survivalratio ~ eggdev, data=third, family = binomial)
summary(res3hplot)
#Extract the predicted values and standard error
pred = predict(res3hplot,type = "response", se.fit = TRUE)
pred2=interval(res3hplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3h= pred2[,1]
uprres3h = pred2[,3]
lwrres3h = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p17 <- third %>% 
  mutate(pred = predsres3h,
         upr= uprres3h,
         lwr = lwrres3h) %>% 
  ggplot(aes(x=eggdev, y=survivalratio))  + geom_point(aes(y=survivalratio), size=2)+
  xlab("Egg gestation time (days)") + ylab("Clutch survival ratio")+
  geom_line(aes(y=predsres3h), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3h,ymax=uprres3h), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=11), axis.title.x = element_text(size=11),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 13, y = 0.93, label = expression(italic(beta)*"=-26.63,"~italic(p)*"=0.026"), 
           size = 3, color ="blue")


#9relationship between egg dev and development ratio (control for mom age) (nonsig)
#Run simple model to get predicted response and confidence intervals
res3iplot = glm(copedevratio ~ eggdev, data=third, family = binomial)
summary(res3iplot)
#Extract the predicted values and standard error
pred = predict(res3iplot,type = "response", se.fit = TRUE)
pred2=interval(res3iplot, pred, conf.level = 0.95)
#pred2 is an atomic vector that has the fitted values upper and lower ci
#extract those to objects
predsres3i= pred2[,1]
uprres3i = pred2[,3]
lwrres3i = pred2[,2]
#Get the predicted, upper and lower objects into the dataframe used for plotting
p18 <- third %>% 
  mutate(pred = predsres3i,
         upr= uprres3i,
         lwr = lwrres3i) %>% 
  ggplot(aes(x=eggdev, y=copedevratio))  + geom_point(aes(y=copedevratio), size=2)+
  xlab("Egg gestation time (days)") + ylab("Clutch development ratio")+
  geom_line(aes(y=predsres3i), lwd = 1.5)+     #This adds a line using the predicted values
  geom_ribbon(aes(ymin=lwrres3i,ymax=uprres3i), alpha =.15) + # This adds the confidence intervals
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=10), axis.title.x = element_text(size=11),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x = 13, y = 0.93, label = expression(italic(beta)*"=-19.46,"~italic(p)*"=0.096"), 
           size = 3, color ="blue")


#combine 15 through 18 together
figureS6 <- ggarrange(p15, p16, p17, p18, labels = c("A", "B", "C", "D"),  ncol=2, nrow=2)

jpeg(file= "FigS6.jpg", units = "in", width = 7, height =5, res=300)
figureS6
dev.off()