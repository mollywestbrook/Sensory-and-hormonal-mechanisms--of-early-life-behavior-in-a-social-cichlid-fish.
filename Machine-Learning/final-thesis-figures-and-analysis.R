#Behavioral Feature Correlations

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)
library(ggeffects)
library(lme4)
library(r2glmm)

#This section of the script brought in each dataset, and assigned them an associated ID and metadata

featuredf <- read.csv(here("18-05-29-2024_piz1_T2_27dpf_KT_features.csv"))

formatteddf_18 <- featuredf %>%
  select(frame, track, velocity_mmpersec, groupdistance, avgnn, pol, Behavior) %>% #grab only the features we need from the features dataset
  mutate(treatment = "KT", #bind in treatment column
         age = 27, #bind in an age column
         ID = 18) #bind in the ID of this video

#save the subsequent formatted dataset

fwrite(formatteddf_18, file = "18-05-29-2024_piz1_T2_27dpf_KT_final.csv", row.names=FALSE)

beep(sound=11)

#stitch together all of the datasets into a large dataset for aggreggation and final analysis

bigdf <- rbind(formatteddf_1, formatteddf_2, formatteddf_3, formatteddf_4, formatteddf_5, formatteddf_6, formatteddf_7,
               formatteddf_8, formatteddf_9, formatteddf_10, formatteddf_11, formatteddf_12, formatteddf_13, formatteddf_14,
               formatteddf_15, formatteddf_16, formatteddf_17, formatteddf_18)
fwrite(bigdf, file = "bigdf_forcorrelations.csv", row.names=FALSE) #save it

##################################################

bigdf <- read.csv(here("bigdf_forcorrelations.csv")) #load in the big dataset

aggressivebehaviors <- read.csv(here("totalcountofbehaviors.csv"))


#two fixes: we'll not include df 14, not worth it
#and we need to modify df12age
bigdf <- bigdf %>%
  filter(!ID == 14) %>%
  mutate(age = case_when(ID == 12 ~ 18, .default = age))

#This chunk makes a summary figure to show all of the datasets and metadata of each of them. 
sleapinfo <- bigdf %>%
  group_by(ID) %>%
  summarize(Treatment = unique(treatment),
            Age = unique(age))
ggplot(sleapinfo, aes(x=as.factor(Age), y=Treatment, color=Treatment, size=as.factor(Age)))+
  geom_point(alpha=0.7, position=position_jitterdodge(jitter.height=0.3))+
  theme_classic()+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  theme(text = element_text(size=18),
        legend.position="none")+
  ylab("Treatment")+
  xlab("Age (DPF)")

#Figure 2a: Velocity correlated with aggressive behavior count

#summarize across the 5 groups:
bigdf_averagegroupvelocity <- bigdf %>%
  group_by(ID, frame, age, treatment) %>%
  summarize(averagevelocity = mean(velocity_mmpersec, na.rm=T),
            sdvelocity = sd(velocity_mmpersec, na.rm=T))
#supplemental figure 2
ggplot(bigdf_averagegroupvelocity, aes(x=as.factor(ID), y=averagevelocity, fill=treatment, color=age))+
  geom_violin()+
  theme_classic()+
  theme(text = element_text(size = 15))+ 
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  coord_cartesian(ylim=c(0, 35))+
  xlab("TreatmentID")+
  ylab("Average Velocity (mm/sec)")

#aggregate velocity + correlated with number of aggressive behaviors
velocitysummary <- bigdf_averagegroupvelocity %>%
  group_by(ID, age, treatment) %>%
  summarize(medianvelocity = median(averagevelocity, na.rm=T),
            sdvelocity = sd(averagevelocity, na.rm=T))
#merge in total counts of behavior
velocitysummary <- merge(velocitysummary, aggressivebehaviors, by=c("ID"), all.x=T)

#model:
velocitymodel <- lmer(aggressive_behaviors ~ 1 + medianvelocity + (1|individual), data=velocitysummary)
summary <- summary(velocitymodel)

# Linear mixed model fit by REML ['lmerMod']
# Formula: aggressive_behaviors ~ 1 + medianvelocity + (1 | individual)
# Data: velocitysummary
# 
# REML criterion at convergence: 163.2
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.37602 -0.38816  0.09994  0.42429  1.19631 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# individual (Intercept) 1287     35.88   
# Residual                808     28.43   
# Number of obs: 17, groups:  individual, 13
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)     -56.597     23.057  -2.455
# medianvelocity   22.978      3.224   7.127
# 
# Correlation of Fixed Effects:
#   (Intr)
# medianvlcty -0.846

#model residuals
velocitymodel$residuals
hist(velocitymodel$residuals)
# 
#model heteroskedacity
res <- velocitymodel$residuals
plot(fitted(velocitymodel), res)
abline(0,0)

#calculate R2
r2beta(velocitymodel, data=velocitysummary)

# Effect   Rsq upper.CL lower.CL
# 1          Model 0.712    0.874    0.487
# 2 medianvelocity 0.712    0.874    0.487

#spearman correlation
cor.test(velocitysummary$medianvelocity, velocitysummary$aggressive_behaviors, method = 'spearman')

# Spearman's rank correlation rho
# 
# data:  velocitysummary$medianvelocity and velocitysummary$aggressive_behaviors
# S = 153.19, p-value = 7.46e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8122706 

#predict model:
predicteddata <- predict_response(velocitymodel, terms = c("medianvelocity"))
ggplot(data=predicteddata, aes(x=x, y=predicted))+
  geom_ribbon(mapping=aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(data=velocitysummary, mapping=aes(x=medianvelocity, y=aggressive_behaviors, color=treatment, size=age))+
  theme_classic()+
  scale_size_continuous(range = c(4,8))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  theme(text = element_text(size = 15),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Median Group Velocity (mm/sec)")

#mean value of medians: 6.03 mm/sec

####################################################################################

#Figure 2b: Average Nearest Neighbor correlated with Total Counts of Aggressive Behavior

#average nearest neighbor:
bigdf_averagenearestneighbor <- bigdf %>%
  group_by(ID, frame, age, treatment) %>%
  summarize(averagenearestneighbor = mean(avgnn, na.rm=T))
#another summary plot -- this is not used in the document
ggplot(bigdf_averagenearestneighbor, aes(x=as.factor(ID), y=averagenearestneighbor, fill=treatment, color=age))+
  geom_violin()+
  theme_classic()+
  theme(text = element_text(size = 15))+ 
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  xlab("TreatmentID")+
  ylab("Average Nearest Neighbor (mm)")
#aggregate and bind in counts of behavior
avgnnsummary <- bigdf_averagenearestneighbor %>%
  group_by(ID, age, treatment) %>%
  summarize(avgnn = mean(averagenearestneighbor, na.rm=T),
            sdavgnn = sd(averagenearestneighbor, na.rm=T))
avgnnsummary <- merge(avgnnsummary, aggressivebehaviors, by=c("ID"), all.x=T)

nnmodel <- lmer(aggressive_behaviors ~ 1 + avgnn + (1|individual), data=avgnnsummary)
summary <- summary(nnmodel)

# Linear mixed model fit by REML ['lmerMod']
# Formula: aggressive_behaviors ~ 1 + avgnn + (1 | individual)
# Data: avgnnsummary
# 
# REML criterion at convergence: 183
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.4135 -0.6878 -0.2264  1.0700  1.8551 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# individual (Intercept)    0      0.00   
# Residual               6457     80.35   
# Number of obs: 17, groups:  individual, 13
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  265.224    142.098   1.866
# avgnn         -5.355      4.032  -1.328
# 
# Correlation of Fixed Effects:
#   (Intr)
# avgnn -0.991
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')

#calculate R2
r2beta(nnmodel, data=avgnnsummary)

# Effect   Rsq upper.CL lower.CL
# 1  Model 0.104    0.478        0
# 2  avgnn 0.104    0.478        0

cor.test(avgnnsummary$avgnn, avgnnsummary$aggressive_behaviors, method = 'spearman')
# 
# Spearman's rank correlation rho
# 
# data:  avgnnsummary$avgnn and avgnnsummary$aggressive_behaviors
# S = 998.22, p-value = 0.3889
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2233131 

predicteddata <- predict_response(nnmodel, terms = c("avgnn"))
ggplot(data=predicteddata, aes(x=x, y=predicted))+
  geom_ribbon(mapping=aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(data=avgnnsummary, mapping=aes(x=avgnn, y=aggressive_behaviors, color=treatment, size=age))+
  theme_classic()+
  scale_size_continuous(range = c(4,8))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  theme(text = element_text(size = 15),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Average Nearest Neighbor (mm)")

#mean = 34.9053

#################################################################################

#Figure 2c: Group Distance correlated with Total Counts of Aggressive Behavior

#group distance
bigdf_groupdistance <- bigdf %>%
  group_by(ID, frame, age, treatment) %>%
  summarize(avggroupdistance = mean(groupdistance, na.rm=T))
#another summary plot -- this is not used in the document
ggplot(bigdf_groupdistance, aes(x=as.factor(ID), y=avggroupdistance, fill=treatment, color=age))+
  geom_violin()+
  theme_classic()+
  theme(text = element_text(size = 15))+ 
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  xlab("TreatmentID")+
  ylab("Average Group Distance (mm)")

groupdistancesummary <- bigdf_groupdistance %>%
  group_by(ID, age, treatment) %>%
  summarize(avggroupdistance = mean(avggroupdistance, na.rm=T),
            sdgroupdistance = sd(avggroupdistance, na.rm=T))
groupdistancesummary <- merge(groupdistancesummary, aggressivebehaviors, by=c("ID"), all.x=T)

groupdistncemodel <- lmer(aggressive_behaviors ~ 1 + avggroupdistance + (1|individual), data=groupdistancesummary)
summary <- summary(groupdistncemodel)

# Linear mixed model fit by REML ['lmerMod']
# Formula: aggressive_behaviors ~ 1 + avggroupdistance + (1 | individual)
# Data: groupdistancesummary
# 
# REML criterion at convergence: 183.5
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.4728 -0.7056 -0.2598  1.1212  1.7994 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# individual (Intercept)    0      0.0    
# Residual               6130     78.3    
# Number of obs: 17, groups:  individual, 13
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)       285.976    128.849   2.219
# avggroupdistance   -3.283      2.014  -1.630
# 
# Correlation of Fixed Effects:
#   (Intr)
# avggrpdstnc -0.989
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')

#calculate R2
r2beta(groupdistncemodel, data=groupdistancesummary)

# Effect   Rsq upper.CL lower.CL
# 1            Model 0.149    0.527    0.001
# 2 avggroupdistance 0.149    0.527    0.001

cor.test(groupdistancesummary$avggroupdistance, groupdistancesummary$aggressive_behaviors, method = 'spearman')

# Spearman's rank correlation rho
# 
# data:  groupdistancesummary$avggroupdistance and groupdistancesummary$aggressive_behaviors
# S = 1068.3, p-value = 0.2272
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.3092027 

predicteddata <- predict_response(groupdistncemodel, terms = c("avggroupdistance"))
ggplot(data=predicteddata, aes(x=x, y=predicted))+
  geom_ribbon(mapping=aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(data=groupdistancesummary, mapping=aes(x=avggroupdistance, y=aggressive_behaviors, color=treatment, size=age))+
  theme_classic()+
  scale_size_continuous(range = c(4,8))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  theme(text = element_text(size = 15),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Average Group Distance (mm)")

#mean: 63.26412

#Supplemental Figure 3: AVGNN and Group Distance are Correlated

distancecheck <- data.frame(avgnnsummary$avgnn, groupdistancesummary$avggroupdistance)
names(distancecheck) <- c("avgnn", "groupdistance")

ggplot(distancecheck, aes(x=groupdistance, y=avgnn))+
  geom_point()+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Average Nearest Neighbor (mm)")+
  xlab("Group Distance (mm)")

####################################################################################

#Figure 2d: Polarity correlated with Total Counts of Aggressive Behavior

bigdf_pol <- bigdf %>%
  group_by(ID, frame, age, treatment) %>%
  summarize(avgpol = mean(pol, na.rm=T))
#supplementary figure not included in the dataset
ggplot(bigdf_pol, aes(x=as.factor(ID), y=avgpol, fill=treatment, color=age))+
  geom_violin()+
  theme_classic()+
  theme(text = element_text(size = 15))+ 
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  xlab("TreatmentID")+
  ylab("Average Polarity")

polaritysummary <- bigdf_pol %>%
  group_by(ID, age, treatment) %>%
  summarize(avgpol = mean(avgpol, na.rm=T),
            sdpol = sd(avgpol, na.rm=T))
polaritysummary <- merge(polaritysummary, aggressivebehaviors, by=c("ID"), all.x=T)

polaritymodel <- lmer(aggressive_behaviors ~ 1 + avgpol + (1|individual), data=polaritysummary)
summary <- summary(polaritymodel)

# Linear mixed model fit by REML ['lmerMod']
# Formula: aggressive_behaviors ~ 1 + avgpol + (1 | individual)
# Data: polaritysummary
# 
# REML criterion at convergence: 165.5
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.42729 -0.82854 -0.05491  0.84889  1.58169 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# individual (Intercept)    0      0.00   
# Residual               4611     67.91   
# Number of obs: 17, groups:  individual, 13
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  -2008.8      717.2  -2.801
# avgpol        5008.1     1720.6   2.911
# 
# Correlation of Fixed Effects:
#   (Intr)
# avgpol -1.000
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')

#calculate R2
r2beta(polaritymodel, data=polaritysummary)

# Effect   Rsq upper.CL lower.CL
# 1  Model 0.358    0.688    0.055
# 2 avgpol 0.358    0.688    0.055

cor.test(polaritysummary$avgpol, polaritysummary$aggressive_behaviors, method = 'spearman')

# Spearman's rank correlation rho
# 
# data:  polaritysummary$avgpol and polaritysummary$aggressive_behaviors
# S = 404.5, p-value = 0.03899
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5042949 

predicteddata <- predict_response(polaritymodel, terms = c("avgpol"))
ggplot(data=predicteddata, aes(x=x, y=predicted))+
  geom_ribbon(mapping=aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(data=polaritysummary, mapping=aes(x=avgpol, y=aggressive_behaviors, color=treatment, size=age))+
  theme_classic()+
  scale_size_continuous(range = c(4,8))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5", "black"))+
  theme(text = element_text(size = 15),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Average Polarity")

#mean polarity:0.4167377
