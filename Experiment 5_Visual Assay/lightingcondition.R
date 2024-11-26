#Visual Signal

#libraries
library(tidyverse)
library(here)
library(reshape2)
library(ggpattern)
library(lme4)
library(data.table)
library(ggeffects)

trialdata <- read_csv(here("lightingconditionsummary.csv"), col_names = TRUE)

#behaviorsummary:
behaviorsummary <- trialdata %>%
  group_by(lighting, grouped_age) %>%
  summarize(meanbehaviors = mean(total_aggressive_behaviors, na.rm=TRUE))

#overall behaviors:
ggplot(trialdata, aes(x=grouped_age, y=totalaggressivebehaviors_normalized, color=as.factor(lighting), group=as.factor(lighting)))+
  geom_point(size=3, alpha=0.6, position = position_dodge(width = 0.5))+
  geom_path(behaviorsummary, mapping=aes(x=grouped_age, y=meanbehaviors, color=as.factor(lighting)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.125, 0.875),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Lighting Condition"))+
  scale_color_manual(values=c("#EEAD0E", "#6C7B8B"))+
  ylab("Aggressive Behaviors/5 min")+
  xlab("Age (Days Post Fertilization)")


#run the model:

model <- glmer(formula = total_aggressive_behaviors ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#output
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
# ]
# Family: poisson  ( log )
# Formula: 
#   total_aggressive_behaviors ~ 1 + grouped_age + lighting + grouped_age:lighting +  
#   (1 | tank_num)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 821.3    830.2   -405.6    811.3       39 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -5.9659 -2.0460 -0.8092  0.6692 18.3489 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# tank_num (Intercept) 0.8843   0.9404  
# Number of obs: 44, groups:  tank_num, 12
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -5.36177    0.79755  -6.723 1.78e-11 ***
#   grouped_age                0.28245    0.02589  10.908  < 2e-16 ***
#   lightingLight              4.90153    0.91758   5.342 9.20e-08 ***
#   grouped_age:lightingLight -0.13178    0.02752  -4.788 1.69e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) grpd_g lghtnL
# grouped_age -0.864              
# lightngLght -0.869  0.752       
# grpd_g:lghL  0.813 -0.941 -0.794

#bootstrap

set.seed(4523957)
set.seed(7248) #flipping the model so we test dark against lighting
set.seed(3828332)

trialdata <- as.data.table(trialdata)

lightingbootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  
  tryCatch(
    {
      model <- glmer(formula = total_aggressive_behaviors ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), 
                     data = df_sample, 
                     family=poisson
      )
      summary <- summary(model)
      
      rbind(
        Intercept = summary$coefficients[1,1],
        Age_Coef = summary$coefficients[2,1],
        lighting_Coef = summary$coefficients[3,1],
        agelighting_Coef = summary$coefficients[4,1],
        pAge = summary$coefficients[2,4],
        plighting = summary$coefficients[3,4],
        pAgeLight = summary$coefficients[4,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, lightingbootstrap(trialdata, 12)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "lighting_Coef", "agelighting_Coef", "pAge", "plighting", "pAgeLight")

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            plighting = sum(plighting < 0.05, na.rm=T)/length(plighting),
            pAgeLight = sum(pAgeLight < 0.05, na.rm=T)/length(pAgeLight))

#confidence intervals: 
Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
Lighting <- quantile(bootstrapresults$lighting_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
AgeLight <- quantile(bootstrapresults$agelighting_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals <- cbind(Intercept, Lighting, Age, AgeLight)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="#EEAD0E", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=lighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.25)+
  theme_classic()+
  coord_Cartesian()

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="#EEAD0E", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=agelighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.01)+
  theme_classic()

#model visualization
predicteddata <- predict_response(model, terms = c("grouped_age ", "lighting"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(trialdata, mapping=aes(x=grouped_age, y=total_aggressive_behaviors, color=lighting, group=lighting),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#6C7B8B", "#EEAD0E"))+
  scale_fill_manual(values=c("#6C7B8B", "#EEAD0E"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")


fwrite(bootstrappower, file = "bootstrappower_lighting.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_lighting.csv", row.names=FALSE)

#by individual
ggplot(trialdata, aes(x=grouped_age, y=total_aggressive_behaviors, color=as.factor(lighting), group=as.factor(tank_num)))+
  geom_point(size=3, alpha=0.6)+
  geom_path(size=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.125, 0.875),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Lighting Condition"))+
  scale_color_manual(values=c("#6C7B8B", "#EEAD0E"))+
  ylab("Aggressive Behaviors/5 min")+
  xlab("Age (Days Post Fertilization)")

####################################

#shift in behaviors...let me see

trialdata_shift <- trialdata %>%
  select(grouped_age, lighting, chases, bites, lateral_displays)
trialdata_shift <- reshape2::melt(trialdata_shift, id=c("grouped_age", "lighting"))
names(trialdata_shift) <- c("grouped_age", "lighting", "Behavior", "Count")

#behaviorsummary:
behaviorsummary <- trialdata_shift %>%
  group_by(lighting, grouped_age, Behavior, BehaviorID) %>%
  summarize(meanbehaviorcount = mean(Count, na.rm=TRUE))

#cumsum to demonstrate individual behaviors
trialdata_shift <- trialdata_shift %>%
  arrange(lighting, Behavior)
trialdata_shift$BehaviorID <- cumsum(!duplicated(trialdata_shift[c(2,3)]))

#plot
ggplot(trialdata_shift, aes(x=grouped_age, y=Count, color=lighting, group=BehaviorID))+
  geom_bar_pattern(behaviorsummary, mapping=aes(x=grouped_age,  y=meanbehaviorcount, fill=lighting, group=BehaviorID, pattern = Behavior, pattern_angle = Behavior), 
                  pattern_spacing = 0.025,
                   stat="identity", position=position_dodge(width=3.5)) +
  geom_point(size=1.75, position = position_jitterdodge(dodge.width=3.5, jitter.width = 0.2), alpha=0.6)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position.inside = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(fill=guide_legend(title="Lighting Condition"))+
  guides(color=guide_legend(title="Lighting Condition"))+
  scale_fill_manual(values=c("#6C7B8B", "#EEAD0E"))+
  scale_color_manual(values=c("#21262a", "#503b05"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

##################################################################################

#we also want to run the model on just chases and bites to demonstrate the shift in behavioral type:

#chases
model <- glmer(formula = chases ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#output
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: chases ~ 1 + grouped_age + lighting + grouped_age:lighting +      (1 | ID)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 600.1    609.0   -295.1    590.1       39 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -5.0198 -1.2326 -0.3645  0.4858 15.0823 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 0.986    0.993   
# Number of obs: 44, groups:  ID, 12
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -6.13106    1.25905  -4.870 1.12e-06 ***
#   grouped_age                0.26701    0.04427   6.032 1.62e-09 ***
#   lightingLight              4.85776    1.35734   3.579 0.000345 ***
#   grouped_age:lightingLight -0.10192    0.04577  -2.227 0.025963 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) grpd_g lghtnL
# grouped_age -0.935              
# lightngLght -0.928  0.868       
# grpd_g:lghL  0.905 -0.967 -0.895

#model bootstrap :P

set.seed(43284)

trialdata <- as.data.table(trialdata)

chasebootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  
  tryCatch(
    {
      model <- glmer(formula = chases ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), 
                     data = df_sample, 
                     family=poisson
      )
      summary <- summary(model)
      
      rbind(
        Intercept = summary$coefficients[1,1],
        Age_Coef = summary$coefficients[2,1],
        lighting_Coef = summary$coefficients[3,1],
        agelighting_Coef = summary$coefficients[4,1],
        pAge = summary$coefficients[2,4],
        plighting = summary$coefficients[3,4],
        pAgeLight = summary$coefficients[4,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults_chases <- as.data.frame(replicate(1000, chasebootstrap(trialdata, 12)))
bootstrapresults_chases <- t(bootstrapresults_chases)
bootstrapresults_chases <- as.data.frame(bootstrapresults_chases)
names(bootstrapresults_chases) <- c("Intercept", "Age_Coef", "lighting_Coef", "agelighting_Coef", "pAge", "plighting", "pAgeLight")

bootstrappower_chases <- bootstrapresults_chases %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            plighting = sum(plighting < 0.05, na.rm=T)/length(plighting),
            pAgeLight = sum(pAgeLight < 0.05, na.rm=T)/length(pAgeLight))

#confidence intervals: 
Intercept <- quantile(bootstrapresults_chases$Intercept, prob=c(.025,.5,.975), na.rm=T)
Lighting <- quantile(bootstrapresults_chases$lighting_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults_chases$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
AgeLight <- quantile(bootstrapresults_chases$agelighting_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals_chases <- cbind(Intercept, Lighting, Age, AgeLight)
ConfidenceIntervals_chases <- data.frame(ConfidenceIntervals_chases)

ggplot()+
  geom_histogram(bootstrapresults_chases, mapping=aes(x=Intercept), fill="#EEAD0E", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults_chases, mapping=aes(x=lighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults_chases, mapping=aes(x=Age_Coef), fill="#EEAD0E", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults_chases, mapping=aes(x=agelighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.01)+
  theme_classic()

#model visualization
predicteddata <- predict_response(model, terms = c("grouped_age ", "lighting"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(trialdata, mapping=aes(x=grouped_age, y=chases, color=lighting, group=lighting),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#6C7B8B", "#EEAD0E"))+
  scale_fill_manual(values=c("#6C7B8B", "#EEAD0E"))+
  ylab("Chases/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

fwrite(bootstrappower_chases, file = "bootstrappower_lighting_chasesonly.csv", row.names=FALSE)
fwrite(ConfidenceIntervals_chases, file = "confidenceintervals_lighting_chasesonly.csv", row.names=FALSE)

###################################################################################

#now do it again for BITES

model <- glmer(formula = bites ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#output
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: bites ~ 1 + grouped_age + lighting + grouped_age:lighting + (1 |      ID)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 230.3    239.2   -110.1    220.3       39 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.0721 -0.8876 -0.5455  0.4171  5.1747 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 0.4108   0.6409  
# Number of obs: 44, groups:  ID, 12
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -5.87338    0.91055  -6.450 1.12e-10 ***
#   grouped_age                0.28865    0.03227   8.945  < 2e-16 ***
#   lightingLight              1.38137    1.23466   1.119    0.263    
# grouped_age:lightingLight -0.04951    0.04393  -1.127    0.260    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) grpd_g lghtnL
# grouped_age -0.947              
# lightngLght -0.741  0.706       
# grpd_g:lghL  0.698 -0.739 -0.944

set.seed(454269)

trialdata <- as.data.table(trialdata)

bitesbootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  
  tryCatch(
    {
      model <- glmer(formula = bites ~ 1 + grouped_age + lighting + grouped_age:lighting + (1|ID), 
                     data = df_sample, 
                     family=poisson
      )
      summary <- summary(model)
      
      rbind(
        Intercept = summary$coefficients[1,1],
        Age_Coef = summary$coefficients[2,1],
        lighting_Coef = summary$coefficients[3,1],
        agelighting_Coef = summary$coefficients[4,1],
        pAge = summary$coefficients[2,4],
        plighting = summary$coefficients[3,4],
        pAgeLight = summary$coefficients[4,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults_bites <- as.data.frame(replicate(1000, chasebootstrap(trialdata, 12)))
bootstrapresults_bites <- t(bootstrapresults_bites)
bootstrapresults_bites <- as.data.frame(bootstrapresults_bites)
names(bootstrapresults_bites) <- c("Intercept", "Age_Coef", "lighting_Coef", "agelighting_Coef", "pAge", "plighting", "pAgeLight")

bootstrappower_bites <- bootstrapresults_bites %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            plighting = sum(plighting < 0.05, na.rm=T)/length(plighting),
            pAgeLight = sum(pAgeLight < 0.05, na.rm=T)/length(pAgeLight))

#confidence intervals: 
Intercept <- quantile(bootstrapresults_bites$Intercept, prob=c(.025,.5,.975), na.rm=T)
Lighting <- quantile(bootstrapresults_bites$lighting_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults_bites$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
AgeLight <- quantile(bootstrapresults_bites$agelighting_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals_bites <- cbind(Intercept, Lighting, Age, AgeLight)
ConfidenceIntervals_bites <- data.frame(ConfidenceIntervals_bites)

ggplot()+
  geom_histogram(bootstrapresults_bites, mapping=aes(x=Intercept), fill="#EEAD0E", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults_bites, mapping=aes(x=lighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults_bites, mapping=aes(x=Age_Coef), fill="#EEAD0E", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults_bites, mapping=aes(x=agelighting_Coef), fill="#6C7B8B", alpha=0.5, binwidth=0.01)+
  theme_classic()

#model visualization
predicteddata <- predict_response(model, terms = c("grouped_age ", "lighting"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(trialdata, mapping=aes(x=grouped_age, y=bites, color=lighting, group=lighting),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#6C7B8B", "#EEAD0E"))+
  scale_fill_manual(values=c("#6C7B8B", "#EEAD0E"))+
  ylab("Bites/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

fwrite(bootstrappower_bites, file = "bootstrappower_lighting_bitessonly.csv", row.names=FALSE)
fwrite(ConfidenceIntervals_bites, file = "confidenceintervals_lighting_bitesonly.csv", row.names=FALSE)


