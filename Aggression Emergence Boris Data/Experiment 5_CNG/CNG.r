#CNG

library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)

trialdata <- read_csv(here("cngsummary.csv"), col_names = TRUE)

trialdata <- trialdata %>%
  arrange(Experiment, Tank)
trialdata$ID <- cumsum(!duplicated(trialdata[c(1,2)]))

behaviorsummary <- trialdata %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))
ggplot(trialdata, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#009E73", "#D55E00", "#56B4E9"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#stats

model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#output

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [
#   glmerMod]
# Family: poisson  ( log )
# Formula: Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1 | ID)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 1477.8   1494.9   -731.9   1463.8       77 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -10.1176  -1.9989  -0.9732   0.6255  12.9938 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 1.023    1.011   
# Number of obs: 84, groups:  ID, 14
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  -4.09405    0.59596  -6.870 6.44e-12 ***
#   Age                           0.24080    0.01404  17.157  < 2e-16 ***
#   TreatmentCNG D52I2/+          1.72493    1.06851   1.614   0.1065    
# TreatmentCNG D52I2/D52I2      3.07284    0.73903   4.158 3.21e-05 ***
#   Age:TreatmentCNG D52I2/+     -0.04451    0.02352  -1.892   0.0585 .  
# Age:TreatmentCNG D52I2/D52I2 -0.10467    0.01606  -6.518 7.11e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) Age    TCNGD52I2/+ TCNGD52I2/D A:TCNGD52I2/+
#   Age           -0.623                                             
# TCNGD52I2/+   -0.558  0.348                                      
# TCNGD52I2/D   -0.806  0.503  0.450                               
# A:TCNGD52I2/+  0.372 -0.597 -0.597      -0.300                   
# A:TCNGD52I2/D  0.545 -0.873 -0.304      -0.568       0.521     

#bootstrap
set.seed(20387493)

#try grouping 30
set.seed(6222930)

trialdata <- as.data.table(trialdata)

cngbootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  
  tryCatch(
    {
      model <- glmer(
        formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), 
        data = df_sample, 
        family = poisson
      )
      summary <- summary(model)
      
      rbind(
        Intercept = summary$coefficients[1,1],
        Age_Coef = summary$coefficients[2,1],
        Het_Coef = summary$coefficients[3,1],
        Homo_Coef = summary$coefficients[4,1],
        HetAge_Coef = summary$coefficients[5,1],
        HomoAge_Coef = summary$coefficients[6,1],
        pAge = summary$coefficients[2,4],
        pHet = summary$coefficients[3,4],
        pHomo = summary$coefficients[4,4],
        pHetAge = summary$coefficients[5,4],
        pHomoAge = summary$coefficients[6,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, cngbootstrap(trialdata, 14)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "Het_Coef", "Homo_Coef", "HetAge_Coef", "HomoAge_Coef", 
                             "pAge", "pHet", "pHomo", "pHetAge", "pHomoAge")

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            pHet = sum(pHet < 0.05, na.rm=T)/length(pHet),
            pHomo = sum(pHomo < 0.05, na.rm=T)/length(pHomo),
            pHetAge = sum(pHetAge < 0.05, na.rm=T)/length(pHetAge),
            pHomoAge = sum(pHomoAge < 0.05, na.rm=T)/length(pHomoAge))

#confidence intervals: 
Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
Het <- quantile(bootstrapresults$Het_Coef, prob=c(.025,.5,.975), na.rm=T)
Homo <- quantile(bootstrapresults$Homo_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
HetAge <- quantile(bootstrapresults$HetAge_Coef, prob=c(.025,.5,.975), na.rm=T)
HomoAge <- quantile(bootstrapresults$HomoAge_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals <- cbind(Intercept, Het, Homo, Age, HetAge, HomoAge)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

#bootstrap histograms
ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="#009E73", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=Het_Coef), fill="#56B4E9", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=Homo_Coef), fill="#D55E00", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="#009E73", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=HetAge_Coef), fill="#56B4E9", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=HomoAge_Coef), fill="#D55E00", alpha=0.5, binwidth=0.01)+
  theme_classic()

#save all the power, CI etc
fwrite(bootstrappower, file = "bootstrappower_CNG.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_CNG.csv", row.names=FALSE)

#visualization of model:
#step one: predict effects:
predicteddata <- predict_response(model, terms = c("Age", "Treatment"))

ggplot(data=predicteddata, aes(x=x, y=predicted, group=group, color=group))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(trialdata, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#009E73", "#56B4E9", "#D55E00"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#################################################################################

#I will do the same thing but without the het group,as we only have 2 het individuals
#and it just mostly makes things more complicated

trialdata <- read_csv(here("cngsummary_nohet.csv"), col_names = TRUE)

trialdata <- trialdata %>%
  arrange(Experiment, Tank)
trialdata$ID <- cumsum(!duplicated(trialdata[c(1,2)]))

model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#model output:

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [
#   glmerMod]
# Family: poisson  ( log )
# Formula: Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1 | ID)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 1233.7   1245.0   -611.8   1223.7       67 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -10.0790  -1.8213  -0.8364   0.7415  11.8761 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 1.09     1.044   
# Number of obs: 72, groups:  ID, 12
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  -5.08128    0.66341  -7.659 1.87e-14 ***
#   Age                           0.28260    0.01751  16.138  < 2e-16 ***
#   TreatmentCNG D52I2/D52I2      3.01016    0.81382   3.699 0.000217 ***
#   Age:TreatmentCNG D52I2/D52I2 -0.10248    0.02002  -5.119 3.07e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) Age    TCNGD5
# Age         -0.689              
# TCNGD52I2/D -0.815  0.562       
# A:TCNGD52I2  0.603 -0.874 -0.637

set.seed(71723749)

trialdata <- as.data.table(trialdata)

cngbootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  
  tryCatch(
    {
      model <- glmer(
        formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), 
        data = df_sample, 
        family = poisson
      )
      summary <- summary(model)
      
      rbind(
        Intercept = summary$coefficients[1,1],
        Age_Coef = summary$coefficients[2,1],
        Homo_Coef = summary$coefficients[3,1],
        HomoAge_Coef = summary$coefficients[4,1],
        pAge = summary$coefficients[2,4],
        pHomo = summary$coefficients[3,4],
        pHomoAge = summary$coefficients[4,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, cngbootstrap(trialdata, 12)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "Homo_Coef", "HomoAge_Coef", 
                             "pAge", "pHomo","pHomoAge")

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            pHomo = sum(pHomo < 0.05, na.rm=T)/length(pHomo),
            pHomoAge = sum(pHomoAge < 0.05, na.rm=T)/length(pHomoAge))

#confidence intervals: 
Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
Homo <- quantile(bootstrapresults$Homo_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
HomoAge <- quantile(bootstrapresults$HomoAge_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals <- cbind(Intercept, Homo, Age, HomoAge)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

#bootstrap histograms
ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="#009E73", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=Homo_Coef), fill="#D55E00", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="#009E73", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=HomoAge_Coef), fill="#D55E00", alpha=0.5, binwidth=0.01)+
  theme_classic()

#save all the power, CI etc
fwrite(bootstrappower, file = "bootstrappower_CNG_nohet.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_CNG_nohet.csv", row.names=FALSE)

#visualization of model:
#step one: predict effects:
predicteddata <- predict_response(model, terms = c("Age", "Treatment"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(trialdata, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment", names=c("Wild Type", "Homozygous")), fill="none")+
  scale_color_manual(labels=c('Wild Type', 'Homozygous'), values=c("#009E73", "#D55E00"))+
  scale_fill_manual(labels=c('Wild Type', 'Homozygous'), values=c("#009E73", "#D55E00"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")


