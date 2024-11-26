#Thyroid Hormone Stuff

#libraries
library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)

thdf <- read_csv(here("Thyroidhormonesummary.csv"), col_names = TRUE)

behaviorsummary <- thdf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(thdf, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(Treatment)))+
  geom_point(size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#model

thdf <- thdf %>%
  arrange(Brood, Tank)
thdf$ID <- cumsum(!duplicated(thdf[c(1,3)])) #16 total individuals for bootstrap

model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), data = thdf, family=poisson)
summary <- summary(model)

#output
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1 | ID)
# Data: thdf
# 
# AIC      BIC   logLik deviance df.resid 
# 4683.1   4701.4  -2334.6   4669.1       93 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -11.585  -3.881  -2.456   1.015  26.138 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 0.305    0.5522  
# Number of obs: 100, groups:  ID, 16
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -1.097428   0.319776  -3.432 0.000599 ***
#   Age                     0.214170   0.005975  35.846  < 2e-16 ***
#   TreatmentThiourea       0.031833   0.421974   0.075 0.939866    
# TreatmentThyroxine     -1.099699   0.438207  -2.510 0.012089 *  
#   Age:TreatmentThiourea  -0.026662   0.008375  -3.184 0.001454 ** 
#   Age:TreatmentThyroxine  0.040053   0.009996   4.007 6.15e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) Age    TrtmntThr TrtmntThy Ag:TrtmntThr
# Age          -0.499                                        
# TreatmntThr  -0.758  0.378                                 
# TrtmntThyrx  -0.730  0.364  0.553                          
# Ag:TrtmntThr  0.356 -0.713 -0.527    -0.260                
# Ag:TrtmntThy  0.298 -0.598 -0.226    -0.573     0.426      

#bootstrap

set.seed(9267898)

thdf <- as.data.table(thdf)

thbootstrap <- function(df, n) {
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
        TU_Coef = summary$coefficients[3,1],
        TH_Coef = summary$coefficients[4,1],
        TUAge_Coef = summary$coefficients[5,1],
        THAge_Coef = summary$coefficients[6,1],
        pAge = summary$coefficients[2,4],
        pTU = summary$coefficients[3,4],
        pTH = summary$coefficients[4,4],
        pTUAge = summary$coefficients[5,4],
        pTHAge = summary$coefficients[6,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, thbootstrap(thdf, 16)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "TU_Coef", "TH_Coef", "TUAge_Coef", "THAge_Coef", "pAge", "pTU", "pTH", "pTUAge", "pTHAge")

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05, na.rm=T)/length(pAge),
            pTU = sum(pTU < 0.05, na.rm=T)/length(pTU),
            pTH = sum(pTH < 0.05, na.rm=T)/length(pTH),
            pTUAge = sum(pTUAge < 0.05, na.rm=T)/length(pTUAge),
            pTHAge = sum(pTHAge < 0.05, na.rm=T)/length(pTHAge))

#confidence intervals: 
Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
TU <- quantile(bootstrapresults$TU_Coef, prob=c(.025,.5,.975), na.rm=T)
TH <- quantile(bootstrapresults$TH_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
TUAge <- quantile(bootstrapresults$TUAge_Coef, prob=c(.025,.5,.975), na.rm=T)
THAge <- quantile(bootstrapresults$THAge_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals <- cbind(Intercept, TU, TH, Age, TUAge, THAge)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="black", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=TU_Coef), fill="#DC267F", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=TH_Coef), fill="#FE6100", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="black", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=TUAge_Coef), fill="#DC267F", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=THAge_Coef), fill="#FE6100", alpha=0.5, binwidth=0.01)+
  theme_classic()

#visualize model:
predicteddata <- predict_response(model, terms = c("Age", "Treatment"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(thdf, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  scale_fill_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")


#control only
predicteddata_naoh <- predicteddata %>%
  filter(group == "NaOH")

ggplot()+
  geom_ribbon(data=predicteddata_naoh, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata_naoh, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(thdf, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  scale_fill_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

fwrite(bootstrappower, file = "bootstrappower_th.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_th.csv", row.names=FALSE)

#by individual
ggplot(thdf, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(ID)))+
  geom_point(size=3, alpha=0.6)+
  geom_path(linewidth=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")


####################################################################################

#old model:
model <- lm(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment, data = thdf)
summary <- summary(model)
# 
# Call:
#   lm(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment, 
#      data = thdf)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -143.001  -39.357   -2.458   22.206  220.448 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -205.913     54.706  -3.764 0.000291 ***
#   Age                      12.264      2.411   5.086 1.86e-06 ***
#   TreatmentThiourea        88.820     70.653   1.257 0.211818    
# TreatmentThyroxine       20.629     76.193   0.271 0.787176    
# Age:TreatmentThiourea    -5.062      3.115  -1.625 0.107540    
# Age:TreatmentThyroxine   -1.196      3.475  -0.344 0.731541    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 63.04 on 94 degrees of freedom
# (12 observations deleted due to missingness)
# Multiple R-squared:  0.3951,	Adjusted R-squared:  0.3629 
# F-statistic: 12.28 on 5 and 94 DF,  p-value: 3.605e-09


#Weird Swimming Behavior
behaviorsummary <- thdf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(WeirdSwims, na.rm=TRUE))

ggplot(thdf, aes(x=Age, y=WeirdSwims, color=as.factor(Treatment), group=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Weird Swims")+
  xlab("Age (Days Post Fertilization)")

#####################################################################################

#physiology:

#Weight
thphysdata <- read_csv(here("thyroidhormonephysdata.csv"), col_names = TRUE)

ggplot(thphysdata, aes(x=Treatment, y=Weight, group=Treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=2, width=0.1)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  coord_cartesian(ylim=c(0, 20))+
  ylab("Weight (mg)")+
  xlab("Treatment")

ggplot(thphysdata, aes(x=Treatment, y=StandardLength, group=Treatment))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=2, width=0.1)+
  theme_classic()+
  theme(text=element_text(size=15))+
  coord_cartesian(ylim=c(0, 10))+
  ylab("Standard Length (mm)")+
  xlab("Treatment")
  
thphysdatasummary <- thphysdata %>%
  group_by(Treatment) %>%
  summarize(averageweight = mean(Weight, na.rm=T),
            weightsd = sd(Weight, na.rm=T),
            averagelength = mean(StandardLength, na.rm=T),
            lengthsd = sd(StandardLength, na.rm=T))

fwrite(thphysdatasummary, file = "thphysdatasummary.csv", row.names=FALSE)