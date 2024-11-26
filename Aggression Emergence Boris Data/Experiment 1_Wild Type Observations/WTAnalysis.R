#WT Aggression Emergence

library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)
library(paletteer)

trialdata <- read_csv(here("WTaggressionsummary.csv"), col_names = TRUE)

#figure

behaviorsummary <- trialdata %>%
  group_by(Age_DPF) %>%
  summarize(meanbehaviors = mean(TotalBehaviors, na.rm=TRUE))
ggplot(trialdata, aes(x=Age_DPF, y=TotalBehaviors))+
  geom_point(size=2, alpha=0.8)+
  geom_path(behaviorsummary, mapping=aes(x=Age_DPF, y=meanbehaviors), linewidth=1, stat="identity", na.rm=TRUE)+
  geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age_DPF"], linewidth=1)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#figure with individual IDs connected:

ggplot(trialdata, aes(x=Age_DPF, y=TotalBehaviors))+
  geom_point(size=2, alpha=0.8)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

ggplot(trialdata, aes(x=Age_DPF, y=TotalBehaviors, group=ID, color=as.factor(ID)))+
  geom_point(size=2, alpha=0.8)+
  geom_path(linewidth=1, alpha=0.5)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = "none")+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")+
  scale_color_paletteer_d("ggthemes::gdoc") 

#stats

#new model:

model <- glmer(formula = TotalBehaviors ~ 1 + Age_DPF + (1|ID), data = trialdata, family=poisson)
summary <- summary(model)

#results:
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [
#   glmerMod]
# Family: poisson  ( log )
# Formula: TotalBehaviors ~ 1 + Age_DPF + (1 | ID)
# Data: trialdata
# 
# AIC      BIC   logLik deviance df.resid 
# 3161.7   3168.3  -1577.9   3155.7       63 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -16.006  -4.338  -2.724   4.171  16.727 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 0.06963  0.2639  
# Number of obs: 66, groups:  ID, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 0.336802   0.124994   2.695  0.00705 ** 
#   Age_DPF     0.164533   0.003244  50.717  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# Age_DPF -0.651

#calculated power: 1

#now, in order to test the validity of this model, we need to bootstrap and test confidence intervals

set.seed(78574)

trialdata <- as.data.table(trialdata)

wtbootstrap <- function(df, n) {
  samp <- sample(unique(df$ID), n, replace = TRUE)
  setkey(df, "ID")
  df_sample <- df[J(samp), allow.cartesian = TRUE]
  model <- glmer(formula = TotalBehaviors ~ 1 + Age_DPF + (1|ID), data = df_sample, family=poisson)
  summary <- summary(model)
  Intercept <- summary$coefficients[1,1]
  Age_Coef <- summary$coefficients[2,1]
  pAge <- summary$coefficients[2,4]
  coefs <- rbind(Intercept, Age_Coef, pAge)
  return(coefs)
}

bootstrapresults <- as.data.frame(replicate(1000, wtbootstrap(trialdata, 8)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05)/length(pAge))

Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
ConfidenceIntervals <- cbind(Intercept, Age)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

ggplot(bootstrapresults, aes(x=Intercept))+
  geom_histogram()+
  theme_classic()

#save bootstrapresults for ref:
fwrite(bootstrapresults, file = "bootstrapresults_WT.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_wt.csv", row.names=FALSE)


################################################################################

#old model stuff

set.seed(17000)

model <- glmer(TotalBehaviors ~ 1 + Age_DPF + (1|Brood), data=trialdata, family=poisson)
summary <- summary(model)
output <- capture.output(summary(model), file=NULL,append=F)
output_df <- as.data.frame(output)

#model residuals
hist(summary$residuals)

#model heteroskedacity
res <- summary$residuals
plot(fitted(model), res)
abline(0,0)

#try lme:
model <- lme(TotalBehaviors ~ 1 + Age_DPF, random=~1|Brood, data=trialdata)
summary <- summary(model)
plot(model, resid(.) ~ fitted(.))

trialdata_transform <- trialdata %>%
  mutate(behaviortransform = TotalBehaviors^(1/3))
model <- lme(behaviortransform ~ 1 + Age_DPF, random=~1|Brood, data=trialdata_transform)
summary <- summary(model)
hist(model$residuals)
plot(model, resid(.) ~ fitted(.))

hist(trialdata_transform$behaviortransform)

##### Call Confidence

#note: we need to know the average N across the experiment, therefore:
trialdata_countsummary <- trialdata %>%
  count(Age_DPF)

#average is 3

set.seed(1717171)

#bootstrap formula:
wtbootstrap <- function(df, n) {
  df_sample <- df %>%
    group_by(Age_DPF) %>%
    slice_sample(n=n, replace=TRUE)
  model <- glmer(normalizedcount ~ 1 + Age_DPF + (1|Brood), data=trialdata, family=poisson)
  summary <- summary(model)
  pAge <- summary$coefficients[2,4]
  pvals <- rbind(pAge)
  return(pvals)
}

bootstrapresults <- as.data.frame(replicate(1000, wtbootstrap(trialdata, 3)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- "pAge"

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05)/length(pAge))

#model visualizer

predicteddata <- predict_response(model, terms = c("Age_DPF"))

ggplot(data=predicteddata, aes(x=x, y=predicted))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="grey", alpha=0.5)+
  geom_line(linewidth=1)+
  geom_point(trialdata, mapping=aes(x=Age_DPF, y=TotalBehaviors), size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

###############################################################################

#Individual Behavior Breakdown

#chases
chasesummary <- trialdata %>%
  group_by(Age_DPF) %>%
  summarize(meanchases = mean(Chases, na.rm=TRUE))
lateralsummary <- trialdata %>%
  group_by(Age_DPF) %>%
  summarize(meanlateral = mean(LateralDisplays, na.rm=TRUE))
bitesummary <- trialdata %>%
  group_by(Age_DPF) %>%
  summarize(meanbites = mean(Bites, na.rm=TRUE))
ggplot()+
  geom_point(trialdata, mapping=aes(x=Age_DPF, y=Chases, color='Chases'), size=2, alpha=0.6)+
  geom_point(trialdata, mapping=aes(x=Age_DPF, y=LateralDisplays, color='Lateral Displays'), size=2, alpha=0.6)+
  geom_point(trialdata, mapping=aes(x=Age_DPF, y=Bites, color='Bites'), size=2, alpha=0.6)+
  geom_path(chasesummary, mapping=aes(x=Age_DPF, y=meanchases), color='#EE3B3B', linewidth=1, stat="identity", na.rm=TRUE)+
  geom_path(lateralsummary, mapping=aes(x=Age_DPF, y=meanlateral), color='#8EE5EE', linewidth=1, stat="identity", na.rm=TRUE)+
  geom_path(bitesummary, mapping=aes(x=Age_DPF, y=meanbites), color= '#EEAD0E', linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  scale_color_manual(name='Behavioral Ethogram', 
                     breaks=c('Chases', 'Lateral Displays', 'Bites'), 
                     values=c('Chases' = '#EE3B3B', 'Lateral Displays' = '#8EE5EE', 'Bites'='#EEAD0E'))+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.9),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

## reoing the model breakdown of each behavior individually >:[

behaviorbreakdown <- trialdata[c(3, 7, 9, 10, 11, 12)]
names(behaviorbreakdown)[1] <- "Age"
behaviorbreakdown <- reshape2::melt(behaviorbreakdown, id.vars = c("Age", "ID"), variable.name="Behavior", value.name="Count")

model <- glmer(formula = Count ~ 1 + Age + Behavior + (1|ID), data = behaviorbreakdown, family=poisson)
summary <- summary(model)

predicteddata <- predict_response(model, terms = c("Age", "Behavior"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(behaviorbreakdown, mapping=aes(x=Age, y=Count, color=Behavior, group=Behavior),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Behavior"), fill="none")+
  scale_color_manual(labels=c("Chases", "Bites", "Lateral Displays", "Mutual Chases"), values=c("#EE3B3B", "#EEAD0E", "#8EE5EE", "#952323"))+
  scale_fill_manual(labels=c("Chases", "Bites", "Lateral Displays", "Mutual Chases"), values=c("#EE3B3B", "#EEAD0E", "#8EE5EE", "#952323"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

#without chases:

predicteddata_nochase <- predicteddata %>%
  filter(!group == "Chases")
behaviordata_nochase <- behaviorbreakdown %>%
  filter(!Behavior == "Chases")

ggplot()+
  geom_ribbon(data=predicteddata_nochase, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata_nochase, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(behaviordata_nochase, mapping=aes(x=Age, y=Count, color=Behavior, group=Behavior),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Behavior"), fill="none")+
  scale_color_manual(labels=c("Bites", "Lateral Displays", "Mutual Chases"), values=c("#EEAD0E", "#16b5c6", "#952323"))+
  scale_fill_manual(labels=c("Bites", "Lateral Displays", "Mutual Chases"), values=c("#EEAD0E", "#16b5c6", "#952323"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

predicteddata_lats <- predicteddata %>%
  filter(group == "LateralDisplays")
behaviordata_lats <- behaviorbreakdown %>%
  filter(Behavior == "LateralDisplays")

ggplot()+
  geom_ribbon(data=predicteddata_lats, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata_lats, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(behaviordata_lats, mapping=aes(x=Age, y=Count, color=Behavior, group=Behavior),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#16b5c6"))+
  scale_fill_manual(values=c("#16b5c6"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

#Extended Chases
extendedsummary <- trialdata %>%
  group_by(Age_DPF) %>%
  summarize(meanbehaviors = mean(Extended, na.rm=TRUE))%>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!is.na(meanbehaviors))

ggplot(trialdata, aes(x=Age_DPF, y=Extended))+
  geom_point(size=2, alpha=0.8)+
  geom_line(extendedsummary, mapping=aes(x=Age_DPF, y=meanbehaviors), linewidth=1)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")
