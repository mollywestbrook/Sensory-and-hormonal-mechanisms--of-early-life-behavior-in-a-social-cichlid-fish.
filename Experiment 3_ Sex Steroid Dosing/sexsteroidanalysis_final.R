#Sex Steroid Final

#libraries
library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)


#Total Behaviors -- No Model
hormonedf <- read_csv(here("SexSteroidComplete.csv"), col_names = TRUE)

#for this I will filter for just our ages of interest:
hormonedf <- hormonedf %>%
  filter(Age %in% c(15:30))

behaviorsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Normalized_Behaviors, na.rm=TRUE))
ggplot(hormonedf, aes(x=Age, y=Normalized_Behaviors, color=as.factor(Treatment), group=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
 theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#with individual IDs plotted together
ggplot(hormonedf, aes(x=Age, y=Normalized_Behaviors, color=as.factor(Treatment), group=as.factor(ID)))+
  geom_point(size=2, alpha=0.6)+
  geom_path(linewidth=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#Model:

#add an ID column:
hormonedf <- hormonedf %>%
  arrange(Experiment, Tank)
hormonedf$ID <- cumsum(!duplicated(hormonedf[c(1,3)])) #29 total individuals for bootstrap

hormonedfcount <- hormonedf %>%
  group_by(Treatment) %>%
  summarize(count = length(unique(ID)))

model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), data = hormonedf, family=poisson)
summary <- summary(model)

#output:
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1 | ID)
# Data: hormonedf
# 
# AIC      BIC   logLik deviance df.resid 
# 7284.1   7304.9  -3635.1   7270.1      137 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -10.843  -4.404  -2.298   2.526  23.984 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 0.6352   0.797   
# Number of obs: 144, groups:  ID, 29
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     -0.411689   0.324238  -1.270  0.20419    
# Age              0.169181   0.006455  26.211  < 2e-16 ***
#   TreatmentEE      1.268779   0.418250   3.034  0.00242 ** 
#   TreatmentKT      0.339602   0.427849   0.794  0.42735    
# Age:TreatmentEE -0.050991   0.007765  -6.567 5.15e-11 ***
#   Age:TreatmentKT -0.002801   0.007816  -0.358  0.72004    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) Age    TrtmEE TrtmKT Ag:TEE
# Age         -0.480                            
# TreatmentEE -0.775  0.372                     
# TreatmentKT -0.757  0.364  0.588              
# Ag:TrtmntEE  0.399 -0.831 -0.444 -0.302       
# Ag:TrtmntKT  0.397 -0.826 -0.307 -0.439  0.687

#bootstrap
set.seed(56452)

hormonedf <- as.data.table(hormonedf)

sexsteroidbootstrap <- function(df, n) {
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
        EE_Coef = summary$coefficients[3,1],
        KT_Coef = summary$coefficients[4,1],
        EEAge_Coef = summary$coefficients[5,1],
        KTAge_Coef = summary$coefficients[6,1],
        pAge = summary$coefficients[2,4],
        pEE = summary$coefficients[3,4],
        pKT = summary$coefficients[4,4],
        pEEAge = summary$coefficients[5,4],
        pKTAge = summary$coefficients[6,4]
      )
    }, error = \(e) {
      print(e)
      return(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, sexsteroidbootstrap(hormonedf, 29)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "EE_Coef", "KT_Coef",
                             "EEAge_Coef", "KTAge_Coef", "pAge", "pEE", "pKT",
                             "pEEAge", "pKTAge")

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05)/length(pAge),
            pEE = sum(pEE < 0.05)/length(pEE),
            pKT = sum(pKT < 0.05)/length(pKT),
            pEEAge = sum(pEEAge < 0.05)/length(pEEAge),
            pKTAge = sum(pKTAge < 0.05)/length(pKTAge))

#confidence intervals: 
Intercept <- quantile(bootstrapresults$Intercept, prob=c(.025,.5,.975), na.rm=T)
EE <- quantile(bootstrapresults$EE_Coef, prob=c(.025,.5,.975), na.rm=T)
KT <- quantile(bootstrapresults$KT_Coef, prob=c(.025,.5,.975), na.rm=T)
Age <- quantile(bootstrapresults$Age_Coef, prob=c(.025,.5,.975), na.rm=T)
EEAge <- quantile(bootstrapresults$EEAge_Coef, prob=c(.025,.5,.975), na.rm=T)
KTAge <- quantile(bootstrapresults$KTAge_Coef, prob=c(.025,.5,.975), na.rm=T)

ConfidenceIntervals <- cbind(Intercept, EE, KT, Age, EEAge, KTAge)
ConfidenceIntervals <- data.frame(ConfidenceIntervals)

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="#004D40", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=EE_Coef), fill="#D81B60", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=KT_Coef), fill="#1E88E5", alpha=0.5, binwidth=0.25)+
  theme_classic()

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="#004D40", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=EEAge_Coef), fill="#D81B60", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=KTAge_Coef), fill="#1E88E5", alpha=0.5, binwidth=0.01)+
  theme_classic()

#trying to plot our function over the plot:

###Using a custom function
# 
# plotDMSO <- function(y) {
#   exp(summary$coefficients[1,1] + summary$coefficients[2,1]*y)
# }
# plotEE <- function(y) {
#   exp(summary$coefficients[1,1] + summary$coefficients[3,1] + (summary$coefficients[2,1] + summary$coefficients[5,1])*y)
# }
# plotKT <- function(y) {
#   exp(summary$coefficients[1,1] + summary$coefficients[4,1] + (summary$coefficients[2,1] + summary$coefficients[6,1])*y)
# }
# DMSOsubset <- hormonedf %>%
#   arrange(Age) %>%
#   group_by(Treatment) %>%
#   filter(Treatment == "DMSO") %>%
#   mutate(SE = sd(Behaviors, na.rm=T)/length(unique(ID))) %>%
#   mutate(SEhigh = Behaviors + SE) %>%
#   mutate(SElow = Behaviors - SE) %>%
#   mutate(SEhighfun = plotDMSO(SEhigh)) %>%
#   mutate(fun = plotDMSO(Behaviors)) %>%
#   mutate(SElowfun = plotDMSO(SElow))
# EEsubset <- hormonedf %>%
#   group_by(Treatment) %>%
#   filter(Treatment == "EE")
# KTsubset <- hormonedf %>%
#   group_by(Treatment) %>%
#   filter(Treatment == "KT")
# 
# ggplot(hormonedf, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.6, position = position_dodge(width = 0.9))+
#   geom_function(data=DMSOsubset, fun=plotDMSO, color="#004D40", linewidth=1)+
#   geom_function(data=EEsubset, fun=plotEE, color="#D81B60", linewidth=1)+
#   geom_function(data=KTsubset, fun=plotKT, color="#1E88E5", linewidth=1)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
#   ylab("Aggressive Behavior Count")+
#   xlab("Age (Days Post Fertilization)")

#using ggeffect (we're going with this one)

#step one: predict effects:
predicteddata <- predict_response(model, terms = c("Age", "Treatment"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(hormonedf, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

#DMSO curve only
predicteddata_DMSOonly <- predicteddata %>%
  filter(group == "DMSO")

ggplot()+
  geom_ribbon(data=predicteddata_DMSOonly, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata_DMSOonly, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(hormonedf, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  scale_fill_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")

#save bootstrapresults for ref:
fwrite(bootstrappower, file = "bootstrappower_sexsteroids.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_sexsteroids.csv", row.names=FALSE)

fwrite(summary$coefficients, file = "coefficients.csv", row.names=FALSE)

###############################################################################

#Other behaviors in isolation:

#lateral displays
lateralsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(lateralmean = mean(Lateral, na.rm=TRUE))
ggplot(hormonedf, aes(x=Age, y=Lateral, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(lateralsummary, mapping=aes(x=Age, y=lateralmean, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.07, 0.9),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Lateral Display Count")+
  xlab("Age (Days Post Fertilization)")

#extended bouts
extendedsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(extendmean = mean(Extended, na.rm=TRUE))
ggplot(hormonedf, aes(x=Age, y=Extended, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(extendedsummary, mapping=aes(x=Age, y=extendmean, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.07, 0.9),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Extended Bout Count")+
  xlab("Age (Days Post Fertilization)")