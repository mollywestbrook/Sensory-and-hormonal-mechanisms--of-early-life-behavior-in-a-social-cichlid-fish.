#ARD Stuff

library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)

arddf <- read_csv(here("ARDSummary.csv"), col_names = TRUE)

behaviorsummary <- arddf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))
ggplot(arddf, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.6, position = position_dodge(width = 0.9))+
  geom_line(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#004D40", "#D81B60"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#stats

arddf <- arddf %>%
  arrange(Experiment, Tank)
arddf$ID <- cumsum(!duplicated(arddf[c(1,2)]))

ardcount <- arddf %>%
  group_by(Treatment) %>%
  summarize(count = length(unique(ID)))

#GLMM
model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|ID), data = arddf, family=poisson)
summary <- summary(model)

#output
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1 | ID)
# Data: arddf
# 
# AIC      BIC   logLik deviance df.resid 
# 2540.0   2549.9  -1265.0   2530.0       49 
# 
# Scaled residuals: 
#   Min     1Q Median     3Q    Max 
# -8.366 -4.601 -1.552  1.602 24.627 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ID     (Intercept) 2.122    1.457   
# Number of obs: 54, groups:  ID, 12
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        1.945770   0.633247   3.073 0.002121 ** 
#   Age                0.024181   0.007225   3.347 0.000817 ***
#   TreatmentHomo     -2.491102   0.893566  -2.788 0.005306 ** 
#   Age:TreatmentHomo  0.130514   0.010525  12.401  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) Age    TrtmnH
# Age         -0.230              
# TreatmentHm -0.702  0.166       
# Ag:TrtmntHm  0.160 -0.686 -0.247

#try a quick bootstrap with our sad dataset just to test:
set.seed(54512)

#real bootstrap
set.seed(302817)

arddf <- as.data.table(arddf)

ardbootstrap <- function(df, n) {
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
      return(c(NA, NA, NA, NA, NA, NA, NA))
    }
  )
}

bootstrapresults <- as.data.frame(replicate(1000, ardbootstrap(arddf, 12)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)
names(bootstrapresults) <- c("Intercept", "Age_Coef", "Homo_Coef", "HomoAge_Coef",
                             "pAge", "pHomo", "pHomoAge")

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

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Intercept), fill="#004D40", alpha=0.5, binwidth=0.25)+
  geom_histogram(bootstrapresults, mapping=aes(x=Homo_Coef), fill="#D81B60", alpha=0.5, binwidth=0.25)+
  theme_classic()+
  coord_cartesian(xlim=c(-25, 25))

ggplot()+
  geom_histogram(bootstrapresults, mapping=aes(x=Age_Coef), fill="#004D40", alpha=0.5, binwidth=0.01)+
  geom_histogram(bootstrapresults, mapping=aes(x=HomoAge_Coef), fill="#D81B60", alpha=0.5, binwidth=0.01)+
  theme_classic()+
  coord_cartesian(xlim=c(-1, 2))

fwrite(bootstrappower, file = "bootstrappower_ard.csv", row.names=FALSE)
fwrite(ConfidenceIntervals, file = "confidenceintervals_ard.csv", row.names=FALSE)


#model visualization:
#step one: predict effects:
predicteddata <- predict_response(model, terms = c("Age", "Treatment"))

ggplot()+
  geom_ribbon(data=predicteddata, mapping=aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  geom_line(data=predicteddata, mapping=aes(x=x, y=predicted, group=group, color=group), linewidth=1.5)+
  geom_point(arddf, mapping=aes(x=Age, y=Behaviors, color=Treatment, group=Treatment),
             size=3, alpha=0.6, position = position_dodge(width = 0.9))+
  theme_classic()+
  theme(text = element_text(size = 18), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+
  guides(color=guide_legend(title="Treatment"), fill="none")+
  scale_color_manual(values=c("#004D40", "#D81B60"))+
  scale_fill_manual(values=c("#004D40", "#D81B60"))+
  ylab("Aggressive Behaviors/30 Minutes")+
  xlab("Age (Days Post Fertilization)")