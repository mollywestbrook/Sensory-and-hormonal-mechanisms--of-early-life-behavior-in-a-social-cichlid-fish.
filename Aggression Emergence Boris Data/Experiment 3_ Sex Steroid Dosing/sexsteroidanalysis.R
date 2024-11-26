## Fresh Hormone Dosing Analysis

#libraries
library(tidyverse)
library(here)

###step one: organize full boris data into summarized counts

boriscounts <- read_csv(here("HormoneDosing1-6.csv"), col_names = TRUE)

borissummary <- boriscounts %>%
  group_by(ObservationID, Behavior) %>%
  summarize(BehaviorTotal = length(Behavior))
#this will also go faster if we cast these observations into wide format
borissummarywide <- spread(borissummary, Behavior, BehaviorTotal)
borissummarywide[is.na(borissummarywide)] <- 0
#write that out to combine with our summary
write_csv(borissummarywide, "HormoneDoseBorisSummary.csv")

##############################################################

#step 2: analyze complete dataset
hormonedf <- read_csv(here("SexSteroidComplete.csv"), col_names = TRUE)

behaviorsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Normalized_Behaviors, na.rm=TRUE))

ggplot(hormonedf, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#big  
ggplot(hormonedf, aes(x=Age, y=Normalized_Behaviors, color=as.factor(Treatment)))+
    geom_point(size=6, alpha=0.4)+
    geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
              size=3, stat="identity", na.rm=TRUE)+
    theme_classic()+
    theme(text = element_text(size = 28), 
          legend.position = c(0.1, 0.85),
          legend.text = element_text(size=15),
          legend.title = element_text(size=16))+
    guides(color=guide_legend(title="Treatment"))+
    scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behaviors/5 min")+
    xlab("Age (Days Post Fertilization)")
  
  

#lateral displays
lateralsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(lateralmean = mean(Lateral, na.rm=TRUE))
ggplot(hormonedf, aes(x=Age, y=Lateral, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(lateralsummary, mapping=aes(x=Age, y=lateralmean, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
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
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Extended Bout Count")+
  xlab("Age (Days Post Fertilization)")

# #DecemberAnomalystuff
# hormonedf_noanomaly <- read_csv(here("HormoneDosingSummaryAnomalyRemoved.csv"), col_names = TRUE)
# 
# behaviorsummary_noanomaly <- hormonedf_noanomaly %>%
#   group_by(Treatment, Age) %>%
#   summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))
# ggplot(hormonedf_noanomaly, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#  # geom_path(behaviorsummary_noanomaly, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
# #            size=1, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Aggressive Behavior Count")+
#   xlab("Age (Days Post Fertilization)")
# 
# #and again, lateral and extended:
# lateralsummary_noanom <- hormonedf_noanomaly %>%
#   group_by(Treatment, Age) %>%
#   summarize(lateralmean = mean(Lateral, na.rm=TRUE))
# ggplot(hormonedf_noanomaly, aes(x=Age, y=Lateral, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#   geom_path(lateralsummary_noanom, mapping=aes(x=Age, y=lateralmean, color=as.factor(Treatment)),
#             size=1, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Lateral Display Count")+
#   xlab("Age (Days Post Fertilization)")
# 
# extendedsummary_noanom <- hormonedf_noanomaly %>%
#   group_by(Treatment, Age) %>%
#   summarize(extendmean = mean(Extended, na.rm=TRUE))
# ggplot(hormonedf_noanomaly, aes(x=Age, y=Extended, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#   geom_path(extendedsummary_noanom, mapping=aes(x=Age, y=extendmean, color=as.factor(Treatment)),
#             size=1, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Extended Bout Count")+
#   xlab("Age (Days Post Fertilization)")

#and now we do raw behavior counts with all videos include w/5 fish
#including those below 30min
# 
# #full dataset first:
# hormonedose30 <- read_csv(here("HormoneDosing_Sub30VidsIncluded.csv"), col_names = TRUE)
# 
# behaviorsummary <- hormonedose30 %>%
#   group_by(Treatment, Age) %>%
#   summarize(meanbehaviors = mean(BehaviorsNormalizedToVidLength, na.rm=TRUE))
# ggplot(hormonedose30, aes(x=Age, y=BehaviorsNormalizedToVidLength, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#   geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
#             size=1, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Aggressive Behaviors/Minute")+
#   xlab("Age (Days Post Fertilization)")
# 
# #anomaly removed
# hormonedose30_noanomaly <- read_csv(here("HormoneDosing_Sub30VidsIncluded_noanomaly.csv"), col_names = TRUE)
# 
# behaviorsummary_noanomaly <- hormonedose30_noanomaly %>%
#   group_by(Treatment, Age) %>%
#   summarize(meanbehaviors = mean(BehaviorsNormalizedToVidLength, na.rm=TRUE))
# ggplot(hormonedose30_noanomaly, aes(x=Age, y=BehaviorsNormalizedToVidLength, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#   geom_path(behaviorsummary_noanomaly, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
#             size=1, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Aggressive Behaviors/Minute")+
#   xlab("Age (Days Post Fertilization)")

#########################################################

#pulling the most robust data:

#hormonedf <- read_csv(here("HormoneDosingSummaryAllData.csv"), col_names = TRUE)

behaviorsummary <- hormonedf %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

hormonedf$Treatment <- as.factor(hormonedf$Treatment)
model <- lm(formula = Behaviors ~  1 + Age + Treatment + Age:Treatment, data = hormonedf)
summary <- summary(model)
output <- capture.output(summary(model), file=NULL,append=F)
output_df <- as.data.frame(output)
write.csv(output_df, file="modelsummary.csv")

ggplot(hormonedf, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age:TreatmentEE"], color="green", linewidth=1)+
  geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age:TreatmentKT"], color="blue", linewidth=1)+
  geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age"], color="red", linewidth=1)+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=0.8, alpha=0.4, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

# #ok now do it for the data without the anomaly
# hormonedf_noanomaly <- read_csv(here("HormoneDosingSummaryAnomalyRemoved.csv"), col_names = TRUE)
# 
# hormonedf_noanomaly$Treatment <- as.factor(hormonedf_noanomaly$Treatment)
# model <- lm(formula = Behaviors ~  Age + Treatment + Age:Treatment, data = hormonedf_noanomaly)
# summary <- summary(model)
# summary$coefficients
# output <- capture.output(summary(model), file=NULL,append=F)
# output_df <- as.data.frame(output)
# write.csv(output_df, file="modelsummary_noanomaly.csv")
# 
# behaviorsummary_noanomaly <- hormonedf_noanomaly %>%
#   group_by(Treatment, Age) %>%
#   summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))
# ggplot(hormonedf_noanomaly, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
#   geom_point(size=2, alpha=0.4)+
#   geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age:TreatmentEE"], color="green", linewidth=1)+
#   geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age:TreatmentKT"], color="blue", linewidth=1)+
#   geom_abline(intercept = coef(model)["(Intercept)"], slope = coef(model)["Age"], color="red", linewidth=1)+
#   geom_path(behaviorsummary_noanomaly, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
#             linewidth=0.8, alpha=0.4, stat="identity", na.rm=TRUE)+
#   theme_classic()+
#   theme(text = element_text(size = 15), 
#         legend.position = c(0.1, 0.85),
#         legend.text = element_text(size=10),
#         legend.title = element_text(size=11))+
#   guides(color=guide_legend(title="Treatment"))+
#   ylab("Aggressive Behavior Count")+
#   xlab("Age (Days Post Fertilization)")


#### power analysis

#sample from existing data with n=goal# (which can be adjusted)
#bootstrap 1000 times
#run lm extract p-values for age:treatmentEE and age:treatmentKT
#what percentage of time do I reject the null? If it's >95% drop sample size, if not increase

#save model summary as object
#summary$coefficients[x,y] matrix of model coefficients --> select pvalues of interest

#use sub 30 data with anomaly:

#hormonedose30 <- read_csv(here("HormoneDosing_Sub30VidsIncluded.csv"), col_names = TRUE)

hormonedf$Treatment <- as.factor(hormonedf$Treatment)
model <- lm(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment, data = hormonedf)
summary <- summary(model)
summary

bootstrapdf <- as.data.frame(hormonedf[c(4,6,7)])
# bootstrapdf <- bootstrapdf %>%
#   filter(Age %in% c(15, 17, 20, 22, 23, 24, 27, 30))

bootstrapdf_sample <- bootstrapdf %>%
  group_by(Treatment, Age) %>%
  slice_sample(n=10, replace=TRUE)

hormonedosagebootstrap <- function(df, n) {
  df_sample <- df %>%
    group_by(Treatment, Age) %>%
    slice_sample(n=n, replace=TRUE)
  model <- lm(formula = Behaviors ~  Age + Treatment + Age:Treatment, data = df_sample)
  summary <- summary(model)
  pAge <- summary$coefficients[2,4]
  pEE <- summary$coefficients[5,4]
  pKT <- summary$coefficients[6,4]
  pvals <- rbind(pAge, pEE, pKT)
  return(pvals)
}

set.seed(95814583) #1-5
set.seed(374308) #1-6

bootstrapresults <- as.data.frame(replicate(1000, hormonedosagebootstrap(bootstrapdf, 8)))
bootstrapresults <- t(bootstrapresults)
bootstrapresults <- as.data.frame(bootstrapresults)

bootstrappower <- bootstrapresults %>%
  summarize(pAge = sum(pAge < 0.05)/length(pAge),
            pEE = sum(pEE < 0.05)/length(pEE),
            pKT = sum(pKT < 0.05)/length(pKT))

#old model stuff

# #try something real quick
# hormonedf$startofexperiment <- hormonedf$Age - 14
#   
# set.seed(19172)
# 
# model <- glmer(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment + (1|Experiment) + (1|ID), data = hormonedf, family=poisson)
# summary <- summary(model)
# output <- capture.output(summary(model), file=NULL,append=F)
# output_df <- as.data.frame(output)
# 
# #model residuals
# model$residuals
# hist(model$residuals)
# 
# #model heteroskedacity
# res <- model$residuals
# plot(fitted(model), res)
# abline(0,0)
# 
# #bootstrap
# model <- glm(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment, data = hormonedf, family=poisson)
# summary <- summary(model)
# summary
# model <- glm(formula = Behaviors ~ 1 + Age + Treatment, data = hormonedf, family=poisson)
# summary <- summary(model)
# summary
# #formula tester
# hormonedf$Treatment <- as.factor(hormonedf$Treatment)
# model <- lm(formula = Behaviors ~ 1 + Age + Treatment + Age:Treatment, data = hormonedf)
# summary <- summary(model)
# summary
# bootstrapdf <- as.data.frame(hormonedf[c(4,6,8)])
# bootstrapdf_sample <- bootstrapdf %>%
#   group_by(Treatment, StartofExp) %>%
#   slice_sample(n=5, replace=TRUE)
# 
# #actual formula:
# hormonedosagebootstrap <- function(df, n) {
#   df_sample <- df %>%
#     group_by(Treatment, StartofExp) %>%
#     slice_sample(n=n, replace=TRUE)
#   model <- lm(formula = Behaviors ~  StartofExp + Treatment + StartofExp:Treatment, data = df_sample)
#   summary <- summary(model)
#   pAge <- summary$coefficients[2,4]
#   pEE <- summary$coefficients[5,4]
#   pKT <- summary$coefficients[6,4]
#   pvals <- rbind(pAge, pEE, pKT)
#   return(pvals)
# }
# 
# set.seed(4372423) #complete
# 
# set.seed(95814583) #1-5
# set.seed(374308) #1-6
# 
# bootstrapresults <- as.data.frame(replicate(1000, hormonedosagebootstrap(bootstrapdf, 17)))
# bootstrapresults <- t(bootstrapresults)
# bootstrapresults <- as.data.frame(bootstrapresults)
# 
# bootstrappower <- bootstrapresults %>%
#   summarize(pAge = sum(pAge < 0.05)/length(pAge),
#             pEE = sum(pEE < 0.05)/length(pEE),
#             pKT = sum(pKT < 0.05)/length(pKT))


