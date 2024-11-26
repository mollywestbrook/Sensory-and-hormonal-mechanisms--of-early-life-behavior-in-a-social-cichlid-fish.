#Calculating day of aggression emergence

library(tidyverse)
library(here)
library(lme4)
library(data.table)
library(ggeffects)

trialdata <- read_csv(here("AllExperimentsPooled_summary.csv"), col_names = TRUE)

#step one: attach ID

trialdata <- trialdata %>%
  arrange(Experiment, Brood, Tank)
trialdata$ID <- cumsum(!duplicated(trialdata[c(1,2,3)]))

trialdata_filter <- trialdata %>%
  arrange(ID) %>%
  filter(!is.na(Behaviors)) %>%
  filter(!Behaviors == 0) %>%
  group_by(ID) %>%
  slice(1) %>%
  arrange(Age)
  
ggplot(trialdata_filter, aes(x=Age))+
  geom_histogram(binwidth=1, fill="#242424")+
  geom_vline(xintercept=mean(trialdata_filter$Age), color="black", linewidth=2)+
  geom_vline(xintercept=median(trialdata_filter$Age), color="red", linewidth=2)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Count")+
  xlab("Age(DPF) Aggression First Expressed")

#average age = 17.63529
#sd: 2.943302
#median age = 17

  