library(tidyverse)
library(here)

#input data

homedata<- read_csv(here("earlylifeaggression_homeenvironment.csv"), col_names = TRUE)
noveldata<- read_csv(here("earlylifeaggression_novelenvironment.csv"), col_names = TRUE)

#novel environment

noveldatasummary <- noveldata %>%
  group_by(Age) %>%
  summarize(meanbehavior = mean(Count, na.rm=TRUE))

homedatasummary <- homedata %>%
  group_by(Age) %>%
  summarize(meanbehavior = mean(Count, na.rm=TRUE))

homedatasummary$Group <- c("Home")
noveldatasummary$Group <- c("Novel")
combineddatasummary <- rbind(homedatasummary, noveldatasummary)

noveldata$Group <- c("Novel")
homedata$Group <- c("Home")
homedata <- homedata[-c(1)]
combineddata <- rbind(homedata, noveldata)

ggplot(combineddatasummary, aes(x=Age, y=meanbehavior, color=Group))+
  geom_point(size=4, alpha=0.4)+
  geom_path(size=0.5)+
  geom_point(combineddata, mapping=aes(x=Age, y=Count, color=Group), size=2)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (DPF)")

#####

#counts for hormone dosing pt2

hormonedata <- read_csv(here("HormoneResultsPt2.csv"), col_names = TRUE)
hormonedata <- read_csv(here("HormoneResultsPart2Updated.csv"), col_names = TRUE)
hormonedata_lateralonly <- read_csv(here("HormoneResultsPart2Updated.csv"), col_names = TRUE)


##All Behaviors:
#everything together
ggplot(hormonedata, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=3, alpha=0.6)+
#  geom_path(size=0.5)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (DPF)")

#with groups that need to be excluded
hormonedata_goodexclusions <- hormonedata %>%
  group_by(Tank, Experiment, Treatment) %>%
  filter(!(Tank %in% c(4,5) & Experiment %in% c(2))) %>%
  filter(!str_detect(Treatment, "11KT"))

hormonesummary <- hormonedata_goodexclusions %>%
  group_by(Treatment, Age_DPA) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(hormonedata_goodexclusions, aes(x=Age_DPA, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(hormonesummary, mapping=aes(x=Age_DPA, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_color_discrete(name="Treatment:",
                    labels=c("DMSO", "5-7dpf EE", "5-12dpf EE", "10-12dpf EE"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Yolk-Sac Absorbance)")

##Just the Lateral Displays
hormonedata_lateralonly <- read_csv(here("HormoneResultsCombined_borisevents.csv"), col_names = TRUE)

lateraldisplays_excluded <- hormonedata_lateralonly %>%
  group_by(Tank, Experiment, Treatment) %>%
  filter(!(Tank %in% c(4,5) & Experiment %in% c(2))) %>%
  filter(!str_detect(Treatment, "11KT"))

lateralsummary <- lateraldisplays_excluded %>%
  group_by(Treatment, Age_DPA) %>%
  summarize(meandisplays = mean(LateralCount, na.rm=TRUE))

ggplot(lateraldisplays_excluded, aes(x=Age_DPA, y=LateralCount, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(lateralsummary, mapping=aes(x=Age_DPA, y=meandisplays, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_color_discrete(name="Treatment:",
                       labels=c("DMSO", "5-7dpf EE", "5-12dpf EE", "10-12dpf EE"))+
  ylab("Lateral Display Count")+
  xlab("Age (Days Post Yolk-Sa Absorbance)")

##Androgen receptor knockouts

#all behaviors:
ara_allbehaviors <- read_csv(here("ARA_allbehaviors.csv"), col_names = TRUE)

behaviorsummary <- ara_allbehaviors %>%
  group_by(Genotype, Age_DPA) %>%
  summarize(meanbehaviors = mean(Count, na.rm=TRUE))

ggplot(ara_allbehaviors, aes(x=Age_DPA, y=Count, color=as.factor(Genotype)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age_DPA, y=meanbehaviors, color=as.factor(Genotype)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_color_manual(name="Genotype:",
                       values=c("#4D8C8E", "#5C399A", "#969696"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Yolk-Sac Absorbance)")



#lateral displays:
ara_lateralbehaviors <- read_csv(here("ara_lateraldisplaysonly.csv"), col_names = TRUE)

lateralsummary <- ara_lateralbehaviors %>%
  group_by(Genotype, Age_DPA) %>%
  summarize(meanbehaviors = mean(LateralDisplayCount, na.rm=TRUE))

ggplot(ara_lateralbehaviors, aes(x=Age_DPA, y=LateralDisplayCount, color=as.factor(Genotype)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(lateralsummary, mapping=aes(x=Age_DPA, y=meanbehaviors, color=as.factor(Genotype)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_color_manual(name="Genotype:",
                     values=c("#4D8C8E", "#5C399A", "#969696"))+
  ylab("Lateral Display Count")+
  xlab("Age (Days Post Yolk-Sac Absorbance)")

#####

#General WT Figure

wtdata <- read_csv(here("WTAggressionOnly.csv"), col_names = TRUE)

#a bit of filtering -- exclude groups 2 and 5

wtdata_filtered <- wtdata %>%
  group_by(Brood) %>%
  filter(!Brood %in% c(2,5))

wtdata_summary <- wtdata_filtered %>%
  group_by(Age_DPA) %>%
  summarize(meanbehaviors = mean(Count, na.rm=TRUE))

ggplot(wtdata_filtered, aes(x=Age_DPA, y=Count))+
  geom_point(size=2, alpha=0.25)+
  geom_path(wtdata_summary, mapping=aes(x=Age_DPA, y=meanbehaviors),
            linewidth=0.8)+
  theme_classic()+
  theme(text = element_text(size = 15))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Yolk-Sac Absorbtion)")

### 

##11kt Stuff
e11kt <- read_csv(here("11KTExp.csv"), col_names = TRUE)

e11ktsummary <- e11kt %>%
  group_by(KTExposure, Age) %>%
  summarise(meanbehavior= mean(Count, na.rm=TRUE))

ggplot(e11kt, aes(x=Age, y=Count, color=as.factor(KTExposure)))+
  geom_point(size=2, alpha=0.25)+
  geom_path(e11ktsummary, mapping=aes(x=Age, 
                                      y=meanbehavior, 
                                      color=as.factor(KTExposure)),
            linewidth=0.8)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_color_discrete(name="Treatment:")+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")


####################

#ARD

ard_allbehaviors <- read_csv(here("ARDfull.csv"), col_names = TRUE)

ard_allbehaviors <- ard_allbehaviors %>%
  filter(VidLegnth == 30)

behaviorsummary <- ard_allbehaviors %>%
  group_by(Genotype, Age_DPF) %>%
  summarize(meanbehaviors = mean(Count, na.rm=TRUE))

ggplot(ard_allbehaviors, aes(x=Age_DPF, y=Count, color=Genotype))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age_DPF, y=meanbehaviors, color=as.factor(Genotype)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

####################

#Hormone dose 3

hormonedose3 <- read_csv(here("hormonedosing3.csv"), col_names = TRUE)

behaviorsummary <- hormonedose3 %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(hormonedose3, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

####################

#Hormone dose 4

hormonedose4 <- read_csv(here("HormoneDosing4.csv"), col_names = TRUE)

hormonedose4 <- hormonedose4 %>%
  filter(!Treatment == "T7")

behaviorsummary <- hormonedose4 %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(hormonedose4, aes(x=Age, y=Behaviors, color=as.factor(Treatment)))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            size=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#######################################################

#CNG

cng <- read_csv(here("CNG.csv"), col_names = TRUE)

behaviorsummary <- cng %>%
  group_by(Treatment, Age_DPF) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(cng, aes(x=Age_DPF, y=Behaviors, color=Treatment))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age_DPF, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

############3

#thyroidhormone

thyroid <- read_csv(here("thyroidhormonestuff.csv"), col_names = TRUE)

behaviorsummary <- thyroid %>%
  group_by(Treatment, Age) %>%
  summarize(meanbehaviors = mean(Behaviors, na.rm=TRUE))

ggplot(thyroid, aes(x=Age, y=Behaviors, color=Treatment))+
  geom_point(size=2, alpha=0.4)+
  geom_path(behaviorsummary, mapping=aes(x=Age, y=meanbehaviors, color=as.factor(Treatment)),
            linewidth=1, stat="identity", na.rm=TRUE)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

