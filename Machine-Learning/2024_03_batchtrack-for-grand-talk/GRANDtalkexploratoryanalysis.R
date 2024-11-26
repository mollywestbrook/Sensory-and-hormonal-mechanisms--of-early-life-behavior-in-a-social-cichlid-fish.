### Exploratory Analaysis for Hormone Dosing Data

#library Import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#Borrowing things from featureextraction and datapipeline as needed
#as these datasets have known problems and will need to be reanalyzed anyway

trialdata <- read.csv(here("formatted_01-21-2024_piz1_T1.csv"))

#to see if these are even remotely close to accurate let's do a visualization
#of the velocity distributions of each fish
#this will help me determine the threshold to set velcoity to NA for now

velocitycalc <- function(a){
  velocity = a - lag(a)
}
velocityframe <- sapply(trialdata_tracks, velocitycalc)
quantile(velocityframe, probs=c(0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999), na.rm=T)


ggplot(featuredf, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_violin()+
  theme_classic()+
  scale_y_continuous(trans='log10') 

ggplot(featuredf, aes(x=velocity_mmpersec))+
  geom_histogram(binwidth=10)+
  theme_classic()

quantile(featuredf$velocity_mmpersec, probs=c(0.1, 0.25, 0.5, 0.75, 0.99), na.rm = T)

#Went back and forth between distributions...going to do a threshold of 500

#here we'll strap on age and treatment:
featuredf$Treatment <- "KT"
featuredf$Age <- "Post-Aggression"
featuredf$ID <- "11"
featuredf11 <- featuredf

rm(featuredf, head_df, nearestneighbordf, NNDF, trialdata)

#note: seconds are 25fps, replace with 30fps

bigoldf <- rbind(featuredf1, featuredf2, featuredf3, featuredf4, featuredf5,
                 featuredf6, featuredf7, featuredf8, featuredf9, featuredf10, 
                 featuredf11)

write.csv(bigoldf, "alldatasetscombined.csv")


########################################################################

#going to do a very straightforward thing, which is just generate the velocity 
#violin plots for each group
ggplot(featuredf, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_violin()+
  theme_classic()+
  scale_y_continuous(trans='log10')


##########################################################################

#so the best way to make these I think. Is going to be a bootstrap
#so we'll combine whatever replicates we have together
#sample 53000 points from the combined dataset for each age and treatment
#and plot those...as much as I hate to do it
#I may just make a massive dataframe for these analyses -_-


#the coordinates to get rid of for this dataset:
x = 360-390
y = 320-350

##############################################################################

#now onto my chase figure...this will be taken from chaseanalysis.R in some regard

#the last time I did this I filtered between the top two velocities...
#based on our existing definition, a more accurate way of doing this
#would probably be to filter for the maximum acceleration
#and then minimum acceleration of the relative nearest neighbor of the max -_-]

#do it on the whole damn thing!
bdf <- read.csv(here("alldatasetscombined_velocityaccelerationfilter.csv"))

bdf <- bdf %>%
  group_by(ID, track) %>%
  mutate(velocity_mmpersec_filter = case_when(velocity_mmpersec > 600 ~ NA,
                                              velocity_mmpersec < -600 ~ NA,
                                              .default = velocity_mmpersec))

bdf <- bdf %>%
  group_by(ID, track) %>%
  mutate(acceleration_mmpsps_filter = case_when(acceleration_mmpsps > 600 ~ NA,
                                                acceleration_mmpsps < -600 ~ NA,
                                              .default = acceleration_mmpsps))

write.csv(bigoldf_behaviorannotated, "alldatasetscombined_filterandbehaviorannotated.csv")

bdf <- bigoldf

#######################################################

#chase filter!
bdf <- read.csv(here("alldatasetscombined_velocityaccelerationfilter.csv"))

df_11 <- bigoldf %>%
  filter(ID == 11)

#stuff for saving these
annotatedbehavior_11 <- annotatedbehavior
bdf_11 <- bdf

#################################################################

chasenn <- bdf_11 %>%
  group_by(BehaviorID)%>%
  slice_min(order_by=relativedist, n=1) %>%
  mutate(selecttrack = 1) %>%
  select(BehaviorID, frame, track, selecttrack) %>%
  arrange(BehaviorID, frame)

bdf_labelled <- merge(bdf_11, chasenn, by=c("BehaviorID", "track"), all.x=T)
bdf_labelled <- bdf_labelled[-c(24)]
bdf_labelled <- bdf_labelled %>% rename_at(3, ~'frame')

chasemaxs <- bdf_labelled %>%
  group_by(BehaviorID, track) %>%
  filter(selecttrack == 1) %>%
  slice_max(order_by=acceleration_mmpsps, n=1) 
chasemins <- bdf_labelled %>%
  group_by(BehaviorID, track) %>%
  filter(selecttrack == 1) %>%
  slice_min(order_by=acceleration_mmpsps, n=1)

chaselengths <- rbind(chasemaxs, chasemins)
chaselengths <- chaselengths %>%
  arrange(frame)
chaselengths <- chaselengths %>%
  group_by(BehaviorID, Treatment, Age) %>%
  summarize(chaselength = max(frame)-min(frame)) %>%
  mutate(chaselength_sec = chaselength/30)

#idk let's just see how this goes:
chaselengths_11 <- chaselengths

rm(chasemaxs, chasemins, chasenn, bdf_labelled, chaselengths)

######################

#chase length figure stuff

chaselengths_1 <- chaselengths_1 %>%
  mutate(ID == "1")

chaselengthall <- rbind(chaselengths_1, chaselengths_2, chaselengths_3, chaselengths_4, chaselengths_5,
                     chaselengths_6, chaselengths_7, chaselengths_8, chaselengths_9, chaselengths_10,
                     chaselengths_11)

chaselength_summary <- chaselengthall

### whip up these histograms:
bdf <- read.csv(here("alldatasetscombined.csv"))

bdf_velocityfilter <- bdf %>%
  group_by(ID, track) %>%
  filter(velocity_mmpersec < 600)

DMSO <- bdf_velocityfilter %>%
  group_by(ID, track) %>%
  filter(Treatment == "DMSO")

ggplot(DMSO, aes(x=velocity_mmpersec, fill=Age))+
  geom_histogram(binwidth=1, alpha=0.6)+
  theme_classic()+
  aes(y = after_stat(count)/sum(after_stat(count))) + 
  scale_y_continuous(labels = scales::percent)

EE <- bdf_velocityfilter %>%
  group_by(ID, track) %>%
  filter(Treatment == "EE")

ggplot(EE, aes(x=velocity_mmpersec, fill=Age))+
  geom_histogram(binwidth=1, alpha=0.6)+
  theme_classic()+
  aes(y = after_stat(count)/sum(after_stat(count))) + 
  scale_y_continuous(labels = scales::percent)