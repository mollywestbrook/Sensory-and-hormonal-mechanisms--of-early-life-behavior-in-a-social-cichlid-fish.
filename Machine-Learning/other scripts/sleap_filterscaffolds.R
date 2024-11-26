## Behavior Visualizer overall calculations

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)

#bring in the dataset:

featuredf <- read.csv(here("12-21-2023_piz2_DMSO_23dpf_features.csv")) #featuredf
featuredf <- trialdata

#so it turns out velocity is pretty bad at predicting the chasers. So we're gonna try with absolute nearest neighbor instead

### VELOCITY

#Let's start with a visualization of velocity and nearest neighbor for each of the chases...
#Let's just plot the quantiles of each fish per chase...
#there are a few chases that have obscene maximums, so instating a filter for everything <1m/sec

bdf_chasequantiles <- bdf %>%
  group_by(BehaviorID, track) %>%
  reframe('0' = quantile(velocity_mmpersec, probs = c(0), na.rm = T),
            '25' = quantile(velocity_mmpersec, probs = c(0.25), na.rm = T),
            '50' = quantile(velocity_mmpersec, probs = c(0.5), na.rm = T),
            '75' = quantile(velocity_mmpersec, probs = c(0.75), na.rm = T),
            '100' = quantile(velocity_mmpersec, probs = c(1), na.rm = T)) %>%
  gather(key=quantile, value=velocity_mmpersec, '0', '25', '50', '75', '100') %>%
  group_by(BehaviorID, track) %>%
  filter(velocity_mmpersec < 1000)
bdf_chasequantiles$quantile <- as.numeric(bdf_chasequantiles$quantile)

ggplot(bdf_chasequantiles, aes(x=BehaviorID, y=velocity_mmpersec))+
  geom_point(size=1)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

#also do this as a scatter plot grouped by quantile instead of behaviorID as well
ggplot(bdf_chasequantiles, aes(x=velocity_mmpersec, y=quantile))+
  geom_point(size=1)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

#finally, Scott asked for all chases on top of each other...so it will be messy, but I will do it anyway:

#I think still with the <1k filter:
bdf_filter <- bdf %>%
  group_by(BehaviorID, track) %>%
  filter(velocity_mmpersec < 1000)

#trying to filter for ONLY the fish involved in a chase
#this will not be perfect but I will do my best
#I think that using maxvel will work fine for now, rather than NN as NN is still wonky
chasemaximums_tmp <- bdf %>%
  group_by(BehaviorID, track)%>%
  slice_max(order_by=velocity_mmpersec, n=1)
chasemaximums <- chasemaximums_tmp %>%
  group_by(BehaviorID) %>%
  slice_max(order_by=velocity_mmpersec, n=2) %>%
  mutate(fishinchase = T)
rm(chasemaximums_tmp)
#lose the columns in chasemaximums that don't matter...all we need is 
#track behaviorid
chasemaximums <- chasemaximums[,c(2,25,27)]
#then we merge this back in by chaseID, and we should have the fish who are in the chase
bdf_filter_2fish <- merge(bdf_filter, chasemaximums, by=c("BehaviorID", "track"), all.x=T)
bdf_filter_2fish <- arrange(bdf_filter_2fish, track, BehaviorID)
#I think this will work for now!

#so now, we can filter for the two fish with T...let's see if it works:
bdf_fiter_2fish_filter <- bdf_filter_2fish %>%
  group_by(BehaviorID) %>%
  filter(fishinchase == T)

#plot  
ggplot(bdf_fiter_2fish_filter, aes(x=ID, y=velocity_mmpersec, group=BehaviorID, color=BehaviorID))+
  geom_point(size=0.3)+
  geom_line(linewidth=0.2)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

#now he also wants samples that are NOT chases
bdf_filter_2fish <- bdf_filter_2fish %>%
  mutate(chasing = T)
#we must merge this back into annotated behavior so we can filter by whether or not chases
annotatedbehavior_chasemerge <- merge(annotatedbehavior, bdf_filter_2fish, by=c("frame", "track", "nodetrackid", "nodes", "x", "y",
                                                              "second", "distance", "velocity_mmpersec", "acceleration_mmpsps",
                                                              "vx", "vy", "theta", "relativeNN", "relativedist", "truenn",
                                                              "truenndist", "trueheadingangle", "angletofish0",
                                                              "angletofish1", "angletofish2", "angletofish3", "angletofish4"), all.x=T)
annotatedbehavior_chasemerge <- arrange(annotatedbehavior_chasemerge, track, frame)
annotatedbehavior_chasemerge <- annotatedbehavior_chasemerge[-c(24:25)]
names(annotatedbehavior_chasemerge)[c(24:25)] <- c("BehaviorID", "Behavior")
#so now what we need to do is divide this into no chases
#randomly assign points 
#sample around those points using the behavior extraction
#and generate the exact same chaselike figure
bdf_nochases <- annotatedbehavior_chasemerge %>%
  filter(is.na(chasing))

#random draw: 127 numbers
set.seed(498123)
frame <- round(runif(127, min=1, max=53779), digits=0)
nochases <- as.data.frame(frame)
nochases <- nochases %>%
  mutate(Behavior = "F") %>%
  arrange(frame)
#merge back in
bdf_nochases <- merge(bdf_nochases, nochases, by="frame", all.x=T)
#now we take this back over to behaviorextraction etc.
#little bit of name cleanup:
bdf_nochases <- bdf_nochases[-c(24:26)]
names(bdf_nochases)[26] <- "Behavior"

#yay. let us plot!
#do the same filter as before, anything <1000m/s
bdf_nochase_filter <- bdf %>%
  group_by(BehaviorID, track) %>%
  filter(velocity_mmpersec < 1000)
ggplot(bdf_nochase_filter, aes(x=ID, y=velocity_mmpersec, group=BehaviorID, color=BehaviorID))+
  geom_point(size=0.3)+
  geom_line(linewidth=0.2)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

#can I make this a violin plot? Scott sure would love it...
#need to arrange the nochase df so the columns are in the same order as 2fish:
bdf_filter_2fish <- bdf_filter_2fish[c("frame", "track", "nodetrackid", "nodes", "x", "y", "second", "distance", "velocity_mmpersec", "acceleration_mmpsps", "vx", "vy", "theta", "relativeNN", "relativedist", "truenn","truenndist", "trueheadingangle", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4", "fishinchase", "chasing", "Behavior", "BehaviorID", "ID")]
chasecomparison <- rbind(bdf_filter_2fish, bdf_nochase_filter)

#let me also filter for the '2 fish chasing':
bdf_filter_2fish_filter <- bdf_filter_2fish %>%
  filter(fishinchase == T)
chasecomparison_2fish <- rbind(bdf_filter_2fish_filter, bdf_nochase_filter)

#right ok. Now the figgy:

ggplot(chasecomparison_2fish, aes(x=as.factor(track), y=velocity_mmpersec, fill=chasing))+
  geom_violin(position = "dodge", draw_quantiles=T)+
  theme_classic()+
  theme(text=element_text(size=15))+
  coord_cartesian(ylim=c(0, 150))+
  ylab("Velocity (mm/sec)")+
  xlab("Fish Number")


### NEAREST NEIGHBOR

#we also want nearest neighbor stuff. 
#let's just start with the absolute minimum distance per chance. regardless of fish combo:

bdf_nn <- bdf %>%
  group_by(BehaviorID, track) %>%
  slice_min(truenndist) %>%
  select(BehaviorID, track, truenn, truenndist, frame, second, ID)

ggplot(bdf_nn, aes(x=BehaviorID, y=truenndist, color=truenn))+
  geom_point(size=2.5)+
  theme_classic()+
  theme(text=element_text(size=15))
#this will be helpful when compared with a corrected dataset...I went ahead and grabbed the old wt set we've been using

#I think also useful to just do nn quantiles as well, like we did with chases.

bdf_nnquantiles <- bdf %>%
  group_by(BehaviorID, track) %>%
  reframe('0' = quantile(truenndist, probs = c(0), na.rm = T),
          '25' = quantile(truenndist, probs = c(0.25), na.rm = T),
          '50' = quantile(truenndist, probs = c(0.5), na.rm = T),
          '75' = quantile(truenndist, probs = c(0.75), na.rm = T),
          '100' = quantile(truenndist, probs = c(1), na.rm = T)) %>%
  gather(key=quantile, value=truenndist, '0', '25', '50', '75', '100')
bdf_nnquantiles$quantile <- as.numeric(bdf_nnquantiles$quantile)

ggplot(bdf_nnquantiles, aes(x=truenndist, y=quantile))+
  geom_point(size=1)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

#### ACCELERATION

#I am thinking I should probably do this for acceleration

#run the filter again lol
bdf_filter <- bdf %>%
  group_by(BehaviorID, track) %>%
  filter(acceleration_mmpsps < 1000 & acceleration_mmpsps > -1000)

bdf_accelerationquantiles <- bdf_filter %>%
  group_by(BehaviorID, track) %>%
  reframe('0' = quantile(acceleration_mmpsps, probs = c(0), na.rm = T),
          '25' = quantile(acceleration_mmpsps, probs = c(0.25), na.rm = T),
          '50' = quantile(acceleration_mmpsps, probs = c(0.5), na.rm = T),
          '75' = quantile(acceleration_mmpsps, probs = c(0.75), na.rm = T),
          '100' = quantile(acceleration_mmpsps, probs = c(1), na.rm = T)) %>%
  gather(key=quantile, value=acceleration_mmpsps, '0', '25', '50', '75', '100')
bdf_accelerationquantiles$quantile <- as.numeric(bdf_accelerationquantiles$quantile)

ggplot(bdf_accelerationquantiles, aes(x=acceleration_mmpsps, y=quantile))+
  geom_point(size=1)+
  theme_classic()+
  facet_wrap(~track, ncol=1)

####### HEADING ANGLE

#and also some way to summarize heading angle, bc they do seem to be helpful but also I am lost on those...

#maybe I need to plot heading angle with nndist on polar coordinates...


#########################################################################################

#playing around with my older filters...

#origial filter I used (peak velocity to peak velocity)
chasemaximums_tmp <- bdf %>%
  group_by(BehaviorID, track)%>%
  slice_max(order_by=velocity_mmpersec, n=1)
chasemaximums <- chasemaximums_tmp %>%
  group_by(BehaviorID) %>%
  slice_max(order_by=velocity_mmpersec, n=2)
rm(chasemaximums_tmp)

ggplot(chasemaximums, aes(x=ID, y=BehaviorID, group=BehaviorID, color=as.factor(track)))+
  geom_point(size=3)+  
  geom_path(linewidth=0.3, color="black", alpha=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Chase Number")+
  xlab("Time Reset (frame)")

#then for my grand talk, I did acceleration peak to acceleration trough:
chasenn <- bdf %>%
  group_by(BehaviorID)%>%
  slice_min(order_by=relativedist, n=1) %>%
  mutate(selecttrack = 1) %>%
  select(BehaviorID, frame, track, selecttrack) %>%
  arrange(BehaviorID, frame)
bdf_labelled <- merge(bdf, chasenn, by=c("BehaviorID", "track"), all.x=T)
bdf_labelled <- bdf_labelled[-c(27)]
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

ggplot(chaselengths, aes(x=ID, y=BehaviorID, group=BehaviorID, color=as.factor(track)))+
  geom_point(size=3)+  
  geom_path(linewidth=0.3, color="black", alpha=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Chase Number")+
  xlab("Time Reset (frame)")

#so, in order to compare these. I think the key thing to do is to plot like 5 chase examples with these points
#highlighted on a velocity/acceleration curve...

#subset of chases...125, 123, 121, 111, 106
bdf_testsubset <- bdf %>%
  filter(BehaviorID %in% c(106, 111, 121, 123, 125))

chasemaximums_tmp <- bdf_testsubset %>%
  group_by(BehaviorID, track)%>%
  slice_max(order_by=velocity_mmpersec, n=1)
chasemaximums <- chasemaximums_tmp %>%
  group_by(BehaviorID) %>%
  slice_max(order_by=velocity_mmpersec, n=2)
rm(chasemaximums_tmp)

chase106 <- bdf %>%
  filter(BehaviorID == 106)
chase106max <- chasemaximums %>%
  filter(BehaviorID == 106)

ggplot(chase106, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
   geom_line()+
   geom_point(chase106max, mapping=aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
   theme_classic()+
   facet_wrap(~track, ncol=1)

#and the same examples with acceleration peak to trough:
chasenn <- bdf_testsubset %>%
  group_by(BehaviorID)%>%
  slice_min(order_by=relativedist, n=1) %>%
  mutate(selecttrack = 1) %>%
  select(BehaviorID, frame, track, selecttrack) %>%
  arrange(BehaviorID, frame)
bdf_labelled <- merge(bdf, chasenn, by=c("BehaviorID", "track"), all.x=T)
bdf_labelled <- bdf_labelled[-c(27)]
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
chaselengths_max <- chaselengths %>%
  arrange(frame) %>%
  group_by(BehaviorID) %>%
  slice_head()
chaselengths_min <- chaselengths %>%
  arrange(frame) %>%
  group_by(BehaviorID) %>%
  slice_tail()
chaselengths <- rbind(chaselengths_max, chaselengths_min)
chaselengths <- chaselengths %>%
  arrange(frame)

chase125 <- bdf_testsubset %>%
  filter(BehaviorID == 125)
chase125max <- chaselengths %>%
  filter(BehaviorID == 125)

ggplot(chase125, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
  geom_line()+
  geom_point(chase125max, mapping=aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
  theme_classic()+
  facet_wrap(~track, ncol=1)


#I think a filter would look something like this...
# step one: identify points using the below filter -- where nn is in bottom 25% and velocity is in top 75%
# step two: lag to the first acceleration point before that stamp that is in the top 10%? 75% we'll see
# step three: lead to the first point where acceleration of the chaser returns to zero

annotatedbehavior_chasemarked <- annotatedbehavior %>%
  group_by(BehaviorID, track) %>%
  mutate(filtermarker = case_when(truenndist < quantile(truenndist, probs=c(0.25), na.rm=T) & 
           velocity_mmpersec > quantile(velocity_mmpersec, probs=c(0.75), na.rm=T) ~ 1))

#note...when done with the current figure, use a case_when to compare times highlighted by this filter
#vs manual annotation...

#count of times highlighted by chasemarked
#count of times highlighted by behavior
#count of overlapping times

#also, a count of number of behavioral instances that contain at least one instance of this filter...

#hang on, we need to merge bdf back into annotatedbehavior real quick to do this:

annotatedbehavior_merge <- merge(annotatedbehavior, bdf, by=c("frame", "track", "nodetrackid", "nodes", "x", "y",
                                                              "second", "distance", "velocity_mmpersec", "acceleration_mmpsps",
                                                              "vx", "vy", "theta", "relativeNN", "relativedist", "truenn",
                                                              "truenndist", "trueheadingangle", "angletofish0",
                                                              "angletofish1", "angletofish2", "angletofish3", "angletofish4"), all.x=T)
annotatedbehavior_merge <- annotatedbehavior_merge[-c(24:25)]
colnames(annotatedbehavior_merge)[24:25] <- c("Behavior", "BehaviorID")

annotatedbehavior_chasemarked <- annotatedbehavior_merge %>%
  group_by(BehaviorID, track) %>%
  mutate(filtermarker = case_when(truenndist < quantile(truenndist, probs=c(0.25), na.rm=T) & 
                                    velocity_mmpersec > quantile(velocity_mmpersec, probs=c(0.75), na.rm=T) ~ 1)) %>%
  arrange(track, frame)

annotatedbehavior_chasemarked <- annotatedbehavior_chasemarked %>%
  arrange(track, frame)
filtertest <- annotatedbehavior_chasemarked %>%
  group_by(track) %>%
  reframe(count_chasemarked = sum(!is.na(filtermarker)),
          count_behavior = sum(!is.na(Behavior)),
          count_overlap = sum(!is.na(Behavior) & !is.na(filtermarker)))
filtertest_addendum_summary <- annotatedbehavior_chasemarked %>%
  group_by(track, BehaviorID) %>%
  summarize(filtermarkercount = sum(!is.na(filtermarker))) %>%
  group_by(track) %>%
  reframe(missedbyfilter = sum(filtermarkercount == 0),
          hitbyfilter = sum(!filtermarkercount == 0))
falsehits <- annotatedbehavior_chasemarked %>%
  group_by(track, BehaviorID) %>%
  summarize(falsealarm = sum(!is.na(filtermarker))) %>%
  group_by(track) %>% 
  filter(is.na(BehaviorID))
filtertest_summary <- merge(filtertest_addendum_summary, falsehits, by=c("track"))
filtertest_summary <- filtertest_summary[-c(4)]
fwrite(filtertest_summary, file = "filtertest-try2.csv", row.names=FALSE)



