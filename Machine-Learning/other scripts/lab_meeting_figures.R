#Just a little script for formatting figures for lab meeting

#importing data from behavior extraction

#step one: organize all the chase3 figures so they look neat:

chase3 <- bdf %>%
  filter(BehaviorID == 3)

#we would like a nice velocity plot:

ggplot(chase3, aes(x=second, y=velocity_mmpersec))+
  geom_path(linewidth=0.5)+
  xlab("Time (Seconds)")+
  ylab("Velocity (mm/sec)")+
  theme_classic()+
  theme(text = element_text(size = 15))+
  facet_wrap(~track, ncol=1)

#tail angle
ggplot(chase3, aes(x=second, y=theta))+
  geom_path(linewidth=0.5)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle (deg)")+
  theme_classic()+
  theme(text = element_text(size = 15))+
  facet_wrap(~track, ncol=1)

#nearest neighbor
ggplot(chase3, aes(x=second, y=truenn, fill=truenndist))+
  geom_tile()+
  xlab("Time (Seconds)")+
  ylab("Nearest Neighbor (mm)")+
  theme_classic()+
  theme(text = element_text(size = 15))

#we need just a random sample of some swimming to demonstrate tail angle:
#random selection: 36654

tailbeat_sample <- featuredf %>%
  group_by(track) %>%
  filter(frame %in% c(36529:36779))

ggplot(tailbeat_sample, aes(x=second, y=theta))+
  geom_path(linewidth=0.5)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle (deg)")+
  theme_classic()+
  theme(text = element_text(size = 15))+
  facet_wrap(~track, ncol=1)

##nearest neighbor to go with lateral display

lateraldisplay1 <- bdf %>%
  filter(BehaviorID == 1)

ggplot(lateraldisplay1, aes(x=second, y=theta))+
  geom_path(linewidth=0.5)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle (deg)")+
  theme_classic()+
  theme(text = element_text(size = 15))+
  facet_wrap(~track, ncol=1)

ggplot(lateraldisplay1, aes(x=second, y=velocity_mmpersec))+
  geom_path(linewidth=0.5)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle (deg)")+
  theme_classic()+
  theme(text = element_text(size = 15))+
  facet_wrap(~track, ncol=1)

ggplot(lateraldisplay1, aes(x=second, y=truenn, fill=truenndist))+
  geom_tile()+
  xlab("Time (Seconds)")+
  ylab("Nearest Neighbor (mm)")+
  theme_classic()+
  theme(text = element_text(size = 15))

############################################################################

#scott wants subject-directed chases
#this will likely be useful for feature-extraction in the end

#using bdf, begin with identifying the chaser and the chasee:
tmp <- bdf %>%
  group_by(BehaviorID, track)%>%
  slice_max(order_by=velocity_mmpersec, n=1)
#step two, of the 5, now select the top two:
chasemaximums <- tmp %>%
  group_by(BehaviorID) %>%
  slice_max(order_by=velocity_mmpersec, n=2)
rm(tmp)

#the first of these gives us a very good estimation of our chaser
#which we then bind back to bdf, I think
chaser <- chasemaximums %>%
  group_by(BehaviorID) %>%
  slice(1) %>%
  mutate(chaser = track) #this chunk isolates the chaser
chaser <- chaser[c("frame", "chaser")] #take only what we need to merge
bdf <- merge(x = bdf, y = chaser, by = c("frame"), all.x=T) #merge chaser id back to bdf
bdf <- arrange(bdf, track, frame) #arrange by frame
bdf <- bdf %>%
  group_by(track, BehaviorID) %>%
  fill(chaser, .direction=c("downup")) #fill out chaser column
#great. now we figure out how to make this damn figure
#we need all the times the animals aren't chasing...
#merging bdf back into featuredf...
bdf$chaseornochase <- "chase"
bdf_to_merge <- bdf[c("frame", "track", "BehaviorID", "chaser", "chaseornochase")]
featuredf <- merge(x = featuredf, y = bdf_to_merge, by = c("frame", "track"), all.x=T)
featuredf <- arrange(featuredf, track, frame)
featuredf$chaseornochase[is.na(featuredf$chaseornochase)] <- "nochase"

nochase <- featuredf %>%
  group_by(track) %>%
  filter(chaseornochase == "nochase")
nochase$chaser <- nochase$track

chase <- featuredf %>%
  group_by(track) %>%
  filter(chaseornochase == "chase")

nochase <- nochase[c("frame", "track", "chaser")]
chase <- chase[c("frame", "track", "chaser")]
chaserfull <- rbind(chase, nochase)
chaserfull <- arrange(chaserfull, track, frame)
#this is so stupid, but we have the effing column so:
featuredf$chaser <- chaserfull$chaser

#now below for the actual figure:
featuredf$velocity_mmpersec_log10 <- log10(featuredf$velocity_mmpersec)

featuredf_summary_chaser <- featuredf %>%
  filter(chaseornochase == "chase") %>%
  group_by(chaser) %>%
  summarize(velocity_mmpersec_log10 = median(velocity_mmpersec_log10)) %>%
  mutate(chaseornochase = "chase")
featuredf_summary_nochase <- featuredf %>%
  filter(chaseornochase == "nochase") %>%
  group_by(chaser) %>%
  summarize(velocity_mmpersec_log10 = median(velocity_mmpersec_log10, na.rm=T)) %>%
  mutate(chaseornochase = "nochase")
featuredf_summary <- rbind(featuredf_summary_chaser, featuredf_summary_nochase)
rm(featuredf_summary_chaser, featuredf_summary_nochase, nochase, chaserfull, chase, bdf_to_merge)

ggplot(featuredf, aes(x=as.factor(chaser), y=velocity_mmpersec_log10, fill=as.factor(chaseornochase)))+
  geom_violin(alpha=0.2, position=position_dodge())+
  geom_crossbar(featuredf_summary, mapping=aes(x=as.factor(chaser), 
                                                        ymin=velocity_mmpersec_log10, 
                                                        ymax=velocity_mmpersec_log10, 
                                                        color=as.factor(chaseornochase)), 
                position=position_dodge(0.85), width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15), 
        legend.position=c(.9, .15), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))+
  ylab("Log10 Velocity (mm/sec)")+
  xlab("Fish Identity")+
  labs(fill = "Is a Chase \n Happening?", color = "Is a Chase \n Happening?")

#okay and I want to do the same thing, but with theta instead


featuredf_summary_chaser <- featuredf %>%
  filter(chaseornochase == "chase") %>%
  group_by(chaser) %>%
  summarize(theta = median(theta)) %>%
  mutate(chaseornochase = "chase")
featuredf_summary_nochase <- featuredf %>%
  filter(chaseornochase == "nochase") %>%
  group_by(chaser) %>%
  summarize(theta = median(theta, na.rm=T)) %>%
  mutate(chaseornochase = "nochase")
featuredf_summary <- rbind(featuredf_summary_chaser, featuredf_summary_nochase)
rm(featuredf_summary_chaser, featuredf_summary_nochase)

ggplot(featuredf, aes(x=as.factor(chaser), y=theta, fill=as.factor(chaseornochase)))+
  geom_violin(alpha=0.2, position=position_dodge())+
  geom_crossbar(featuredf_summary, mapping=aes(x=as.factor(chaser), 
                                               ymin=theta, 
                                               ymax=theta, 
                                               color=as.factor(chaseornochase)), 
                position=position_dodge(0.85), width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15), 
        legend.position=c(.9, .15), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))+
  ylab("Tail Angle (degrees)")+
  xlab("Fish Identity")+
  labs(fill = "Is a Chase \n Happening?", color = "Is a Chase \n Happening?")

