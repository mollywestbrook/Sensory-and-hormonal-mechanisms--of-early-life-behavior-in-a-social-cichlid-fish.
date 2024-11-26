# Draft for sleap data pipeline

#library import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#bring in data
trialdata <- read.csv(here("formatted_05-29-2023_group5_R.csv"))

#####################################################################################
#assemble head point df for calcs where we only need one point:

#filter for head:
head_df <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

###################################################################333333

#Chase stuff

#seq along for each instance of chase, for grouping later:
tmp <- head_df %>%
  group_by(track) %>%
  filter(Behavior == "C") %>%
  mutate(ChaseID = seq_along(Behavior))
#then, merge the labels back in for analysis
annotatedchasedf <- merge(x = head_df, 
                  y = tmp, 
                  by = c("frame", "nodetrackid", "nodes", "track", "x", "y", "Behavior", 
                         "distance", "velocity_mmpersec", "second"),
                  all.x = TRUE)
annotatedchasedf <- annotatedchasedf %>%
  arrange(track, frame)

#this will select +/- two seconds before and after each chase
#sec increased to +/- 4 seconds, as we don't really know how long chases tend to last on average
leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}
annotatedchasedf_chases <- annotatedchasedf[leadlag(annotatedchasedf$Behavior == "C", 100, 100),]
#and then from here, we can fill out the rest
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  fill(ChaseID, .direction=c("downup"))
#and mutate a time column that resets for each bout:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track,ChaseID) %>%
  mutate(timereset = seq_along(cumsum(ChaseID)))

#really quick, a df of just blanket categories of whether or not a chase is happening
annotatedchasedf_nonchases <- annotatedchasedf[!leadlag(annotatedchasedf$Behavior == "C", 100, 100),]
annotatedchasedf_nonchases <- annotatedchasedf_nonchases %>%
  mutate(chaseorna = "NoChase")
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  mutate(chaseorna = "Chases")
annotatedchasedf <- rbind(annotatedchasedf_nonchases, annotatedchasedf_chases)
#Right, so. Now we can actually do some visualizations!

ggplot(annotatedchasedf_chases, aes(x=timereset, y=velocity_mmpersec))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.3)+
  facet_wrap(~track, ncol=1)

#ok. we're gonna do a boxplot of each chase, and facet by individual. Just to give a sense of speed for each one
ggplot(annotatedchasedf_chases, aes(x=as.factor(ChaseID), y=velocity_mmpersec))+
  geom_boxplot()+
  facet_wrap(~track, ncol=1)

#let's try doing all chases compared to everything that is not a chase, per individual...
ggplot(annotatedchasedf, aes(x=chaseorna, y=velocity_mmpersec))+
  geom_boxplot(outlier.shape=NA)+
  coord_cartesian(ylim=c(0, 1000))+
  facet_wrap(~track, ncol=1)

#calculate some stats for this...
chasestats <- annotatedchasedf_chases %>%
  group_by(track, ChaseID) %>%
  summarize(minimum_vel = min(velocity_mmpersec),
            quantile25_vel = quantile(velocity_mmpersec, probs=c(.25)),
            quantile50_vel = quantile(velocity_mmpersec, probs=c(.50)),
            quantile75_vel = quantile(velocity_mmpersec, probs=c(.75)),
            maximum_vel = max(velocity_mmpersec),
            mean_vel = mean(velocity_mmpersec)
            )
write.csv(chasestats, file="chasestats.csv")

####################################################################################

#Vector field of chases

#let's use a clean head_df so we can calculate the vector components
#and we'll also clean up the behaviors as well
head_df_vectors <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(vx = ((x - lag(x))*25)) %>%
  mutate(vy = ((y - lag(y))*25)) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}

chase_vectors <- head_df_vectors[leadlag(head_df_vectors$Behavior == "C", 10, 10),]

#chase 1:
chase1 <- chase_vectors %>%
  group_by(track) %>%
  slice(1:21)
chase1plot <- ggplot(chase1, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase2 <- chase_vectors %>%
  group_by(track) %>%
  slice(22:42)
chase2plot <-  ggplot(chase2, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase3 <- chase_vectors %>%
  group_by(track) %>%
  slice(43:63) %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>%
  mutate(acceleration = velocity_mmpersec - lag(velocity_mmpersec))

chase3plot <- ggplot(chase3, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=20))
#to go with chase three I'd like to do a facet of acceleration + velocity with our fish of interest:
chase3_velocity <- ggplot(chase3, aes(x=second, y=velocity_mmpersec))+
  geom_point()+
  geom_path()+
  theme_classic()+
  ylab("Velocity (mm/sec)")+
  xlab("Time (Seconds)")+
  theme(text=element_text(size=18))+
  facet_wrap(~track, ncol=1)
ggsave("chase3_velocity.png", plot=chase3_velocity, width=4, height=8, unit="in")
chase3_acceleration <- ggplot(chase3, aes(x=second, y=acceleration))+
  geom_point()+
  geom_path()+
  theme_classic()+
  ylab("Acceleration (mm/sec/sec)")+
  xlab("Time (Seconds)")+
  theme(text=element_text(size=18))+
  facet_wrap(~track, ncol=1)
ggsave("chase3_acceleration.png", plot=chase3_acceleration, width=4, height=8, unit="in")


chase4 <- chase_vectors %>%
  group_by(track) %>%
  slice(64:84)
chase4plot <- ggplot(chase4, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase5 <- chase_vectors %>%
  group_by(track) %>%
  slice(85:105)
chase5plot <- ggplot(chase5, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

ggsave(filename="chase1.png", plot=chase1plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase2.png", plot=chase2plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase3.png", plot=chase3plot, 
       width = 7, height = 5, units = "in", path=here("HotFigures"))
ggsave(filename="chase4.png", plot=chase4plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase5.png", plot=chase5plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

####################################################################################
#Make the nearest neighbor calculation more flexible!

#once again we want to focus on the head df. 
#let's assemble just the essential components of this df, as most of this we don't need

nndf <- head_df[-c(2,3,7:10)]
nndf$xy <- paste(nndf$x, nndf$y, sep='-')
nndf$track <- paste("track", nndf$track, sep='')
nndf <- nndf[-c(3,4)]

#################################################################################

#Tail Beat Dynamics:
#in order to model the tail beat, we'll need these:
#vector between nose and head (this will be our 'base' to measure against)
#and the vector between spine and caudal (this will be the thing that rotates)

#therefore, we'll arrange the dataset with the things we need first.
#this will require some casting.
tailbeats <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "head", "spine1", "caudal"))
tailbeats <- tailbeats[-c(2)]

tailbeats <- unite(tailbeats, col='x-y', c('x', 'y'), sep='-') #group x and y together for manipulation
tailbeats <- spread(tailbeats, key=nodes, value="x-y") #spread out so that nodes are in their own column
tailbeats <- tailbeats %>% #put it in order
  arrange(track,frame)

#spread it back out so columns have their proper labels
tailbeats <- separate(tailbeats, col=caudal, into=c('caudal.x', 'caudal.y'), sep='-')
tailbeats <- separate(tailbeats, col=head, into=c('head.x', 'head.y'), sep='-')
tailbeats <- separate(tailbeats, col=nose, into=c('nose.x', 'nose.y'), sep='-')
tailbeats <- separate(tailbeats, col=spine1, into=c('spine1.x', 'spine1.y'), sep='-')
#and for some reason there was NA nonsense and str problems, so:
tailbeats[tailbeats == "NA"] <- NA
tailbeats$caudal.x <- as.numeric(tailbeats$caudal.x)
tailbeats$caudal.y <- as.numeric(tailbeats$caudal.y)
tailbeats$head.x <- as.numeric(tailbeats$head.x)
tailbeats$head.y <- as.numeric(tailbeats$head.y)
tailbeats$nose.x <- as.numeric(tailbeats$nose.x)
tailbeats$nose.y <- as.numeric(tailbeats$nose.y)
tailbeats$spine1.x <- as.numeric(tailbeats$spine1.x)
tailbeats$spine1.y <- as.numeric(tailbeats$spine1.y)

################################################

calculatetailangle <- function(ny, hy, nx, hx, cy, sy, cx, sx) {
  slope1 = (ny-hy)/(nx-hx)
  slope2 = (cy-sy)/(cx-sx)
  angle1rad = atan((slope1-slope2)/(1+slope1*slope2))
  angle1 = rad2deg(angle1rad)
  angle1
}

#we'll need the vector between nose and head, and the vector between spine 1 and caudal:
tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(theta = calculatetailangle(nose.y, head.y, nose.x, head.x, caudal.y, spine1.y, caudal.x, spine1.x))

#and now the plots.

#just for visualization's sake, I will filter this a little bit heh.

tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

individ <- unique(tailbeats$track)
figuretitle <- paste('Individual', individ, sep="_")
thetalist <- split(tailbeats, tailbeats$track)

thetafacet <-  
  ggplot(tailbeats, aes(y=theta, x=frame))+
  geom_point(mapping=aes(size=Behavior, color=Behavior, group=Behavior))+
  geom_path(linewidth=0.3)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle Relative to Head (Radians)")+
  scale_size_manual(values= c(2,2,2,2,0.5))+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)
thetafacet

ggsave(filename="thetafacet.png",
       plot=thetafacet, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

plottheta <- function(n) {
  ggplot(n, aes(y=thetadeg, x=frame))+
    geom_point(mapping=aes(color=Behavior, size=Behavior, group=Behavior))+
    geom_path(linewidth=0.3)+
    xlab("Time (Seconds)")+
    ylab("Tail Angle Relative to Head (Radians)")+
    scale_size_manual(values= c(5,5,5,5))+
    theme_classic()+
    theme(text = element_text(size = 18))  
}

thetaplotlist <- lapply(thetalist, plottheta)
thetaplotlist

names(thetaplotlist) <- figuretitle
lapply(names(thetaplotlist), 
       function(x) ggsave(filename=paste(x,"_tailangle_colorcoded",".png",sep=""),
                          plot=thetaplotlist[[x]], width = 10, 
                          height = 5, units = "in",
                          path=here("HotFigures")))

#of just the lateral displays:
lateraldisplays <- tailbeats[leadlag(tailbeats$Behavior == "L", 20, 20),]

lateraldisplayfacet <-  
  ggplot(lateraldisplay1, aes(y=theta, x=frame))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.3)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle Relative to Head (Degrees)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)
lateraldisplayfacet

ggsave("lateraldisplayfacet.png", plot=lateraldisplayfacet, width=7, height=5,
       unit="in", path=here("HotFigures"))

#a little bit of filtering to just look at the lateral displays we have, using the same method as chasing
tmp <- head_df %>%
  group_by(track) %>%
  filter(Behavior == "L") %>%
  mutate(LatID = seq_along(Behavior))
#then, merge the labels back in for analysis
annotateddisplaydf <- merge(x = head_df, 
                          y = tmp, 
                          by = c("frame", "nodetrackid", "nodes", "track", "x", "y", "Behavior", 
                                 "distance", "velocity_mmpersec", "second"),
                          all.x = TRUE)
annotateddisplaydf <- annotateddisplaydf %>%
  arrange(track, frame)
annotateddisplaydf <- annotateddisplaydf[-c(8,9)]
annotateddisplaydf$theta <- tailbeats$theta
#
lateraldisplays <- annotateddisplaydf[leadlag(annotateddisplaydf$Behavior == "L", 50, 50),]
#annotation will simply be done by hand as fill is not working for me:
lateraldisplays <- lateraldisplays %>%
  group_by(track) %>%
  mutate(LatID = case_when(frame %in% c(12040:12140) ~ 1,
                          frame %in% c(26029:26129) ~ 2,
                          frame %in% c(30521:30621) ~ 3))
#and mutate a time column that resets for each bout:
lateraldisplays <- lateraldisplays %>%
  group_by(track,LatID) %>%
  mutate(timereset = seq_along(cumsum(LatID)))

#now the actual plot for each display
latlist <- split(lateraldisplays, lateraldisplays$LatID)
plotdisplay <- function(n){
  ggplot(n, aes(x=second, y=theta, group=LatID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
latplotlist <- lapply(latlist, plotdisplay)
latnames <- paste("lateraldisplay", unique(lateraldisplays$LatID), sep='-')
names(latplotlist) <- latnames
lapply(names(latplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""), width=7, height=7, units=c("in"), plot=latplotlist[[x]], path=here("WT031623_26dpf_Figures")))

#We're also gonna do this for chases:
#grabbing code from up above
tmp <- head_df %>%
  group_by(track) %>%
  filter(Behavior == "C") %>%
  mutate(ChaseID = seq_along(Behavior))
#then, merge the labels back in for analysis
annotatedchasedf <- merge(x = head_df, 
                          y = tmp, 
                          by = c("frame", "nodetrackid", "nodes", "track", "x", "y", "Behavior", 
                                 "distance", "velocity_mmpersec", "second"),
                          all.x = TRUE)
annotatedchasedf <- annotatedchasedf %>%
  arrange(track, frame)
#bind in theta:
annotatedchasedf$theta <- tailbeats$theta
#pull out chases
annotatedchasedf_chases <- annotatedchasedf[leadlag(annotatedchasedf$Behavior == "C", 25, 25),]
#and then from here, we need to expand the numbers to their corresponding chases
#this fill works!!!
annotatedchasedf_chases$Behavior_ChaseID <- paste(annotatedchasedf_chases$Behavior, annotatedchasedf_chases$ChaseID, sep='-')
annotatedchasedf_chases$Behavior_ChaseID[annotatedchasedf_chases$Behavior_ChaseID == "NA-NA"] <- NA
threshold <- 25
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  mutate(Group_ID = cumsum(!is.na(Behavior_ChaseID))) %>%
  group_by(Group_ID) %>%
  mutate(ID = row_number() - 1) %>%
  mutate(Behavior_ChaseID = ifelse(ID <= threshold, first(Behavior_ChaseID), NA_character_)) %>%
  ungroup() %>%
  fill(Behavior_ChaseID, .direction=c("up"))
#little bit of df cleanup:
annotatedchasedf_chases <- annotatedchasedf_chases[-c(7, 11, 14, 15)]
annotatedchasedf_chases[c('Behavior', 'ChaseID')] <- str_split_fixed(annotatedchasedf_chases$Behavior_ChaseID, '-', 2)
annotatedchasedf_chases <- annotatedchasedf_chases[-c(11,12)]
annotatedchasedf_chases$ChaseID <- as.numeric(annotatedchasedf_chases$ChaseID)
#and mutate a time column that resets for each bout:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track,ChaseID) %>%
  mutate(timereset = seq_along(cumsum(ChaseID)))

chaselist <- split(annotatedchasedf_chases, annotatedchasedf_chases$ChaseID)
plottheta <- function(n){
  ggplot(n, aes(x=second, y=theta, group=ChaseID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
chaseplotlist <- lapply(chaselist, plottheta)
chasenames <- paste("chase", unique(annotatedchasedf_chases$ChaseID), sep='-')
names(chaseplotlist) <- chasenames
lapply(names(chaseplotlist), 
       function(x) ggsave(filename=paste(x,"theta", ".png",sep=""), width=7, height=7, units=c("in"), plot=chaseplotlist[[x]], path=here("WT031623_26dpf_Figures")))

chaselist <- split(annotatedchasedf_chases, annotatedchasedf_chases$ChaseID)
plottheta <- function(n){
  ggplot(n, aes(x=second, y=theta, group=ChaseID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
chaseplotlist <- lapply(chaselist, plottheta)
chasenames <- paste("chase", unique(annotatedchasedf_chases$ChaseID), sep='-')
names(chaseplotlist) <- chasenames
lapply(names(chaseplotlist), 
       function(x) ggsave(filename=paste(x,"theta", ".png",sep=""), width=7, height=7, units=c("in"), plot=chaseplotlist[[x]], path=here("WT031623_26dpf_Figures")))

#now that we've corrected our grouping, let's also redo the chase velocity spikes as well:
plotvelocity <- function(n){
  ggplot(n, aes(x=second, y=velocity_mmpersec, group=ChaseID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
chaseplotlist <- lapply(chaselist, plotvelocity)
chasenames <- paste("chase", unique(annotatedchasedf_chases$ChaseID), sep='-')
names(chaseplotlist) <- chasenames
lapply(names(chaseplotlist), 
       function(x) ggsave(filename=paste(x,"velocity", ".png",sep=""), plot=chaseplotlist[[x]], path=here("WT031623_26dpf_Figures")))

#idk if 'average theta' is going to tell us anything...I think theta range is a better statistic
thetasummary <- annotatedchasedf_chases %>%
  group_by(track, ChaseID) %>%
  summarize(min = min(theta, na.rm=TRUE), max = max(theta, na.rm=T), thetarange = max - min)
thetaaveragerange <- thetasummary %>%
  group_by(track) %>%
  summarize(thetarange = mean(thetarange))

ggplot(thetasummary, aes(x=track, y=thetarange, group=track))+
  geom_violin()+
  geom_jitter(width=0.1, alpha=0.3, size=0.8)+
  geom_crossbar(thetaaveragerange, mapping=aes(x=track, ymin=thetarange, ymax=thetarange), color="red", width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Range of Tail Angle (Degrees)")+
  xlab("Fish IDentity")


#for lab meeting next week, I want average nearest neighbor during a chase
#I will summarize for absolute minimum

nnsummary <- bdf %>%
  select(BehaviorID, truenn, truenndist) %>%
  group_by(BehaviorID) %>%
  reframe(absoluteminimum = min(truenndist), truenn = first(truenn))

nnsummaryaverage <- nnsummary %>%
  group_by(truenn) %>%
  summarize(absoluteminimum = mean(absoluteminimum))

#now I can graph:

ggplot(nnsummary, aes(x=truenn, y=absoluteminimum))+
  geom_bar(nnsummaryaverage, mapping=aes(x=truenn, y=absoluteminimum, fill=truenn), stat="identity")+
  geom_jitter(width=0.1, alpha=0.5)+
  ylab("Nearest Neighbor Distance (mm)")+
  xlab("Fish Identities of the Nearest Neighbors")+
  theme_classic()+
  theme(text=element_text(size=15))

#well, this makes me ask herm, the body length of the fish?
#so, let's do that. it will take some fancy casting again but w/e
bodylength <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "caudal")) %>%
  select(frame, nodes, track, x, y) %>%
  unite(col='x-y', c('x', 'y'), sep='-') %>%
  spread(key=nodes, value="x-y") %>%
  separate(col=nose, into=c('nose.x', 'nose.y'), sep='-') %>%
  separate(col=caudal, into=c('caudal.x', 'caudal.y'), sep='-') %>%
  arrange(track, frame) %>%
  mutate_at(c("caudal.x", "caudal.y", "nose.x", "nose.y"), as.numeric) %>%
  mutate(bodylength = distance_function(nose.x, nose.y, caudal.x, caudal.y))

#kinda want to see this on a jitter plot...
ggplot(bodylength, aes(x=track, y=bodylength))+
  geom_jitter(width=0.2)

###############################################################################

#Heading angle drafts

#We'll need to use the same strategy as nearest neighbor
#as we'll need relative heading angle of each fish
#code taken from the feature extraction script from nearest neighbor
head_df <- trialdata %>%
  filter(nodes == "head")
bearingangledf <- head_df[c(1,4,5,6)]
bearingangledf$xy <- paste(bearingangledf$x, bearingangledf$y, sep="-")
bearingangledf <- bearingangledf[-c(3,4)]
bearingangledf <- reshape2::dcast(bearingangledf, frame~track)
names(bearingangledf) <- c("frame", "track0xy", "track1xy", "track2xy", "track3xy", "track4xy")
bearingangledf[c('track0x', 'track0y')] <- str_split_fixed(bearingangledf$track0xy, '-', 2)
bearingangledf[c('track1x', 'track1y')] <- str_split_fixed(bearingangledf$track1xy, '-', 2)
bearingangledf[c('track2x', 'track2y')] <- str_split_fixed(bearingangledf$track2xy, '-', 2)
bearingangledf[c('track3x', 'track3y')] <- str_split_fixed(bearingangledf$track3xy, '-', 2)
bearingangledf[c('track4x', 'track4y')] <- str_split_fixed(bearingangledf$track4xy, '-', 2)
bearingangledf <- bearingangledf[-c(2:6)]
names(bearingangledf) <- c("frame", "X0", "Y0", "X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4")
bearingangledf[] <- lapply(bearingangledf, as.numeric)
frame <- bearingangledf$frame

calculatebearingangle <- function(X1, Y1, X2, Y2) {
  bearing_rad = atan2(X2 - X1, Y1 - Y2)
  bearing = rad2deg(bearing_rad)
  bearing
}

#So we'll need to do each fish to each other fish
#WITH repeats -_-
F0_F0 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X0, bearingangledf$Y0)
F0_F1 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X1, bearingangledf$Y1)
F0_F2 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X2, bearingangledf$Y2)
F0_F3 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X3, bearingangledf$Y3)
F0_F4 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X4, bearingangledf$Y4)
F1_F0 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X0, bearingangledf$Y0)
F1_F1 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X1, bearingangledf$Y1)
F1_F2 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X2, bearingangledf$Y2)
F1_F3 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X3, bearingangledf$Y3)
F1_F4 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X4, bearingangledf$Y4)
F2_F0 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X0, bearingangledf$Y0)
F2_F1 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X1, bearingangledf$Y1)
F2_F2 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X2, bearingangledf$Y2)
F2_F3 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X3, bearingangledf$Y3)
F2_F4 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X4, bearingangledf$Y4)
F3_F0 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X0, bearingangledf$Y0)
F3_F1 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X1, bearingangledf$Y1)
F3_F2 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X2, bearingangledf$Y2)
F3_F3 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X3, bearingangledf$Y3)
F3_F4 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X4, bearingangledf$Y4)
F4_F0 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X0, bearingangledf$Y0)
F4_F1 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X1, bearingangledf$Y1)
F4_F2 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X2, bearingangledf$Y2)
F4_F3 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X3, bearingangledf$Y3)
F4_F4 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X4, bearingangledf$Y4)

f0df <- data.frame(frame,F0_F0, F0_F1, F0_F2, F0_F3, F0_F4)
f1df <- data.frame(frame,F1_F0, F1_F1, F1_F2, F1_F3, F1_F4)
f2df <- data.frame(frame,F2_F0, F2_F1, F2_F2, F2_F3, F2_F4)
f3df <- data.frame(frame,F3_F0, F3_F1, F3_F2, F3_F3, F3_F4)
f4df <- data.frame(frame,F4_F0, F4_F1, F4_F2, F4_F3, F4_F4)
names(f0df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f1df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f2df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f3df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f4df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
bearingangle_featuredf <- rbind(f0df, f1df, f2df, f3df, f4df)
rm(F0_F0, F0_F1, F0_F2, F0_F3, F0_F4,
   F1_F0, F1_F1, F1_F2, F1_F3, F1_F4,
   F2_F0, F2_F1, F2_F2, F2_F3, F2_F4,
   F3_F0, F3_F1, F3_F2, F3_F3, F3_F4,
   F4_F0, F4_F1, F4_F2, F4_F3, F4_F4)

#that should do it folks!

#################################################################################

#just a couple of xy plots for meeting today

trialdata <- tracks_final

first200frames <- trialdata %>%
  group_by(track) %>%
  filter(frame %in% c(0:200))

ggplot(first200frames, aes(x=x, y=y, color=as.factor(track)))+
  geom_point(size=1, alpha=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))

