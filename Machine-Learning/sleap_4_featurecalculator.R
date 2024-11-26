#Feature Extraction

#library Import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(beepr)
library(swaRmverse)

#output of this is name of video and _sleap4output

#import dataset:
trialdata <- read.csv(here("05-29-2024_piz2_T4_27dpf_DMSO_features.csv"))

featuredf <- trialdata

#set csv name output:
name = "05-29-2024_piz2_T4_27dpf_DMSO_sleap4output.csv"

#set FPS (important only for polarization)
FPS = "0.033S" #for 30 fps
#FPS = "0.04S" #for 25fps

step2time = 0.033 #for 30 fps
#step2time = 0.04 #for 25fps

#or rename from previous sheet:
#trialdata <- dataframe_full

########velocity, acceleration, vector x, vector y, cleanup ####################
featuredf <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(second = frame/30) %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>% #frame-by-frame instantaneous vel calculation
  mutate(velocity_mmpersec = distance*30) %>% #instantaneous vel conversion to mmpersec
  mutate(acceleration_mmpsps = velocity_mmpersec - lag(velocity_mmpersec)) %>% #instantaneous acceleration calculation
  mutate(vx = ((x - lag(x))*30)) %>% #velocity vector x component
  mutate(vy = ((y - lag(y))*30)) %>% #velocity vector y component
  mutate(Behavior = na_if(Behavior, "S")) %>% #set behaviors we annotated for sleap but don't need for analysis to NA
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#idk what was up with repeted rows, but if that happens again run this code:
# trialdata_tmp <- trialdata %>%
#   select(frame, nodetrackid, nodes, track, x, y, second, Behavior) %>%
#   group_by(nodetrackid, frame) %>%
#   slice(1)
# trialdata <- trialdata_tmp
# trialdata <- data.frame(trialdata)

#note: in one single instance there were randomly duplicated rows in 5 random rows. Look out for this in the future...

####### tail angle #####################################################

tailbeats <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "head", "spine1", "caudal"))
tailbeats <- tailbeats[-c(2)]

tailbeats <- unite(tailbeats, col='x-y', c('x', 'y'), sep='-') #group x and y together for manipulation
tailbeats <- spread(tailbeats, key=nodes, value="x-y") #spread out so that nodes are in their own column
tailbeats <- tailbeats %>% #put it in order
  arrange(track,frame)
#string separations
tailbeats <- separate(tailbeats, col=caudal, into=c('caudal.x', 'caudal.y'), sep='-')
tailbeats <- separate(tailbeats, col=head, into=c('head.x', 'head.y'), sep='-')
tailbeats <- separate(tailbeats, col=nose, into=c('nose.x', 'nose.y'), sep='-')
tailbeats <- separate(tailbeats, col=spine1, into=c('spine1.x', 'spine1.y'), sep='-')
#and for some reason there was NA nonsense and str problems, so:
tailbeats[tailbeats == "NA"] <- NA
tailbeats[] <- lapply(tailbeats, as.numeric)

calculatetailangle <- function(ny, hy, nx, hx, cy, sy, cx, sx) {
  slope1 = (ny-hy)/(nx-hx)
  slope2 = (cy-sy)/(cx-sx)
  angle1rad = atan((slope1-slope2)/(1+slope1*slope2))
  angle1 = rad2deg(angle1rad)
  angle1
}

tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(theta = calculatetailangle(nose.y, head.y, nose.x, head.x, caudal.y, spine1.y, caudal.x, spine1.x))

featuredf$theta <- tailbeats$theta

######### Nearest Neighbor ############################################

#note: this will be more useful if relative nearest neighbor is numerical

#to calculate nearest neighbor, we'll need the head dataframe once more:
head_df <- trialdata %>%
  filter(nodes == "head")
#all we really need is frame, track and xy position coordinates, thus:
nearestneighbordf <- head_df[c(1,4,5,6)]
#this does some data manipulation to organize the dataframe:
nearestneighbordf$xy <- paste(nearestneighbordf$x, nearestneighbordf$y, sep="-")
nearestneighbordf <- nearestneighbordf[-c(3,4)]
nearestneighbordf <- nearestneighbordf %>% pivot_wider(names_from = track, values_from = xy)
nearestneighbordf <- data.frame(nearestneighbordf)
names(nearestneighbordf) <- c("frame", "track0xy", "track1xy", "track2xy", "track3xy", "track4xy")
nearestneighbordf[c('track0x', 'track0y')] <- str_split_fixed(nearestneighbordf$track0xy, '-', 2)
nearestneighbordf[c('track1x', 'track1y')] <- str_split_fixed(nearestneighbordf$track1xy, '-', 2)
nearestneighbordf[c('track2x', 'track2y')] <- str_split_fixed(nearestneighbordf$track2xy, '-', 2)
nearestneighbordf[c('track3x', 'track3y')] <- str_split_fixed(nearestneighbordf$track3xy, '-', 2)
nearestneighbordf[c('track4x', 'track4y')] <- str_split_fixed(nearestneighbordf$track4xy, '-', 2)
nearestneighbordf <- nearestneighbordf[-c(2:6)]
names(nearestneighbordf) <- c("frame", "X0", "Y0", "X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4")
nearestneighbordf[] <- lapply(nearestneighbordf, as.numeric)

#now our distance function:
distance_function <- function(x, y, n, p) {
  sqrt((n-x)^2+(p-y)^2)
}

#grab the frame vector, as we'll need it for each individual matrix:
frame <- nearestneighbordf$frame

#calculate our distance matrix and melt it to long form
F0_F1 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X1, nearestneighbordf$Y1)
F0_F2 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X2, nearestneighbordf$Y2)
F0_F3 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X3, nearestneighbordf$Y3)
F0_F4 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X4, nearestneighbordf$Y4)
F1_F2 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X2, nearestneighbordf$Y2)
F1_F3 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X3, nearestneighbordf$Y3)
F1_F4 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X4, nearestneighbordf$Y4)
F2_F3 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X3, nearestneighbordf$Y3)
F2_F4 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X4, nearestneighbordf$Y4)
F3_F4 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X4, nearestneighbordf$Y4)
nnearestneighbor_fulldist<- data.frame(frame, F0_F1, F0_F2, F0_F3, F0_F4, F1_F2, F1_F3, F1_F4, F2_F3, F2_F4, F3_F4)
names(nnearestneighbor_fulldist) <- c("frame", "F0_F1", "F0_F2", "F0_F3", "F0_F4", "F1_F2", "F1_F3", "F1_F4", "F2_F3", "F2_F4", "F3_F4")
rm(F0_F1, F0_F2, F0_F3, F0_F4, F1_F2, F1_F3, F1_F4, F2_F3, F2_F4, F3_F4)
nnearestneighbor_fulldist <- reshape2::melt(nnearestneighbor_fulldist, id="frame")

#right here we want to grab group distance from Tang 2020 paper
#mean pairwise distances across whole group -- they summarize it with median
#but we just want it for the entire dataframe

groupdistance <- nnearestneighbor_fulldist %>%
  group_by(frame) %>%
  summarize(groupdistance = mean(value, na.rm=T)) %>%
  mutate(groupdistance = na_if(groupdistance, NaN))

featuredf <- merge(featuredf, groupdistance, by=c("frame"), all.x=T)
featuredf <- featuredf %>%
  arrange(track, frame)

#this calculates true nearest neighbor
nnearestneighbor_fulldist <- nnearestneighbor_fulldist %>%
  group_by(frame) %>%
  slice(which.min(value))
names(nnearestneighbor_fulldist) <- c("frame", "truenn", "truenndist")

#in order to get actual nearest neighbor, we'll have to be a little funky with our distance matrix.
#We'll calculate a distance matrix for each animal
#then select the animal that is the closest to that particular individual at any particular frame:

frame <- data.frame(frame) #for frame length preservation

#fish one:
F1 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish0_NN <- data.frame(frame, F1, F2, F3, F4)
rm(F1, F2, F3, F4)
Fish0_NN <- reshape2::melt(Fish0_NN, id="frame")
Fish0_NN <- Fish0_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish0_NN) <- c("frame", "relativeNN", "relativedist")
Fish0_NN <- merge(frame, Fish0_NN, by=c("frame"), all.x=TRUE) #preserve frame length
Fish0_NN$track <- 0

#fish two
F0 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X0, nearestneighbordf$Y0)
F2 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish1_NN <- data.frame(frame, F0, F2, F3, F4)
rm(F0, F2, F3, F4)
Fish1_NN <- reshape2::melt(Fish1_NN, id="frame")
Fish1_NN <- Fish1_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish1_NN) <- c("frame", "relativeNN", "relativedist")
Fish1_NN <- merge(frame, Fish1_NN, by=c("frame"), all.x=TRUE)
Fish1_NN$track <- 1

#fish three
F0 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X1, nearestneighbordf$Y1)
F3 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish2_NN <- data.frame(frame, F0, F1, F3, F4)
rm(F0, F1, F3, F4)
Fish2_NN <- reshape2::melt(Fish2_NN, id="frame")
Fish2_NN <- Fish2_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish2_NN) <- c("frame", "relativeNN", "relativedist")
Fish2_NN <- merge(frame, Fish2_NN, by=c("frame"), all.x=TRUE)
Fish2_NN$track <- 2

#fish four
F0 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X2, nearestneighbordf$Y2)
F4 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish3_NN <- data.frame(frame, F0, F1, F2, F4)
rm(F0, F1, F2, F4)
Fish3_NN <- reshape2::melt(Fish3_NN, id="frame")
Fish3_NN <- Fish3_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish3_NN) <- c("frame", "relativeNN", "relativedist")
Fish3_NN <- merge(frame, Fish3_NN, by=c("frame"), all.x=TRUE)
Fish3_NN$track <- 3

#fish five
F0 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X3, nearestneighbordf$Y3)
Fish4_NN <- data.frame(frame, F0, F1, F2, F3)
rm(F0, F1, F2, F3)
Fish4_NN <- reshape2::melt(Fish4_NN, id="frame")
Fish4_NN <- Fish4_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish4_NN) <- c("frame", "relativeNN", "relativedist")
Fish4_NN <- merge(frame, Fish4_NN, by=c("frame"), all.x=TRUE)
Fish4_NN$track <- 4


#final df:
NNDF <- rbind(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN)
NNDF <- merge(NNDF, nnearestneighbor_fulldist, by=c("frame"), all.x=T)
NNDF <- NNDF %>%
  arrange(track, frame)

#we also want average nearest neighbor, therefore:
avgnn <- NNDF %>%
  group_by(frame) %>%
  summarize(avgnn = mean(relativedist, na.rm=T)) %>%
  mutate(avgnn = na_if(avgnn, NaN))

featuredf <- merge(featuredf, avgnn, by=c("frame"), all.x=T)

featuredf <- merge(featuredf, NNDF, by=c("frame", "track"), all.x=T) #this line has not been checked thoroughly; first glance is working
featuredf <- featuredf %>%
  arrange(track, frame)

#cleanup:
rm(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN, nnearestneighbor_fulldist, frame)
rm(avgnn, groupdistance, NNDF)

########## true heading angle (relative to (0,0) in rad) ##############################

featuredf$trueheadingangle <- atan2(featuredf$y, featuredf$x) 
featuredf$trueheadingangle <- rad2deg(featuredf$trueheadingangle)
featuredf$trueheadingangle <- 90 - featuredf$trueheadingangle #make it relative to yaxis

#then we want to calculate polarity (again from that Tang 2020 paper)

#in Tungstrom 2013 they define it as the absolute value of the mean heading angle 
#we should not have to worry about negative values here but I'll include it anyway
#instead of suffering anymore, we're just gonna try and calculate polarization with
#swarmverse

#need a succinct dataset for this one:
polarizationdf <- set_data_format(raw_x = featuredf$x,
                                  raw_y = featuredf$y,
                                  raw_t = featuredf$frame,
                                  raw_id = featuredf$track,
                                  period = "0.04S")
polarizationdf <- add_velocities(polarizationdf)

groupmetrics <- group_metrics(polarizationdf[[1]], geo=F, step2time = step2time)
#so this spits out a df, all we need is the polarization column
#to bind back into the featuredf, therefore:

pol <- groupmetrics$pol
polarization <- data.frame(pol, seq(pol))
names(polarization) <- c("pol", "frame")

#and merge to the featuredf:
featuredf <- merge(featuredf, polarization, by=c("frame"), all.x=T)
featuredf <- featuredf %>%
  arrange(track, frame)

########## relative heading angle in deg #############################################

# #old calculation:
# calculatebearingangle <- function(X1, Y1, X2, Y2) {
#   bearing_rad = atan2(X2 - X1, Y1 - Y2)
#   bearing = rad2deg(bearing_rad)
#   bearing
# }

#redoing this calculation: 
headingslopes <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("head", "spine1"))
headingslopes <- headingslopes[-c(2)]

#note: grab the frame vector, we'll ened it later:
frame <- nearestneighbordf$frame

headingslopes <- unite(headingslopes, col='x-y', c('x', 'y'), sep='_') #group x and y together for manipulation
headingslopes <- spread(headingslopes, key=nodes, value="x-y") #spread out so that nodes are in their own column
headingslopes <- headingslopes %>% #put it in order
  arrange(track,frame)
#separate head and spine1 xy points:
headingslopes <- separate(headingslopes, col=head, into=c('head.x', 'head.y'), sep='_')
headingslopes <- separate(headingslopes, col=spine1, into=c('spine1.x', 'spine1.y'), sep='_')
headingslopes[headingslopes == "NA"] <- NA
headingslopes[] <- lapply(headingslopes, as.numeric)
#calculate the head-spine x vector and head-spine y vector:
headingslopes$xvector <- headingslopes$spine1.x - headingslopes$head.x
headingslopes$yvector <- headingslopes$spine1.y - headingslopes$head.y
#now, we can calculate relative angles between each of these vectors...
#drop the head and spine xy points:
headingslopes <- headingslopes[c(1:4, 9, 10)]
#and we do the same exact manipulation as the nearest neighbor df:
headingslopes$xy <- paste(headingslopes$xvector, headingslopes$yvector, sep="_")
headingslopes <- headingslopes[-c(5,6)]

headingslopes <- spread(headingslopes, key=track, value="xy")


names(headingslopes) <- c("frame", "second", "Behavior", "track0xy", "track1xy", "track2xy", "track3xy", "track4xy")
headingslopes[c('track0x', 'track0y')] <- str_split_fixed(headingslopes$track0xy, '_', 2)
headingslopes[c('track1x', 'track1y')] <- str_split_fixed(headingslopes$track1xy, '_', 2)
headingslopes[c('track2x', 'track2y')] <- str_split_fixed(headingslopes$track2xy, '_', 2)
headingslopes[c('track3x', 'track3y')] <- str_split_fixed(headingslopes$track3xy, '_', 2)
headingslopes[c('track4x', 'track4y')] <- str_split_fixed(headingslopes$track4xy, '_', 2)
headingslopes <- headingslopes[-c(4:8)]
names(headingslopes) <- c("frame", "second", "Behavior", "X0", "Y0", "X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4")
headingslopes[] <- lapply(headingslopes, as.numeric)
#excellent. Now for the actual calculation:

calculatebearingangle <- function(X1, Y1, X2, Y2) {
  bearing_1 = atan2(X1, Y1)
  bearing_2 = atan2(X2, Y2)
  bearingangle_rad = bearing_1 - bearing_2
  bearingangle_deg_premodulo = rad2deg(bearingangle_rad)
  bearingangle_deg = bearingangle_deg_premodulo + 360 %% 360
  bearingangle_deg
}

#the calculation for each fish pair with repeats
F0_F0 <- calculatebearingangle(headingslopes$X0, headingslopes$Y0, headingslopes$X0, headingslopes$Y0)
F0_F1 <- calculatebearingangle(headingslopes$X0, headingslopes$Y0, headingslopes$X1, headingslopes$Y1)
F0_F2 <- calculatebearingangle(headingslopes$X0, headingslopes$Y0, headingslopes$X2, headingslopes$Y2)
F0_F3 <- calculatebearingangle(headingslopes$X0, headingslopes$Y0, headingslopes$X3, headingslopes$Y3)
F0_F4 <- calculatebearingangle(headingslopes$X0, headingslopes$Y0, headingslopes$X4, headingslopes$Y4)
F1_F0 <- calculatebearingangle(headingslopes$X1, headingslopes$Y1, headingslopes$X0, headingslopes$Y0)
F1_F1 <- calculatebearingangle(headingslopes$X1, headingslopes$Y1, headingslopes$X1, headingslopes$Y1)
F1_F2 <- calculatebearingangle(headingslopes$X1, headingslopes$Y1, headingslopes$X2, headingslopes$Y2)
F1_F3 <- calculatebearingangle(headingslopes$X1, headingslopes$Y1, headingslopes$X3, headingslopes$Y3)
F1_F4 <- calculatebearingangle(headingslopes$X1, headingslopes$Y1, headingslopes$X4, headingslopes$Y4)
F2_F0 <- calculatebearingangle(headingslopes$X2, headingslopes$Y2, headingslopes$X0, headingslopes$Y0)
F2_F1 <- calculatebearingangle(headingslopes$X2, headingslopes$Y2, headingslopes$X1, headingslopes$Y1)
F2_F2 <- calculatebearingangle(headingslopes$X2, headingslopes$Y2, headingslopes$X2, headingslopes$Y2)
F2_F3 <- calculatebearingangle(headingslopes$X2, headingslopes$Y2, headingslopes$X3, headingslopes$Y3)
F2_F4 <- calculatebearingangle(headingslopes$X2, headingslopes$Y2, headingslopes$X4, headingslopes$Y4)
F3_F0 <- calculatebearingangle(headingslopes$X3, headingslopes$Y3, headingslopes$X0, headingslopes$Y0)
F3_F1 <- calculatebearingangle(headingslopes$X3, headingslopes$Y3, headingslopes$X1, headingslopes$Y1)
F3_F2 <- calculatebearingangle(headingslopes$X3, headingslopes$Y3, headingslopes$X2, headingslopes$Y2)
F3_F3 <- calculatebearingangle(headingslopes$X3, headingslopes$Y3, headingslopes$X3, headingslopes$Y3)
F3_F4 <- calculatebearingangle(headingslopes$X3, headingslopes$Y3, headingslopes$X4, headingslopes$Y4)
F4_F0 <- calculatebearingangle(headingslopes$X4, headingslopes$Y4, headingslopes$X0, headingslopes$Y0)
F4_F1 <- calculatebearingangle(headingslopes$X4, headingslopes$Y4, headingslopes$X1, headingslopes$Y1)
F4_F2 <- calculatebearingangle(headingslopes$X4, headingslopes$Y4, headingslopes$X2, headingslopes$Y2)
F4_F3 <- calculatebearingangle(headingslopes$X4, headingslopes$Y4, headingslopes$X3, headingslopes$Y3)
F4_F4 <- calculatebearingangle(headingslopes$X4, headingslopes$Y4, headingslopes$X4, headingslopes$Y4)

#build the df relative to each fish to line up with our feature extraction
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
bearingangle_featuredf <- rbind(f0df, f1df, f2df, f3df, f4df) #bind dataframe together
bearingangle_featuredf <- bearingangle_featuredf[-c(1)]
rm(F0_F0, F0_F1, F0_F2, F0_F3, F0_F4,
   F1_F0, F1_F1, F1_F2, F1_F3, F1_F4,
   F2_F0, F2_F1, F2_F2, F2_F3, F2_F4,
   F3_F0, F3_F1, F3_F2, F3_F3, F3_F4,
   F4_F0, F4_F1, F4_F2, F4_F3, F4_F4,
   f0df, f1df, f2df, f3df, f4df)

featuredf <- cbind(featuredf, bearingangle_featuredf)

featuredf <- featuredf %>%
  group_by(frame, track) %>%
  arrange(track, frame)

#one more arrangement for stuff: make things exactly equal to zero na:
featuredf[, 20:24][featuredf[, 20:24] == 0] <- NA


######## body length ############################################

# #this may actually be useful for a potential frame alerter
# bodylength <- trialdata %>%
#   group_by(track) %>%
#   filter(nodes %in% c("nose", "caudal")) %>%
#   select(frame, nodes, track, x, y) %>%
#   unite(col='x-y', c('x', 'y'), sep='-') %>%
#   spread(key=nodes, value="x-y") %>%
#   separate(col=nose, into=c('nose.x', 'nose.y'), sep='-') %>%
#   separate(col=caudal, into=c('caudal.x', 'caudal.y'), sep='-') %>%
#   arrange(track, frame) %>%
#   mutate_at(c("caudal.x", "caudal.y", "nose.x", "nose.y"), as.numeric) %>%
#   mutate(bodylength = distance_function(nose.x, nose.y, caudal.x, caudal.y))

#when you're ready to move to analysis, run this line to clean up your environment:
#rm(bearingangle_featuredf, bearingangledf, nearestneighbordf, NNDF, tailbeats, frame)

write_csv(featuredf, name)

beep(sound =11)

