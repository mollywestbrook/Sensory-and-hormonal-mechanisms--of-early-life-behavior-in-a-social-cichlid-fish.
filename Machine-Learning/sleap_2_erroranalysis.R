## Plug n Play Error Analysis: Single Model Work

#libraries
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(reticulate)
library(gridExtra)
library(grobblR)
library(ggplotify)

#output of this is a pdf summary of error statistics. 

#file import
#change the file name in quotes to your current .h5 file
modeldata <- "inference_conditions_t8.v001.003_HiQualCamTestVid2.analysis.h5"
figuretitle <- "HiQualCamTestVid2_test"

modeldata_bits <- h5ls(here(modeldata))

tmp <- as.data.frame(h5read(here(modeldata), "/tracks"))

#overall model histograms

#instance scores
modeldata_instancescores <- as.data.frame(h5read(here(modeldata), "/instance_scores"))
modeldata_instancescores <- modeldata_instancescores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_instancescores <- mutate(modeldata_instancescores, ID = row_number())
modeldata_instancescores <- melt(modeldata_instancescores, id=c("ID"))
instancescores <- ggplot(modeldata_instancescores, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Instance Scores For", figuretitle, sep=" "), x="Instance Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
instancescores

#point scores
modeldata_pointscores <- as.data.frame(h5read(here(modeldata), "/point_scores"))
modeldata_pointscores <- modeldata_pointscores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_pointscores <- mutate(modeldata_pointscores, ID = row_number())
modeldata_pointscores <- melt(modeldata_pointscores, id=c("ID"))
pointscores <- ggplot(modeldata_pointscores, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Point Scores For", figuretitle, sep=" "), x="Point Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
pointscores

#tracking scores
modeldata_trackingscores <- as.data.frame(h5read(here(modeldata), "/tracking_scores"))
modeldata_trackingscores <- modeldata_trackingscores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_trackingscores <- mutate(modeldata_trackingscores, ID = row_number())
modeldata_trackingscores <- melt(modeldata_trackingscores, id=c("ID"))
trackingscores <- ggplot(modeldata_trackingscores, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Tracking Scores For", figuretitle, sep=" "), x="Tracking Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
trackingscores

lowtrackscores <- modeldata_trackingscores %>%
  filter(value < 0.7)
write.csv(lowtrackscores, file=paste(figuretitle, "lowtrackscores", ".csv", sep=""))

####

#node-by-node point scores
#organize the data first
modeldata_pointscores <- as.data.frame(h5read(here(modeldata), "/point_scores")) #regrab pointscores
modeldata_nodes <- t(as.data.frame(h5read(here(modeldata), "/node_names"))) #transpose node names frame
modeldata_numtracks <- length(as.vector(h5read(here(modeldata), "/track_names"))) #calculate # tracks
modeldata_nodes <- rep(modeldata_nodes, modeldata_numtracks) #repeat the number of nodes for each number of tracks
modeldata_pointscores <- as.data.frame(modeldata_pointscores) #turn the pointscores back into a data frame
names(modeldata_pointscores) <- c(modeldata_nodes) #change the name of the pointscores to the node names
modeldata_pointscores$frame <- seq_along(modeldata_pointscores[,1]) #create a frame column
modeldata_pointscores_melted <- reshape2::melt(modeldata_pointscores, id=c("frame")) #melt to the frame column so data is in wide format and can be filtered

#nose
modeldata_pointscores_noses <- modeldata_pointscores_melted %>%
  filter(variable == "nose")
nosescores <- ggplot(modeldata_pointscores_noses, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Nose Scores For", figuretitle, sep=" "), x="Nose Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
nosescores

#eyeL
modeldata_pointscores_eyeL <- modeldata_pointscores_melted %>%
  filter(variable == "eyeL")
eyeLscores <- ggplot(modeldata_pointscores_eyeL, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("eyeL Scores For", figuretitle, sep=" "), x="eyeL Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
eyeLscores


#eyeR
modeldata_pointscores_eyeR <- modeldata_pointscores_melted %>%
  filter(variable == "eyeR")
eyeRscores <- ggplot(modeldata_pointscores_eyeR, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("eyeR Scores For", figuretitle, sep=" "), x="eyeR Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
eyeRscores

#head
modeldata_pointscores_head <- modeldata_pointscores_melted %>%
  filter(variable == "head")
headscores <- ggplot(modeldata_pointscores_head, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Head Scores For", figuretitle, sep=" "), x="Head Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
headscores

#spine1
modeldata_pointscores_spine1 <- modeldata_pointscores_melted %>%
  filter(variable == "spine1")
spine1scores <- ggplot(modeldata_pointscores_spine1, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Spine1 Scores For", figuretitle, sep=" "), x="Spine1 Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
spine1scores

#spine2
modeldata_pointscores_spine2 <- modeldata_pointscores_melted %>%
  filter(variable == "spine2")
spine2scores <- ggplot(modeldata_pointscores_spine2, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Spine2 Scores For", figuretitle, sep=" "), x="Spine2 Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
spine2scores

#caudal

modeldata_pointscores_caudal <- modeldata_pointscores_melted %>%
  filter(variable == "caudal")
caudalscores <- ggplot(modeldata_pointscores_caudal, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Caudal Scores For", figuretitle, sep=" "), x="Caudal Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
caudalscores

#tail

modeldata_pointscores_tail <- modeldata_pointscores_melted %>%
  filter(variable == "tail")
tailscores <- ggplot(modeldata_pointscores_tail, aes(x=value, y=stat(count)/sum(stat(count))))+
  geom_histogram(binwidth=0.01)+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(xlim=c(0, 1.2), clip="off")+
  labs(title=paste("Tail Scores For", figuretitle, sep=" "), x="Tail Scores", y="Frequency")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=15), 
        plot.title=element_text(size=12)) 
tailscores

####

#MLE Score summary table

#Hi scores
instancescoresum <- sum(modeldata_instancescores$value > 0.9, na.rm=T)/length(modeldata_instancescores$value)*100
pointscoresum <- sum(modeldata_pointscores_melted$value > 0.9, na.rm=T)/length(modeldata_pointscores_melted$value)*100
trackingscoresum <- sum(modeldata_trackingscores$value > 0.9, na.rm=T)/length(modeldata_trackingscores$value)*100
nosesum <- sum(modeldata_pointscores_noses$value > 0.9, na.rm=T)/length(modeldata_pointscores_noses$value)*100
eyeLsum <- sum(modeldata_pointscores_eyeL$value > 0.9, na.rm=T)/length(modeldata_pointscores_eyeL$value)*100
eyeRsum <- sum(modeldata_pointscores_eyeR$value > 0.9, na.rm=T)/length(modeldata_pointscores_eyeR$value)*100
headsum <- sum(modeldata_pointscores_head$value > 0.9, na.rm=T)/length(modeldata_pointscores_head$value)*100
spine1sum <- sum(modeldata_pointscores_spine1$value > 0.9, na.rm=T)/length(modeldata_pointscores_spine1$value)*100
spine2sum <- sum(modeldata_pointscores_spine2$value > 0.9, na.rm=T)/length(modeldata_pointscores_spine2$value)*100
caudalsum <- sum(modeldata_pointscores_caudal$value > 0.9, na.rm=T)/length(modeldata_pointscores_caudal$value)*100
tailsum <- sum(modeldata_pointscores_tail$value > 0.9, na.rm=T)/length(modeldata_pointscores_tail$value)*100

names <- c("instancescore", "pointscore", "trackingscore", "nose", "head", "eyeL", "eyeR", "spine1", "spine2", "caudal", "tail")
hivalues <- c(instancescoresum, pointscoresum, trackingscoresum, nosesum, headsum, eyeLsum, eyeRsum, spine1sum, spine2sum, caudalsum, tailsum)

hiscoresummary <- data.frame(names, hivalues)

#lo scores
instancescoresum <- sum(modeldata_instancescores$value < 0.5, na.rm=T)/length(modeldata_instancescores$value)*100
pointscoresum <- sum(modeldata_pointscores_melted$value < 0.5, na.rm=T)/length(modeldata_pointscores_melted$value)*100
trackingscoresum <- sum(modeldata_trackingscores$value < 0.5, na.rm=T)/length(modeldata_trackingscores$value)*100
nosesum <- sum(modeldata_pointscores_noses$value < 0.5, na.rm=T)/length(modeldata_pointscores_noses$value)*100
eyeLsum <- sum(modeldata_pointscores_eyeL$value < 0.5, na.rm=T)/length(modeldata_pointscores_eyeL$value)*100
eyeRsum <- sum(modeldata_pointscores_eyeR$value < 0.5, na.rm=T)/length(modeldata_pointscores_eyeR$value)*100
headsum <- sum(modeldata_pointscores_head$value < 0.5, na.rm=T)/length(modeldata_pointscores_head$value)*100
spine1sum <- sum(modeldata_pointscores_spine1$value < 0.5, na.rm=T)/length(modeldata_pointscores_spine1$value)*100
spine2sum <- sum(modeldata_pointscores_spine2$value < 0.5, na.rm=T)/length(modeldata_pointscores_spine2$value)*100
caudalsum <- sum(modeldata_pointscores_caudal$value < 0.5, na.rm=T)/length(modeldata_pointscores_caudal$value)*100
tailsum <- sum(modeldata_pointscores_tail$value < 0.5, na.rm=T)/length(modeldata_pointscores_tail$value)*100

names <- c("instancescore", "pointscore", "trackingscore", "nose", "head", "eyeL", "eyeR", "spine1", "spine2", "caudal", "tail")
lovalues <- c(instancescoresum, pointscoresum, trackingscoresum, nosesum, headsum, eyeLsum, eyeRsum, spine1sum, spine2sum, caudalsum, tailsum)

loscoresummary <- data.frame(names, lovalues)

####

##NA count

#data organization
#note: this is all taken from lines 58-104 for a detailed breakdown of what each of these lines do
modeldata_tracks <- as.data.frame(h5read(here(modeldata), "/tracks")) #grab the tracks
modeldata_track_names <- as.data.frame(h5read(here(modeldata), "/track_names"))
names(modeldata_track_names) <- "tracks"
modeldata_track_names$tracks <- as.numeric(str_replace(modeldata_track_names$tracks, "track_", ""))
modeldata_node_names <- as.data.frame(h5read(here(modeldata), "/node_names"))
names(modeldata_node_names) <- "nodes"
track_names <- rep(modeldata_node_names$nodes, count(modeldata_track_names)*2)
track_names <- reshape2::melt(track_names)
names(track_names) <- "nodes"
tracks <- rep(modeldata_track_names, count(modeldata_node_names)*2)
tracks <- reshape2::melt(tracks)
tracks$value <- sort(tracks$value)
names(tracks) <- c("tracks, label")
track_labels <- as.data.frame(cbind(track_names$nodes, tracks$tracks))
names(track_labels) <- c("nodes", "track")
coordinate_mini <- c(rep("x", count(modeldata_node_names)), rep("y", count(modeldata_node_names)))
coordinate <- rep(coordinate_mini, count(modeldata_track_names))
track_labels <- cbind(track_labels, coordinate)
rm(track_names, tracks, modeldata_node_names, modeldata_track_names) #cleanup
rep_track_labels <- as.data.frame(lapply(track_labels, rep, each=nrow(modeldata_tracks)))
modeldata_tracks_tmp <- reshape2::melt(modeldata_tracks, value.name="data")
tracks_final <- cbind(rep_track_labels, modeldata_tracks_tmp)
tracks_final$track <- as.numeric(tracks_final$track)
tracks_final <- tracks_final %>%
  arrange(track, nodes)
tracks_final$nodetrackid <- cumsum(!duplicated(tracks_final[c(1,2)]))
tracks_final <- tracks_final[-c(4)]
tracks_final <- tracks_final %>%
  group_by(nodetrackid, coordinate) %>%
  mutate(frame = row_number())
tracks_final <- reshape2::dcast(tracks_final, nodetrackid + nodes + track + frame  ~ coordinate, value.var = "data")
rm(modeldata_tracks, modeldata_tracks_tmp, rep_track_labels, track_labels) #cleanup

nodeorder <- c("nose", "eyeL", "eyeR", "head", "spine1", "spine2", "caudal", "tail")
tracks_nacount <- tracks_final %>%
  group_by(track, nodes) %>%
  summarize(nacount = sum(is.na(x))) %>%
  arrange(factor(nodes, levels = nodeorder)) %>%
  arrange(track)
#then we split into 5 mini tables for output:
table0 <- tracks_nacount %>%
  filter(track == 0)
table1 <- tracks_nacount %>%
  filter(track == 1)
table2 <- tracks_nacount %>%
  filter(track == 2)
table3 <- tracks_nacount %>%
  filter(track == 3)
table4 <- tracks_nacount %>%
  filter(track == 4)

table0grob <- tableGrob(table0)
table1grob <- tableGrob(table1)
table2grob <- tableGrob(table2)
table3grob <- tableGrob(table3)
table4grob <- tableGrob(table4)

hiscoresummarygrob <- tableGrob(hiscoresummary)
loscoresummarygrob <- tableGrob(loscoresummary)

title1 <- "Overall Model MLE Scores, compared to Training Dataset"
title2 <- "Node MLE Scores, compared to Training Dataset"
title3 <- "NA Count for each Fish, for each Node"
title4 <- "MLE Score Summary"

page1 = grob_layout(
  grob_row(height=20, grob_col(title1, border=T)),
  grob_row(grob_col(instancescores)),
  grob_row(grob_col(pointscores)),
  grob_row(grob_col(trackingscores)),
  width=200,
  height=200
)
page2 = grob_layout(
  grob_row(height=20, grob_col(title2, border=T)),
  grob_row(grob_col(nosescores)),
  grob_row(grob_col(eyeLscores)),
  grob_row(grob_col(eyeRscores)),
  width=200,
  height=200
)
page3 = grob_layout(
  grob_row(grob_col(headscores)),
  grob_row(grob_col(spine1scores)),
  grob_row(grob_col(spine2scores)),
  width=200,
  height=200
)
page4 = grob_layout(
  grob_row(grob_col(caudalscores)),
  grob_row(grob_col(tailscores)),
  width=200,
  height=200
)

page5 = grob_layout(
  grob_row(height=20, grob_col(title3, border=T)),
  grob_row(grob_col(table0grob), grob_col(table1grob), grob_col(table2grob)),
  grob_row(grob_col(table3grob), grob_col(table4grob)),
  width=200,
  height=200
)

page6 = grob_layout(
  grob_row(height=20, grob_col(title4, border=T)),
  grob_row(grob_col(hiscoresummarygrob), grob_col(loscoresummarygrob)),
  width=200,
  height=200
)

grob_to_pdf(
  page1, page2, page3, page4, page5, page6,
  file_name = paste(figuretitle, "summary", ".pdf", sep=""),
  meta_data_title = "Test PDF"
)
