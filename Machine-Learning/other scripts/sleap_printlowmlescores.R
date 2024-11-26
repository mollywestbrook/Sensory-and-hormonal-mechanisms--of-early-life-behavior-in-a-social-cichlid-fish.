#Excerpts from Error Analysis
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
library(stringr)

#We just want to generate a list of moments of low MLE for visual analysis
modeldata <- "20231030_1074framemodel_batchtrainingtest.002_08-02-2023_piz3_T6.analysis.h5"
figuretitle <- "08-02-2023_piz3_T6_1074framemodel"

modeldata_bits <- h5ls(here(modeldata))

#we'll need point scores, thus:
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

#then grab timestamps where MLE is <.8 I think to start
#excluding the zeros

#taking data organization from plug-n-play:
modeldata_pointscores <- as.data.frame(h5read(here(modeldata), "/point_scores")) #regrab pointscores
modeldata_nodes <- t(as.data.frame(h5read(here(modeldata), "/node_names"))) #transpose node names frame
modeldata_numtracks <- as.vector(h5read(here(modeldata), "/track_names"))
modeldata_numtracks <- str_remove(modeldata_numtracks, "track_")
modeldata_numtracks_length <- length(as.vector(h5read(here(modeldata), "/track_names"))) #calculate # tracks
modeldata_nodes <- rep(modeldata_nodes, modeldata_numtracks_length) #repeat the number of nodes for each number of tracks
modeldata_nodes <- paste(modeldata_nodes, modeldata_numtracks, sep="-")
modeldata_pointscores <- as.data.frame(modeldata_pointscores) #turn the pointscores back into a data frame
names(modeldata_pointscores) <- c(modeldata_nodes) #change the name of the pointscores to the node names
modeldata_pointscores$frame <- seq_along(modeldata_pointscores[,1]) #create a frame column
modeldata_pointscores_melted <- reshape2::melt(modeldata_pointscores, id=c("frame")) #melt to the frame column so data is in wide format and can be filtered
modeldata_pointscores_melted <- modeldata_pointscores_melted %>%
  separate(variable, c("node", "track"), "-")
#print table of low pointscores and their respective frame

lowpointscores <- modeldata_pointscores_melted %>%
  group_by(frame, track) %>%
  filter(!node %in% c("tail", "caudal")) %>%
  filter(value < 0.7 & !value == 0) %>%
  arrange(track, frame)

#this gives me a conservative list to look through
write.csv(lowpointscores, file="lowpointscores_08-07-2023_piz3_T6.csv")

#just kidding we ACTUALLY need the track scores!
#this will be way easier

#from plugnplay:
modeldata_trackingscores <- as.data.frame(h5read(here(modeldata), "/tracking_scores"))
modeldata_trackingscores <- modeldata_trackingscores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_numtracks <- as.vector(h5read(here(modeldata), "/track_names"))
names(modeldata_trackingscores) <- modeldata_numtracks
modeldata_trackingscores <- mutate(modeldata_trackingscores, ID = row_number())
modeldata_trackingscores <- melt(modeldata_trackingscores, id=c("ID"))

lowtrackscores <- modeldata_trackingscores %>%
  filter(value < 0.7)
write.csv(lowtrackscores, file="lowtrackscores_03-08-2023_group1_R_1010framemodel.csv")
