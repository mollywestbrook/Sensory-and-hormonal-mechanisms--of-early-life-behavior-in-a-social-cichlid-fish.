#sleap autosort followtimes from boris observations

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)

#output of this file is the name of the video and _sleap3output

#bring in file (csv form, ready for analysis)
borisdata <- read.csv(here("10-02-2023_group6_L_25dpf_WT_boris.csv")) #boris

trialdata <- read.csv(here("10-02-2023_group6_L_25dpf_WT.csv")) #sleap (new input)
#trialdata <- tracks_final #sleap (importing from sleap_1)

#some organizing stuff for the boris data
borisdata <- reshape2::melt(borisdata, id.vars = c("Observation.id", "Behavior"), measure.vars = c("Start"))
borisdata <- borisdata[order(borisdata$value),]
borisdata$value <- borisdata$value*30

#round the boris data to the nearest whole number:
borisdata$value <- round(borisdata$value, digits=0)
#and change 'value' to frames to perform a left join:
names(borisdata)[names(borisdata) == "value"] <- "frame"

#now, merge the boris dataframe into the trial dataframe by timestamp:
dataframe_full <- merge(x=trialdata, y=borisdata, by="frame", all=TRUE)
#reorder it
dataframe_full <- dataframe_full %>%
  arrange(nodetrackid, frame)
#drop unneeded columns
dataframe_full <- dataframe_full[-c(8,10)]

#and now we write our csv:
fwrite(dataframe_full, file = "10-02-2023_group6_L_25dpf_WT_sleap3output.csv", row.names=FALSE)

#########################################################################################
#OLD ATTEMPTS

# the old omr method does not realy wor here, so commenting this out for archives
# #remember, we need to set the start time and end time to bind properly:
# borisdata <- rbind(c("20230331_WT031623_26dpf_25fps", "NA", "Start", "1"), borisdata)
# borisdata <- rbind(borisdata, c("20230331_WT031623_26dpf_25fps", "NA", "Start", "44983"))
# borisdata <- borisdata %>%
#   mutate(difference = as.numeric(value) - lag(as.numeric(value)))
# borisdata$difference <- lead(borisdata$difference)
# visualobservations <- unlist(lapply(seq_along(borisdata$difference),
#                                     function(i) rep((i - 1)%% 2, each = borisdata$difference[i])))
# visualobservations <- unlist(lapply(seq_along(borisdata$difference), 
#                                     function(i) rep(borisdata$Behavior, each = borisdata$difference[i])))