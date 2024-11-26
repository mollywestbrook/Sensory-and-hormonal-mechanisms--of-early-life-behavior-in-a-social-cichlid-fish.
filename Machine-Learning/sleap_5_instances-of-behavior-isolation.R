## Sleap Data Pipeline: Individual Behavioral Analysis

#library Import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#######

#bring in feature df if needed:
featuredf <- read.csv(here("12-21-2023_piz2_DMSO_23dpf_features.csv")) #featuredf
featuredf <- bdf

behaviorofinterest <- "C" #alpha code from ethogram
behaviorwindow <- 60 #generally +/- two seconds

#calculate a unique id for each behavioral instance
tmp <- featuredf %>%
  group_by(track) %>%
  filter(Behavior == behaviorofinterest) %>%
  mutate(BehaviorID = seq_along(Behavior)) %>%
  select(frame, track, BehaviorID)

annotatedbehavior <- merge(x = featuredf, 
                          y = tmp, 
                          by = c("frame", "track"),
                          all.x=T)
annotatedbehavior <- annotatedbehavior %>%
  arrange(track, frame)

leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}

#leadlag select (generally 25)
bdf <- annotatedbehavior[leadlag(annotatedbehavior$Behavior == behaviorofinterest, behaviorwindow, behaviorwindow),]

bdf$Behavior_BehaviorID <- paste(bdf$Behavior, bdf$BehaviorID, sep='-')
bdf$Behavior_BehaviorID[bdf$Behavior_BehaviorID == "NA-NA"] <- NA
bdf$Behavior_BehaviorID[!str_detect(bdf$Behavior_BehaviorID, behaviorofinterest)] <- NA
threshold <- behaviorwindow
bdf <- bdf %>%
  mutate(Group_ID = cumsum(!is.na(Behavior_BehaviorID))) %>%
  group_by(Group_ID) %>%
  mutate(ID = row_number() - 1) %>%
  mutate(Behavior_BehaviorID = ifelse(ID <= threshold, first(Behavior_BehaviorID), NA_character_)) %>%
  ungroup() %>%
  fill(Behavior_BehaviorID, .direction=c("up"))
bdf[c('Behavior', 'BehaviorID')] <- str_split_fixed(bdf$Behavior_BehaviorID, '-', 2)
bdf$BehaviorID <- as.numeric(bdf$BehaviorID)

bdf <- subset(bdf, select=-c(Behavior_BehaviorID, Group_ID))

bdf <- bdf %>%
  arrange(track, frame)

#for some reason the ID column is a little fucked up, so we'll do a re-do real quick seq-along:
bdf <- bdf %>%
  group_by(track, BehaviorID) %>%
  mutate(ID = seq_along(BehaviorID))

#excellent!

#then, also if you want to merge the behavioral windows back into the big df:

annotatedbehavior_merge <- merge(annotatedbehavior, bdf, by=c("frame", "track", "nodetrackid", "nodes", "x", "y",
                                                              "second", "distance", "velocity_mmpersec", "acceleration_mmpsps",
                                                              "vx", "vy", "theta", "relativeNN", "relativedist", "truenn",
                                                              "truenndist", "trueheadingangle", "angletofish0",
                                                              "angletofish1", "angletofish2", "angletofish3", "angletofish4"), all.x=T)
annotatedbehavior_merge <- annotatedbehavior_merge[-c(24:25)]
colnames(annotatedbehavior_merge)[24:25] <- c("Behavior", "BehaviorID")
annotatedbehavior_merge <- annotatedbehavior_merge %>%
  arrange(track, frame)
rm(annotatedbehavior, annotatedbehavior_merge, tmp)

#write_csv(annotatedbehavior_merge, "formatted_dmso-23dpf_12-21-2023_piz2_T4_behavioralwindow.csv")
