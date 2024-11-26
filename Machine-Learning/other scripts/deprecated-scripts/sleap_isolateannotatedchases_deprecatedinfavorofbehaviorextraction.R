#Characterize a chase

#library import
library(rhdf5)
library(here)
library(tidyverse)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#data import, if necessary:
trialdata <- read.csv(here("WT031623_26dpf_readyforanalysis.csv"))

#filter for chase moments, 2 seconds before and after:

#filter for head:
head_df <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

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

#chases are actually happening really close together to each other, so we'll do +/- 1 second (for now)
leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}
annotatedchasedf_chases <- annotatedchasedf[leadlag(annotatedchasedf$Behavior == "C", 25, 25),]
#and then from here, we can fill out the rest

roll <- function(x, n) {
  idx <- which(!is.na(x))
  for (i in idx) x[pmax(i - n, 0):pmin(i + n, length(x))] <- x[i]
  return(x)
}
annotatedchasedf_chases$ChaseID <- roll(annotatedchasedf_chases$ChaseID, 25)

#and mutate a time column that resets for each bout:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track,ChaseID) %>%
  mutate(timereset = seq_along(cumsum(ChaseID)))

write.csv(annotatedchasedf_chases, file="isolatedchases.csv")

#Right, so. Now we can actually do some visualizations!

annotatedchasedf_chases <- read.csv(here("isolatedchases.csv"))