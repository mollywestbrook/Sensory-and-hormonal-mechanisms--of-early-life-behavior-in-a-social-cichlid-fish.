#Sleap Error Analysis

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(reticulate)

#python notes for python analysis:
#note: python installation is done via
#reticulate::repl_python()
#py_install("lasio",pip=TRUE)
#conda_list()

#activating sleap stuff:
#conda activate sleap
#sleap-label for the gui (should not need this)

filename <- "sleapmodel1_labels.v001.007_20230209_test2_5larva_15dpf_25fps_195frames.analysis.h5"
trialdata_bits <- h5ls(here(filename))

#for each of these files there are three score metrics: instance, point and track
#let's start by generating a histogram of the scores of each of these...
#theoretically they should be at 1

trialdata_instance_scores <- as.data.frame(h5read(here(filename), "/instance_scores"))
trialdata_instance_scores <- trialdata_instance_scores %>% mutate_all(~ifelse(is.nan(.), NA, .))

trialdata_instance_scores <- as.matrix(trialdata_instance_scores)
hist(trialdata_instance_scores, breaks=50, xlab="Instance Scores", ylab="Frequency", 
     main="Instance Scores for 214 frame model")

trialdata_point_scores <- as.data.frame(h5read(here(filename), "/point_scores"))
trialdata_point_scores <- trialdata_point_scores %>% mutate_all(~ifelse(is.nan(.), NA, .))

trialdata_point_scores <- as.matrix(trialdata_point_scores)
hist(trialdata_point_scores, breaks=50, xlab="Point Scores", ylab="Frequency", 
     main="Point Scores for 195 frame model")

trialdata_tracking_scores <- as.data.frame(h5read(here(filename), "/tracking_scores"))
trialdata_tracking_scores <- trialdata_tracking_scores %>% mutate_all(~ifelse(is.nan(.), NA, .))

trialdata_tracking_scores <- as.matrix(trialdata_tracking_scores)
hist(trialdata_tracking_scores, breaks=50, xlab="Point Scores", ylab="Frequency", 
     main="Tracking Scores for 214 frame model")
