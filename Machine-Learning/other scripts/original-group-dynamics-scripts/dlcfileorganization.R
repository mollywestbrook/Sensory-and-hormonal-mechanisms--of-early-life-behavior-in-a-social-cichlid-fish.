#DLC Organization



library(tidyverse)
library(here)
library(reshape2)
library(data.table)
library(gridExtra)
library(zoo)

#Step 1: Import CSV and modify 

trialdata <- read_csv(here("test3DLC_dlcrnetms5_Aggression TestApr25shuffle1_100000_el.csv"), col_names=FALSE)

#Due to DLC's syntax, the best way of keeping this code flexible is to create our own titles based on the first three rows of information
#note, this will require use of grepl to filter, but we don't have to devote a stupid amount of code to renaming things

colnames(trialdata) <- paste(trialdata[1,], trialdata[2,], trialdata[3,], sep="_") #new column names
trialdata <- trialdata[-c(1:3), ] #drop the first three rows, as they're now our colnames
colnames(trialdata)[1] <- "frame" #rename 1st column frame

#probably use likelihood to get a sense of error on each animal, and then we can just get rid of them

#[will come back to this eventually]

#here I think it's best to just do a quick visual check of each of the animals and see if we're losing identities throughout the run:

ggplot(trialdata, aes(x=fish1_nose_x, y=fish1_nose_y))+
  geom_point(size=1)+
  geom_path()+
  theme_classic()


