## Tracking csv conversion and basic figures for DLC

library(tidyverse)
library(here)
library(schoolmath)
library(reshape2)
library(data.table)
library(gridExtra)
library(nortest)
library(swfscMisc)
library(devtools)
library(ggpmisc)
library(viridis)
library(RColorBrewer)

#for bonsai
trialdata<- read_csv(here("20220425_collectivebehaviortest_15larva_partitionpanicroom_2_partition2_bonsai.csv"), col_names = FALSE)

#for IDtracker
trialdata<- read_csv(here("20220519_15juvenile_6mo_panicroompartition_30min_partition2.csv"), col_names = TRUE)

### data organization for bonsai

names(trialdata) <- c("Fish1", "Fish2", "Fish3", "Fish4", "Fish5")

trialdata$Fish1 <- str_replace(trialdata$Fish1, "[(]|[)]", "")
trialdata$Fish2 <- str_replace(trialdata$Fish2, "[(]|[)]", "")
trialdata$Fish3 <- str_replace(trialdata$Fish3, "[(]|[)]", "")
trialdata$Fish4 <- str_replace(trialdata$Fish4, "[(]|[)]", "")
trialdata$Fish5 <- str_replace(trialdata$Fish5, "[(]|[)]", "")
trialdata$Fish1 <- str_replace(trialdata$Fish1, "[(]|[)]", "")
trialdata$Fish2 <- str_replace(trialdata$Fish2, "[(]|[)]", "")
trialdata$Fish3 <- str_replace(trialdata$Fish3, "[(]|[)]", "")
trialdata$Fish4 <- str_replace(trialdata$Fish4, "[(]|[)]", "")
trialdata$Fish5 <- str_replace(trialdata$Fish5, "[(]|[)]", "")

trialdata[,6:7] <- str_split_fixed(trialdata$Fish1, ",", 2)
trialdata[,8:9] <- str_split_fixed(trialdata$Fish2, ",", 2)
trialdata[,10:11] <- str_split_fixed(trialdata$Fish3, ",", 2)
trialdata[,12:13] <- str_split_fixed(trialdata$Fish4, ",", 2)
trialdata[,14:15] <- str_split_fixed(trialdata$Fish5, ",", 2)

trialdata <- trialdata[,6:15]

names(trialdata) <- c("X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4", "X5", "Y5")

#everything is annoyingly a character, so:

trialdata$X1 <- as.numeric(trialdata$X1)
trialdata$Y1 <- as.numeric(trialdata$Y1)
trialdata$X2 <- as.numeric(trialdata$X2)
trialdata$Y2 <- as.numeric(trialdata$Y2)
trialdata$X3 <- as.numeric(trialdata$X3)
trialdata$Y3 <- as.numeric(trialdata$Y3)
trialdata$X4 <- as.numeric(trialdata$X4)
trialdata$Y4 <- as.numeric(trialdata$Y4)
trialdata$X5 <- as.numeric(trialdata$X5)
trialdata$Y5 <- as.numeric(trialdata$Y5)

#

### Data organization for IDTracker

#all we need to do is lose the probID columns:

trialdata <- trialdata[-c(3,6,9,12,15)]

#

### Convert to mm

#dish conversion to mm:
trialdata$X1 <- trialdata$X1*0.196
trialdata$Y1 <- trialdata$Y1*0.196
trialdata$X2 <- trialdata$X2*0.196
trialdata$Y2 <- trialdata$Y2*0.196
trialdata$X3 <- trialdata$X3*0.196
trialdata$Y3 <- trialdata$Y3*0.196
trialdata$X4 <- trialdata$X4*0.196
trialdata$Y4 <- trialdata$Y4*0.196 
trialdata$X5 <- trialdata$X5*0.196
trialdata$Y5 <- trialdata$Y5*0.196

#tank conversion to mm
trialdata$X1 <- trialdata$X1*0.357
trialdata$Y1 <- trialdata$Y1*0.357
trialdata$X2 <- trialdata$X2*0.357
trialdata$Y2 <- trialdata$Y2*0.357
trialdata$X3 <- trialdata$X3*0.357
trialdata$Y3 <- trialdata$Y3*0.357
trialdata$X4 <- trialdata$X4*0.357
trialdata$Y4 <- trialdata$Y4*0.357 
trialdata$X5 <- trialdata$X5*0.357
trialdata$Y5 <- trialdata$Y5*0.357

#Extra: bonsai partition 
trialdata$X1 <- trialdata$X1*0.397
trialdata$Y1 <- trialdata$Y1*0.397
trialdata$X2 <- trialdata$X2*0.397
trialdata$Y2 <- trialdata$Y2*0.397
trialdata$X3 <- trialdata$X3*0.397
trialdata$Y3 <- trialdata$Y3*0.397
trialdata$X4 <- trialdata$X4*0.397
trialdata$Y4 <- trialdata$Y4*0.397 
trialdata$X5 <- trialdata$X5*0.397
trialdata$Y5 <- trialdata$Y5*0.397

#Slap a frame on there:

trialdata <- trialdata %>%
  mutate(Frame = row_number())

trialdata <- trialdata %>%
  mutate(Second = Frame/25)

################################################################################################################################
# 
# #basic visualization to see if this even worked:
# 
# ggplot(trialdata, aes(x=X1, y=Y1, color=Frame))+
#   geom_point(size=1)+
#   #coord_cartesian(xlim=c(0,800), ylim=c(0,800))+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=X2, y=Y2, color=Frame))+
#   geom_point(size=1)+
#   #coord_cartesian(xlim=c(0,800), ylim=c(0,800))+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=X3, y=Y3, color=Frame))+
#   geom_point(size=1)+
#   #coord_cartesian(xlim=c(0,800), ylim=c(0,800))+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=X4, y=Y4, color=Frame))+
#   geom_point(size=1)+
#   #coord_cartesian(xlim=c(0,800), ylim=c(0,800))+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=X5, y=Y5, color=Frame))+
#   geom_point(size=1)+
#   #coord_cartesian(xlim=c(0,800), ylim=c(0,800))+
#   theme_classic()
# 
# ggplot()+
#   geom_point(trialdata, mapping=aes(X1, Y1), color="red", size=1, alpha=0.4)+
#   geom_point(trialdata, mapping=aes(X2, Y2), color="blue", size=1, alpha=0.4)+
#   geom_point(trialdata, mapping=aes(X3, Y3), color="darkgreen", size=1, alpha=0.4)+
#   geom_point(trialdata, mapping=aes(X4, Y4), color="purple", size=1, alpha=0.4)+
#   geom_point(trialdata, mapping=aes(X5, Y5), color="orange", size=1, alpha=0.4)+
#   theme_classic()
# 
# #divide up into chunks:
# 
# trialdata_5min <- trialdata %>%
#   filter(Frame %in% c(0:7500))
# trialdata_10min <- trialdata %>%
#   filter(Frame %in% c(7501:15000))
# trialdata_15min <- trialdata %>%
#   filter(Frame %in% c(15001:22500))
# trialdata_20min <- trialdata %>%
#   filter(Frame %in% c(225001:30000))
# trialdata_25min <- trialdata %>%
#   filter(Frame %in% c(30001:37500))
# 
# #divide into point sampling:
# 
# pointsample <- trialdata %>%
#   filter(Frame %in% c(1, 5200, 7849, 12190, 17880, 21500, 27300, 30000))
# 
# ggplot()+
#   geom_point(pointsample, mapping=aes(x=X1, y=Y1, color=as.factor(Frame)), size=2)+
#   #geom_line(pointsample, mapping=aes(x=X1, y=Y1))+
#   geom_point(pointsample, mapping=aes(x=X2, y=Y2, color=as.factor(Frame)), size=2)+
#   #geom_line(pointsample, mapping=aes(x=X2, y=Y2))+
#   geom_point(pointsample, mapping=aes(x=X3, y=Y3, color=as.factor(Frame)), size=2)+
#   #geom_line(pointsample, mapping=aes(x=X3, y=Y3))+
#   geom_point(pointsample, mapping=aes(x=X4, y=Y4, color=as.factor(Frame)), size=2)+
#   #geom_line(pointsample, mapping=aes(x=X4, y=Y4))+
#   geom_point(pointsample, mapping=aes(x=X5, y=Y5, color=as.factor(Frame)), size=2)+
#   #geom_line(pointsample, mapping=aes(x=X5, y=Y5))+
#   scale_color_manual(values=c("blue", "brown2", "bisque", "darkgoldenrod1", "cyan", "green", "darkorange", "darkorchid1"))+
#   theme_classic()

##########################################################################################3

#Velocity Calculations

trialdata <- trialdata %>%
  mutate(Fish1Velocity = sqrt((X1 - lag(X1))^2 + (Y1 - lag(Y1))^2),
         Fish2Velocity = sqrt((X2 - lag(X2))^2 + (Y2 - lag(Y2))^2),
         Fish3Velocity = sqrt((X3 - lag(X3))^2 + (Y3 - lag(Y3))^2),
         Fish4Velocity = sqrt((X4 - lag(X4))^2 + (Y4 - lag(Y4))^2),
         Fish5Velocity = sqrt((X5 - lag(X5))^2 + (Y5 - lag(Y5))^2),
         )

#also will be useful just to have a dedicated df for only velocity:

velocitydf <- trialdata[c(13:17)]

#Quick visualization

# ggplot(trialdata, aes(x=Second, y=Fish1Velocity))+
#   geom_line(size=0.4)+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=Second, y=Fish2Velocity))+
#   geom_line(size=0.4)+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=Second, y=Fish3Velocity))+
#   geom_line(size=0.4)+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=Second, y=Fish4Velocity))+
#   geom_line(size=0.4)+
#   theme_classic()
# 
# ggplot(trialdata, aes(x=Second, y=Fish5Velocity))+
#   geom_line(size=0.4)+
#   theme_classic()

#calling outliers those above 20

trialdata$Fish1Velocity <- replace(trialdata$Fish1Velocity, trialdata$Fish1Velocity > 20, NA)
trialdata$Fish2Velocity <- replace(trialdata$Fish2Velocity, trialdata$Fish2Velocity > 20, NA)
trialdata$Fish3Velocity <- replace(trialdata$Fish3Velocity, trialdata$Fish3Velocity > 20, NA)
trialdata$Fish4Velocity <- replace(trialdata$Fish4Velocity, trialdata$Fish4Velocity > 20, NA)
trialdata$Fish5Velocity <- replace(trialdata$Fish5Velocity, trialdata$Fish5Velocity > 20, NA)

#now we have to calculate a few other things:

#quantiles we need to calculate: 0.05, 0.25, 0.75

quants <- trialdata %>%
  summarize(cutoff = quantile(c(Fish1Velocity, Fish2Velocity, Fish3Velocity, Fish4Velocity, Fish5Velocity), probs = 0.05, na.rm=T),
            twentyfifth = quantile(c(Fish1Velocity, Fish2Velocity, Fish3Velocity, Fish4Velocity, Fish5Velocity), probs = 0.25, na.rm=T),
            seventyfifth = quantile(c(Fish1Velocity, Fish2Velocity, Fish3Velocity, Fish4Velocity, Fish5Velocity), probs = 0.75, na.rm=T))

#############################################################################################

#Group Centroid Stuff

#centroid is calculated as average of points (https://gis.stackexchange.com/questions/6025/finding-centroid-of-cluster-of-points-using-r)

trialdata$groupcentroidx <- apply(trialdata[c(1,3,5,7,9)], MARGIN =  1, FUN = mean, na.rm = T)
trialdata$groupcentroidy <- apply(trialdata[c(2,4,6,8,10)], MARGIN =  1, FUN = mean, na.rm = T)

trialdata <- trialdata %>%
  mutate(centroidspeed = sqrt((groupcentroidx - lag(groupcentroidx))^2 + (groupcentroidy - lag(groupcentroidy))^2))
