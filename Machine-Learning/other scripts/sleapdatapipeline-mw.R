#SLEAP Data Basic visualizations + calculations

#library import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#this script is designed for use with dataframes formatted by sleaporganization.r
#this dataframe has information about behaviors annotated by boris data

#import dataset:
trialdata <- read.csv(here("WT031623_26dpf_readyforanalysis.csv"))

#############################################################################################################

#Script To-dos

#Easy stuff
  #add a heat map, both for all animals and for each individual's performance
  #make it so that for the 5 individual figures (position, velocity, acceleration, heatmaps) get saved in one big file together for easy visualization

#Nearest Neighbor
  #calculate Nearest Neighbor among all individuals
    #Molly has a script that does this, but I don't like it and want to work on it
  #filter for moments when animals are close (<2 body lengths, possibly <1 bodylength) away from each other

#"Spike Filter"
  #at Juan's suggestion, the other filter we can use is when velocity and acceleration are high
  #I want to do a couple of statistics on chases before I go in this direction
  #as that will help us identify what the 'average chase' looks like

#Incorporate group cohesion
  #this is an old script I wrote that I'd like to use a bit more of

#Tail beat calc
  #this is mostly just checking over the calculation I wrote

#############################################################################################################
###Figures for each individual fish

##individual position across entire assay:
individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")

triallist <- split(trialdata, trialdata$track)

plotposition <- function(n) {
  ggplot(data=n, aes(x=x, y=y, color=as.factor(nodes)))+
    geom_point(alpha=0.4, size=0.5)+
    geom_path(alpha=0.4, linewidth=0.3)+
    #title(paste('Individual', n$track, sep="_"))+
    xlab("X Direction (mm)")+
    ylab("Y Direction (mm)")+
    labs(color="Skeletal Nodes")+
    theme_classic()+
    coord_cartesian(x=c(0, 1250), y=c(0, 1000))+
    theme(text = element_text(size = 15))  
}

#save to designated folder
positionslist <- lapply(triallist, plotposition)
names(positionslist) <- figuretitle
lapply(names(positionslist), 
       function(x) ggsave(filename=paste(x,"_positions",".png",sep=""),
                          plot=positionslist[[x]], width = 13, height = 7, units = "in",
                          path=here("HotFigures")))
#Ideally I would like to stitch all five figures together into a panel as well

##Velocity Over Time

#velocity calculation, using the head node as our centroid
head_df <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>% #this assumes the video is filmed at 25fps
  mutate(second = frame/25) %>% #stop here if no boris data
  mutate(Behavior = na_if(Behavior, "S")) %>% #these filter out extra behaviors noted in boris
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#individual velocity over time:
individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")
headlist <- split(head_df, head_df$track)

plotheadvelocity <- function(n) {
  ggplot(n, aes(y=velocity_mmpersec, x=second))+
    geom_point(size=0.5)+
    geom_path(linewidth=0.4)+
    xlab("Time (Seconds")+
    ylab("Velocity (mm/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

#save to designated folder
headvelocitylist <- lapply(headlist, plotheadvelocity)
names(headvelocitylist) <- figuretitle
lapply(names(headvelocitylist), 
       function(x) ggsave(filename=paste(x,"_velocity",".png",sep=""),
                          plot=headvelocitylist[[x]], width = 13, height = 7, units = "in",
                          path=here("HotFigures")))

#and the version where all figures displayed in the same chart
facetvelocity <- ggplot(head_df, aes(y=velocity_mmpersec, x=second))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.4)+
  xlab("Time (Seconds")+
  ylab("Velocity (mm/sec)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)

ggsave(filename="facetvelocity.png" ,plot=facetvelocity, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

#color coded for behavior, individual plots:
plotheadvelocity_colorcoded <- function(n) {
  ggplot(n, aes(y=velocity_mmpersec, x=second))+
    geom_point(mapping=aes(color=Behavior, size=Behavior))+
    geom_path(alpha=0.5, linewidth=0.4)+
    scale_size_manual(values= c(5,5,5,5))+
    xlab("Time (Seconds)")+
    ylab("Velocity (mm/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

headvelocitylist_colorcoded <- lapply(headlist, plotheadvelocity_colorcoded)
names(headvelocitylist_colorcoded) <- figuretitle
lapply(names(headvelocitylist_colorcoded), 
       function(x) ggsave(filename=paste(x,"_velocity_colorcoded",".png",sep=""),
                          plot=headvelocitylist_colorcoded[[x]], width = 13, 
                          height = 7, units = "in",
                          path=here("HotFigures")))

#color coded for behavior, facet
facetvelocity_colors <-   
  ggplot(head_df, aes(y=velocity_mmpersec, x=second))+
  geom_point(mapping=aes(color=Behavior, size=Behavior))+
  geom_path(alpha=0.5, linewidth=0.4)+
  scale_size_manual(values= c(2,2,2,2))+
  xlab("Time (Seconds)")+
  ylab("Velocity (mm/sec)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)

ggsave(filename="facetvelocity_colors.png" ,plot=facetvelocity_colors, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

##Acceleration

#calculate from velocity data
head_df <- head_df %>%
  group_by(track) %>%
  mutate(acceleration = velocity_mmpersec - lag(velocity_mmpersec))

individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")
headlist <- split(head_df, head_df$track)

plotheadacceleration_colorcoded <- function(n) {
  ggplot(n, aes(y=acceleration, x=second))+
    geom_point(mapping=aes(color=Behavior, size=Behavior))+ #color coded for behavior
    geom_path(alpha=0.5, linewidth=0.4)+
    scale_size_manual(values= c(5,5,5,5))+
    xlab("Time (Seconds)")+
    ylab("Acceleration (mm/sec/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

headaccelerationlist_colorcoded <- lapply(headlist, plotheadacceleration_colorcoded)
headaccelerationlist_colorcoded

#save individual figures to designated folder
names(headaccelerationlist_colorcoded) <- figuretitle
lapply(names(headaccelerationlist_colorcoded), 
       function(x) ggsave(filename=paste(x,"_acceleration_colorcoded",".png",sep=""),
                          plot=headaccelerationlist_colorcoded[[x]], width = 13, 
                          height = 7, units = "in",
                          path=here("HotFigures")))

#############################################################################################################
### Average Performance

##average velocity

#without outliers
averagevelocity_nooutliers <- ggplot(head_df, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Fish Identity")+
  ylab("Median Velocity")+
  theme_classic()+
  coord_cartesian(y=c(0, 500))+
  theme(text=element_text(size=15))
averagevelocity_nooutliers
#save to designated folder:
ggsave(filename="averagevelocity_nooutliers.png", plot=averagevelocity_nooutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

#average velocity with outliers
averagevelocity_withoutliers <- ggplot(head_df, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_boxplot(outlier.size=0.5)+
  xlab("Fish Identity")+
  ylab("Median Velocity")+
  theme_classic()+
  theme(text=element_text(size=15))
averagevelocity_withoutliers
#save to designated folder:
ggsave(filename="averagevelocity_withoutliers.png", plot=averagevelocity_withoutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

## average acceleration

#acceleration with outliers
averageacceleration <- ggplot(head_df, aes(x=track, y=acceleration, group=as.factor(track)))+
  geom_boxplot(outlier.size=0.5)+
  xlab("Fish Identity")+
  ylab("Median Acceleration (mm/sec/sec)")+
  theme_classic()+
  #coord_cartesian(y=c(0, 500))+
  theme(text=element_text(size=15))
averageacceleration
ggsave(filename="averageacceleration.png", plot=averageacceleration, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="averageacceleration_nooutliers.png", plot=averageacceleration_nooutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

#no outliers
averageacceleration_nooutliers <- ggplot(head_df, aes(x=track, y=acceleration, group=as.factor(track)))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Fish Identity")+
  ylab("Median Acceleration (mm/sec/sec)")+
  theme_classic()+
  coord_cartesian(y=c(-10, 10))+
  theme(text=element_text(size=15))
averageacceleration_nooutliers

##############################################################################

#to do: using bdf, plot metrics in one figure for each associated behavior
#methinks faceted velocity, acceleration, theta, nearest neighbor, and vector components
#and then we can do the nick-style stitch them all together into one figure
#just for reference


