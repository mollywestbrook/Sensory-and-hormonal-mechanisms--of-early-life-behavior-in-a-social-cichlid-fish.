#Actual chase analysis breakdown
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)

write.csv(annotatedchasedf_full, file="chasedffull.csv")
featuredf <- read.csv(here("formatted_dmso-23dpf_12-21-2023_piz2_T4_features.csv"))
annotatedchasedf_chases <- bdf

#this figure is not very good so I'll leave it up here for now
#violin plot visualization for Scott, each chase
ggplot(annotatedchasedf_chases, aes(x=ChaseID, y=velocity_mmpersec, group=ChaseID))+
  geom_violin()+
  theme_classic()+
  scale_y_continuous(trans='log10') +
  facet_wrap(~track, ncol=1)

### Chase Visualization

#step one: capture a large enough time window to visualize each chase
#plot each chase, fiddling around with time windows in #isolatechases
chaselist <- split(bdf, bdf$BehaviorID)
plotchase <- function(n){
  ggplot(n, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
chaseplotlist <- lapply(chaselist, plotchase)
chasenames <- paste("chase", unique(bdf$BehaviorID), sep='-')
names(chaseplotlist) <- chasenames
lapply(names(chaseplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""), plot=chaseplotlist[[x]], path=here("12-21-2023_piz1_T2_Figures")))

plotchase_accel <- function(n){
  ggplot(n, aes(x=second, y=acceleration_mmpsps, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
chaseplotlist <- lapply(chaselist, plotchase_accel)
chasenames <- paste("chase-acceleration", unique(bdf$BehaviorID), sep='-')
names(chaseplotlist) <- chasenames
lapply(names(chaseplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""), plot=chaseplotlist[[x]], path=here("12-21-2023_piz1_T2_Figures")))

#seems like chases are still highly variable so
#try and identify the maximum and it's subsequent maximum?

#this will be a two step filter
#step one: pull out maximum velocity across each fish per chase bout
#five maximums total
chasemaximums_tmp <- annotatedchasedf_chases %>%
  group_by(ChaseID, track)%>%
  slice_max(order_by=velocity_mmpersec, n=1)
#step two, of the 5, now select the top two:
chasemaximums <- chasemaximums_tmp %>%
  group_by(ChaseID) %>%
  slice_max(order_by=velocity_mmpersec, n=2)
rm(chasemaximums_tmp)
#haha! very promising :)

#Thinking now about how I want to do a scatter
#try 1: 
#first, a summary so we can have the average:
chasemaximums_summary <- chasemaximums %>%
  group_by(track) %>%
  summarize(velocity_mmpersec = mean(velocity_mmpersec, na.rm=TRUE))

#next we do a violin plot of the maximum velocity of chases
ggplot(chasemaximums, aes(x=track, y=velocity_mmpersec, group=track))+
  geom_violin()+
  geom_jitter(width=0.2)+
  geom_crossbar(chasemaximums_summary, mapping=aes(x=track, ymin=velocity_mmpersec, ymax=velocity_mmpersec), color="red", width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Velocity During Chase (mm/sec)")+
  xlab("Fish IDentity")

#now what i want is color coded to which fish, with some kind of visualization
#of the two timestamps per chase
ggplot(chasemaximums, aes(x=timereset, y=ChaseID, group=ChaseID, color=as.factor(track)))+
  geom_point(size=3)+  
  geom_path(size=0.3, color="black", alpha=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Chase Number")+
  xlab("Time Reset (frame)")

#and to go with this figure, what tends to be the length of chases?
chaselengths <- chasemaximums %>%
  group_by(ChaseID) %>%
  summarize(chaselength = max(timereset)-min(timereset)) %>%
  mutate(chaselength_sec = chaselength/25) %>%
  mutate(group = 1)

chaseaverage <- chaselengths %>%
  group_by(group) %>%
  summarize(chaselength_sec = mean(chaselength_sec, na.rm=TRUE))

ggplot(chaselengths, aes(y=group, x=chaselength_sec))+
  geom_violin()+
  geom_jitter(width=0.1)+
  geom_crossbar(chaseaverage, mapping=aes(xmin=chaselength_sec,
                                          xmax=chaselength_sec,
                                          y=group),
                color="red")+
  theme_classic()+
  theme(text=element_text(size=15))


#This is all awesome

#########################################################################

#going back to chase velocity, let's contextualize the chases within the rest of the dataframe:

#we already made this graph, but I want to just mark it again
#as this is what we'll be focusing on:
ggplot(chasemaximums, aes(x=track, y=velocity_mmpersec, group=track))+
  geom_violin()+
  geom_jitter(width=0.2)+
  geom_crossbar(chasemaximums_summary, mapping=aes(x=track, ymin=velocity_mmpersec, ymax=velocity_mmpersec), color="red", width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Velocity During Chase (mm/sec)")+
  xlab("Fish IDentity")

#in order to contextualize these
#Ill stick with the arbitrary time windows in annotatedchasedf_chases
#and merge them back in to the main df:

annotatedchasedf_full <- left_join(annotatedchasedf, annotatedchasedf_chases,
                                   by=c("frame","nodetrackid","nodes","track","x","y",
                                        "Behavior","distance","velocity_mmpersec","second"))
#and then some cleanup:
annotatedchasedf_full <- annotatedchasedf_full[-c(11)]
names(annotatedchasedf_full)[c(11,12)] <- c("ChaseID", "timereset")
#beaut

#now, we can take a look at average velocity of each fish
#when the chase window is up and at every other time:

annotatedchasedf_full <- annotatedchasedf_full %>%
  mutate(isachase = case_when(!is.na(ChaseID) ~ "chase",
                              is.na(ChaseID) ~ "nochase"))

#and the visualization (let's start messy and refine):
#jitterplot for each fish by chase/nochase, with violins overlaid bc Scott likes them:

#summary first, to make this stupidly complicated:
annotatedchasefull_summary <- annotatedchasedf_full %>%
  group_by(track, isachase) %>%
  summarize(velocity_mmpersec = mean(velocity_mmpersec, na.rm=TRUE)) %>%
  mutate(velocity_mmpersec_log10 = log10(velocity_mmpersec))

#to plot my averages, we must hand log10 transform the velocity first:
annotatedchasedf_full <- annotatedchasedf_full %>%
  mutate(velocity_mmpersec_log10 = log10(velocity_mmpersec))
  

ggplot(annotatedchasedf_full, aes(x=as.factor(track), y=velocity_mmpersec_log10, fill=as.factor(isachase)))+
  geom_violin(alpha=0.2, draw_quantiles = T)+
  geom_crossbar(annotatedchasefull_summary, mapping=aes(x=as.factor(track), 
                                                        ymin=velocity_mmpersec_log10, 
                                                        ymax=velocity_mmpersec_log10, 
                                                        color=as.factor(isachase)), 
                position=position_dodge(0.85), width=0.5)+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Log10 Velocity (mm/sec)")+
  xlab("Fish Identity")

#coming back here after doing a random selection, I think we need the quartiles, really:
annotatedchasefull_quartiles <- annotatedchasedf_full %>%
  group_by(track, isachase) %>%
  summarize(velmin = min(velocity_mmpersec, na.rm=TRUE),
            vel25 = quantile(velocity_mmpersec, probs = 0.25, na.rm=T),
            vel50 = quantile(velocity_mmpersec, probs = 0.5, na.rm=T),
            vel75 = quantile(velocity_mmpersec, probs = 0.75, na.rm=T),
            velmax = max(velocity_mmpersec, na.rm=T),
            )
#save as table maybe?
#this is useful in the future, but I think for now I will save as excel and make table there
#ggtexttable(annotatedchasefull_quartiles, rows = NULL, theme = ttheme("light"))
write.csv(annotatedchasefull_quartiles, file="chasequartiles.csv")
#melt + ggplot visualization
annotatedchasefull_quartiles <- reshape2::melt(annotatedchasefull_quartiles, id=c("track", "isachase"),variable.name="quartile")
ggplot()+
  geom_point(annotatedchasefull_quartiles, mapping=aes(x=as.factor(track), 
                                                       y=value, 
                                                       shape=as.factor(isachase),
                                                       color=as.factor(isachase)),
             size=2, position=position_dodge(width=0.9))+
  geom_violin(annotatedchasedf_full, mapping=aes(x=as.factor(track),
                                              y=velocity_mmpersec,
                                              fill=as.factor(isachase)),
              alpha=0.2, position="dodge")+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Velocity (mm/sec) Quartiles")+
  xlab("Fish Identity")+
  scale_y_continuous(trans='log10')

ggplot(annotatedchasedf_full, aes(x=as.factor(track), y=velocity_mmpersec, fill=as.factor(isachase)))+
  geom_violin(alpha=0.2)+
  geom_point(annotatedchasefull_quartiles, mapping=aes(x=as.factor(track), y=quartile, color=as.factor(isachase)))+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Velocity (mm/sec)")+
  xlab("Fish Identity")
#also, so I think it may be helpful to do another like, random selection of times
#from times the fish are not chasing
#and compare those velocity curves with chases right?

#random list of timestamps for non-chases:
#setseed for random generation
set.seed(1722)
round(runif(n=5, min=1, max=44983), 0) #nochases
set.seed(1717)
round(runif(n=5, min=1, max=81),0) #chases

#nochases are: c(3562,20521,786,31792,17923)
#chases are: c(69, 13, 50, 54, 67)

#here we do filter: 
annotatedchasedf_full_marked <- annotatedchasedf_full
annotatedchasedf_full_marked[annotatedchasedf_full_marked$frame %in% c(3562,20521,786,31792,17923),7] <- "N"
#checked in the terminal, we're marked up. So now we filter for +/- 1 sec before and after random times:
annotatedchasedf_marks <- annotatedchasedf_full_marked[leadlag(annotatedchasedf_full_marked$Behavior == "N", 25, 25),]

#pulling code to annotate chaseids to annotate markids:
annotatedchasedf_marks <- annotatedchasedf_marks[-c(11,12)]
tmp <- annotatedchasedf_marks %>%
  group_by(track) %>%
  filter(Behavior == "N") %>%
  mutate(MarkID = seq_along(Behavior))
annotatedchasedf_marks <- merge(x = annotatedchasedf_marks, 
                          y = tmp, 
                          by = c("frame", "nodetrackid", "nodes", "track", "x", "y", "Behavior", 
                                 "distance", "velocity_mmpersec", "second"),
                          all.x = TRUE)
annotatedchasedf_marks <- annotatedchasedf_marks %>%
  arrange(track, frame)
annotatedchasedf_marks$MarkID <- roll(annotatedchasedf_marks$MarkID, 25)
rm(tmp, annotatedchasedf_full_marked)

#that took some time but now we can spit out deez figures:
marklist <- split(annotatedchasedf_marks, annotatedchasedf_marks$MarkID)
plotmark <- function(n){
  ggplot(n, aes(x=second, y=velocity_mmpersec, group=MarkID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
}
markplotlist <- lapply(marklist, plotmark)
marknames <- paste("mark", unique(annotatedchasedf_marks$MarkID), sep='-')
names(markplotlist) <- marknames
lapply(names(markplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""), plot=markplotlist[[x]], path=here("WT031623_26dpf_Figures")))

#

  