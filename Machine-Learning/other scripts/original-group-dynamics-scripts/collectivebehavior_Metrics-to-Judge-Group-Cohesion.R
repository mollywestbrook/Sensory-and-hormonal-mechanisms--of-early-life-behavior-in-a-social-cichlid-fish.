## Metrics to juge group cohesion
#Again, based in DLC

#note: these are figures and summary statistics only! All other

#Note: for metrics 1-4, we do not need fish identity
#for metrics 5-7, we do (as it's based on speed and speed will be jumbing around if the identities are flip flopping)

trialdata <- read.csv(here("formatted_05-22-2023_group5_L.csv"))

trialdata_head <- trialdata %>%
  filter(nodes == "head")

#to start, we'll do a sampled behavioral trace

trialdata_framesample <- trialdata_head %>%
  filter(frame%%100 == 0)

foursecondtrace <- ggplot(trialdata_framesample, aes(x, y, color=as.factor(track)))+
  geom_point(size=1, alpha=0.3)+
  geom_path(size=1, alpha=0.3)+
  ylab("X Position (mm)")+
  xlab("Y Position (mm)")+
  theme(text = element_text(size = 20))+   
  #coord_flip()+
  theme_classic()
foursecondtrace
ggsave(foursecondtrace, file="foursecondtrace_052223_G5_L.png", path=here("Figures"))

#Single Statistics:

#1: Group Spacing (median of the mean pairwise distances across whole group)

groupspacing <- median(trialdata_coordinatedistance$distavg, na.rm=TRUE)

#maybe normalize to size of container: 

groupspacing <- groupspacing/456.96 #tank
#groupspacing <- groupspacing/157 #dish

#2: Polarization 

unitvelocity

#3: Nearest Neighbor distance (median of the closest neighbor for each fish at each time step)

nnd <- median(trialdata_coordinatedistance$distmin, na.rm=TRUE)

#normalize

nnd <- nnd/456.96 #tank
#nnd <- nnd/157 #dish

#4: Centroid speed (median of centroid of all fish centroids, speed over time)

centroidspeed <- median(trialdata$centroidspeed, na.rm=TRUE)

#5: Fraction of Time Moving (need a moving window average of 10 sec, whenever this value is greater than 0.05 quantile out of whole)

fractiontimemoving <- (sum(velocitydf > quants$cutoff, na.rm=T)) / ((nrow(velocitydf)*ncol(velocitydf)))

#6: Speed while moving (median speed when it's above threshold of #5)

speedwhilemoving <- median(velocitydf > quants$cutoff, na.rm=TRUE)

#7: Speed IQR (difference between 0.75 ad 0.25 quartile speeds of speeds above threshold of #5)

speedIQR <- quants$seventyfifth - quants$twentyfifth

#and assemble into df:

statisticsummary <- data.frame(groupspacing, nnd, fractiontimemoving, centroidspeed, speedwhilemoving, speedIQR)
statisticsummary <- reshape2::melt(statisticsummary, measure.vars = c("groupspacing", "nnd", "fractiontimemoving", "centroidspeed", "speedwhilemoving", "speedIQR"))
fwrite(statisticsummary, file = "statisticsummary.csv", row.names=FALSE)

#and then I also want to make a little bar graph figure that shows the summary statistics:

statisticsummaryfigure <- ggplot()+ 
  geom_bar(data=statisticsummary, mapping=aes(x=value, y=as.factor(variable)), stat = "identity", position= "dodge") + 
  geom_vline(xintercept=1, linetype="dashed")+
  ylab("Statistic")+
  xlab("Value")+
  theme_classic() +
  theme(text = element_text(size = 18))
statisticsummaryfigure
ggsave(statisticsummaryfigure, file="statisticsummaryfigure.png", path=here("Figures"))

##Figures

##Group spacing (mean pairwise distance across whole group)

groupspacingplot <- ggplot(trialdata_coordinatedistance, aes (x=Second, y=distavg))+
  geom_line(size=0.4)+
  ylab("Average Distance Between Animals (mm)")+
  xlab("Time (Seconds)")+
  theme(text = element_text(size = 20))+   
  theme_classic()
groupspacingplot
ggsave(groupspacingplot, file="groupspacingplot.png", path=here("Figures"))

#centroid trace:

centroidtrace <- ggplot(trialdata, aes(x=groupcentroidx, y=groupcentroidy))+
  geom_path(size=0.4)+
  xlab("Centroid X Position (mm)")+
  ylab("Centroid Y Position (mm)")+
  theme(text = element_text(size = 20))+   
  #coord_flip()+
  theme_classic()
centroidtrace
ggsave(centroidtrace, file="centroidtrace.png", path=here("Figures"))

#Centroid speed

centroidspeedovertime <- ggplot(trialdata, aes(x=Second, y=centroidspeed))+
  geom_line(size=0.4)+
  xlab("Time (Seconds)")+
  ylab("Group Centroid Speed (mm/frame)")+
  theme(text = element_text(size = 20))+
  theme_classic()
centroidspeedovertime
ggsave(centroidspeedovertime, file="centroidspeedovertime.png", path=here("Figures"))
  
##NND over time (we'll call this the minimum distance between all members of the total group, not a per individual basis)

minimumdistanceplot <- ggplot(trialdata_coordinatedistance, aes (x=Second, y=distmin))+
  geom_line(size=0.4)+
  ylab("Minimum Distance Between 2 Animals (mm)")+
  xlab("Time (Seconds)")+
  theme(text = element_text(size = 20))+   
  theme_classic()
minimumdistanceplot
ggsave(minimumdistanceplot, file="minimumdistanceplot.png", path=here("Figures"))