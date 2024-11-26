## Chase Analysis: Individual Variation

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)

#based on bfd generated in 'behaviorextraction'

#palette
palette <- c('0' = "#307CAB", '1' = "#C1774B", '2' = "#CCB22B", '3' = "#86628C", '4' = "#8FAF61")

#Big ol figure generator (everything):

behaviorlist <- split(bdf, bdf$BehaviorID)

arrangeplots <- function (n){
  q <- n %>% 
    group_by(BehaviorID) %>%
    summarize(quantile = quantile(velocity_mmpersec, probs=c(0.975)))
  p1 = ggplot(n, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
    geom_line()+
    geom_hline(q, mapping=aes(yintercept=quantile), color="red")+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p2 = ggplot(n, aes(x=second, y=acceleration_mmpsps, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p3 = ggplot(n, aes(x=second, y=theta, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p4 =  ggplot(n)+
      geom_path(mapping=aes(x=second, y=trueheadingangle, group=track, color=as.factor(track)))+
      theme_classic()+
      theme(legend.position="none")+
    scale_color_manual(values=c(palette))
  
  p5 = ggplot(n, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
    geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25, alpha=0.3)+
    geom_path(linewidth=0.2)+
    geom_point(size=0.3)+
    theme_classic()+
    ylab("Vector Y Component")+
    xlab("Vector X Component")+
    labs(color="Individual")+
    scale_color_manual(values=palette)+
    theme(text=element_text(size=10),
          legend.position = c(0.95, 0.9),
          legend.key.height = unit(0.25, 'cm'),
          legend.background = element_rect(fill="transparent"))

  p6 = ggplot(n, aes(x=second, y=truenn, fill=truenndist))+
    geom_tile()+
    theme_classic()+
    ylab("Fish Combination")+
    xlab("Time (Seconds)")+
    labs(fill="Dist Btwn\nFish (mm)")+
    theme(text=element_text(size=15),
          legend.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.position = c(0.95, 0.9),
          legend.key.height = unit(0.25, 'cm'),
          legend.background = element_rect(fill="transparent"))
  
  p7 =  ggplot(n, aes(x=x, y=y, color=as.factor(track)))+
    geom_point(size=1)+
    geom_path(linewidth=0.5)+
    theme_classic()+
    scale_color_manual(values=palette)+
    coord_cartesian(xlim=c(0, 160), ylim=c(0, 125))
  
  p8 = ggplot()+
    geom_line(n, mapping=aes(x=angletofish0, y=truenndist), color="#307CAB")+
    geom_line(n, mapping=aes(x=angletofish1, y=truenndist), color="#C1774B")+
    geom_line(n, mapping=aes(x=angletofish2, y=truenndist), color="#CCB22B")+
    geom_line(n, mapping=aes(x=angletofish3, y=truenndist), color="#86628C")+
    geom_line(n, mapping=aes(x=angletofish4, y=truenndist), color="#8FAF61")+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  finalplot = ggarrange(
    ggarrange(p1, p2, p3, ncol=3, labels = c("Velocity (MMPS)", "Acceleration (MMPSPS)", "Tail Angle (Deg)")), 
    ggarrange(p6, p4, p8, ncol=3, labels =c("Nearest Neighbor Dist (MM)", "True Heading (Deg)", "Relative Heading (Deg)")),
    ggarrange(p5, p7, ncol=2, labels =c("Vector Field", "Coordinate Plot"), widths=1, heights=1),
    nrow=3)
}

arrangedplotlist <- lapply(behaviorlist, arrangeplots)
plotnames <- paste("chase", unique(bdf$BehaviorID), sep='-')
names(arrangedplotlist) <- plotnames
lapply(names(arrangedplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""),
                          plot=arrangedplotlist[[x]], 
                          width = 80,
                          height = 80,
                          unit = "cm",
                          path=here("Figures")))

# #the bare essentials (velocity, acceleration, nn, position and vectors)
# 
# arrangeplots_abbrev <- function (n){
#   p1 = ggplot(n, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
#     geom_line()+
#     theme_classic()+
#     facet_wrap(~track, ncol=1)
#   
#   p2 = ggplot(n, aes(x=second, y=acceleration_mmpsps, group=BehaviorID))+
#     geom_line()+
#     theme_classic()+
#     facet_wrap(~track, ncol=1)
#   
#   p3 = ggplot(n, aes(x=second, y=truenn, fill=truenndist))+
#     geom_tile()+
#     theme_classic()+
#     ylab("Fish Combination")+
#     xlab("Time (Seconds)")+
#     labs(fill="Dist Btwn\nFish (mm)")+
#     theme(text=element_text(size=15),
#           legend.title = element_text(size=10),
#           legend.text = element_text(size=8),
#           legend.position = c(0.95, 0.9),
#           legend.key.height = unit(0.25, 'cm'),
#           legend.background = element_rect(fill="transparent"))
#   
#   p4 =  ggplot(n, aes(x=x, y=y, color=as.factor(track)))+
#     geom_point(size=1)+
#     geom_path(linewidth=0.5)+
#     theme_classic()+
#     scale_color_manual(values=palette)+
#     coord_cartesian(xlim=c(0, 160), ylim=c(0, 125))
#   
#   p5 = ggplot(n, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
#     geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25, alpha=0.3)+
#     geom_path(linewidth=0.2)+
#     geom_point(size=0.3)+
#     theme_classic()+
#     ylab("Vector Y Component")+
#     xlab("Vector X Component")+
#     labs(color="Individual")+
#     scale_color_manual(values=palette)+
#     theme(text=element_text(size=10),
#           legend.position = c(0.95, 0.9),
#           legend.key.height = unit(0.25, 'cm'),
#           legend.background = element_rect(fill="transparent"))
#   
#   finalplot = ggarrange(
#     ggarrange(p1, p2, p3, ncol=3, labels = c("Velocity (MMPS)", "Acceleration (MMPSPS)", "Relative Nearest Neighbor (MM)")),
#     ggarrange(p4,p5, ncol=2, labels =c("Coordinate Plot", "Vector Field")),
#     nrow=2)
# }
# arrangedplotlist <- lapply(behaviorlist, arrangeplots_abbrev)
# plotnames <- paste("chase-essentials", unique(bdf$BehaviorID), sep='-')
# names(arrangedplotlist) <- plotnames
# lapply(names(arrangedplotlist), 
#        function(x) ggsave(filename=paste(x,".png",sep=""),
#                           plot=arrangedplotlist[[x]], 
#                           width = 60,
#                           height = 75,
#                           unit = "cm",
#                           path=here("Figures")))

###################################################################################

#individual checks

arrangeplots(chase26)

chase26 <- bdf %>%
  filter(BehaviorID == 26)

#idk if I'll use it but this is a really nice palette, for future:

palette <- c('0' = "#307CAB", '1' = "#C1774B", '2' = "#CCB22B", '3' = "#86628C", '4' = "#8FAF61")
 
 #just a quick velocity test:
 ggplot(chase26, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
   geom_line()+
   theme_classic()+
   facet_wrap(~track, ncol=1)
 
 #nn
 ggplot(chase26, aes(x=second, y=truenn, fill=truenndist))+
   geom_tile()+
   theme_classic()+
   ylab("Fish Combination")+
   xlab("Time (Seconds)")+
   labs(fill="Dist Btwn\nFish (mm)")+
   theme(text=element_text(size=15),
         legend.title = element_text(size=10),
         legend.text = element_text(size=8),
         legend.position = c(0.95, 0.9),
         legend.key.height = unit(0.25, 'cm'),
         legend.background = element_rect(fill="transparent"))
 
 #xycoordinate:
 ggplot(chase26, aes(x=x, y=y, color=as.factor(track)))+
   geom_point(size=1)+
   geom_path(linewidth=0.5)+
   theme_classic()+
   scale_color_manual(values=palette)+
   coord_cartesian(xlim=c(0, 160), ylim=c(0, 125))
 
 #reworkingheading angle...
 ggplot()+
   geom_line(chase26, mapping=aes(x=angletofish0, y=truenndist), color="#307CAB")+
   geom_line(chase26, mapping=aes(x=angletofish1, y=truenndist), color="#C1774B")+
   geom_line(chase26, mapping=aes(x=angletofish2, y=truenndist), color="#CCB22B")+
   geom_line(chase26, mapping=aes(x=angletofish3, y=truenndist), color="#86628C")+
   geom_line(chase26, mapping=aes(x=angletofish4, y=truenndist), color="#8FAF61")+
   theme_classic()+
   coord_polar()+
   facet_wrap(~track, ncol=1)
 
#vector field component
ggplot(chase11, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
   geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25, alpha=0.3)+
   geom_path(linewidth=0.2)+
   geom_point(size=0.3)+
   theme_classic()+
   ylab("Vector Y Component")+
   xlab("Vector X Component")+
   labs(color="Individual")+
   scale_color_manual(values=palette)+
   theme(text=element_text(size=10),
         legend.position = c(0.95, 0.9),
         legend.key.height = unit(0.25, 'cm'),
         legend.background = element_rect(fill="transparent"))
 
 #still need to get rid of the 0 bin but other than that this is workable.
 ggplot(chase26)+
   geom_histogram(mapping=aes(x=angletofish0, fill="Angle to fish 0"), binwidth=5, fill="#307CAB")+
   geom_histogram(mapping=aes(x=angletofish1, fill="Angle to fish 1"), binwidth=5, fill="#C1774B")+
   geom_histogram(mapping=aes(x=angletofish2, fill="Angle to fish 2"), binwidth=5, fill="#CCB22B")+
   geom_histogram(mapping=aes(x=angletofish3, fill="Angle to fish 3"), binwidth=5, fill="#86628C")+
   geom_histogram(mapping=aes(x=angletofish4, fill="Angle to fish 4"), binwidth=5, fill="#8FAF61")+
   coord_polar()+
   theme_bw()+
   theme(legend.position = c(0.8, 0.28),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x = element_text(size=8))+
   facet_wrap(~track)

 #vector field old:
 ggplot(chase10, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
   geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25)+
   theme_classic()+
   ylab("Vector Y Component")+
   xlab("Vector X Component")+
   labs(color="Individual")+
   theme(text=element_text(size=10),
         legend.position = c(0.95, 0.9),
         legend.key.height = unit(0.25, 'cm'),
         legend.background = element_rect(fill="transparent"))

 #nearest neighbor:
 ggplot(chase10, aes(x=second, y=truenn, fill=truenndist))+
   geom_tile()+
   theme_classic()+
   ylab("Fish Combination")+
   xlab("Time (Seconds)")+
   labs(fill="Dist Btwn\nFish (mm)")+
   theme(text=element_text(size=15),
         legend.title = element_text(size=10),
         legend.text = element_text(size=8),
         legend.position = c(0.95, 0.9),
         legend.key.height = unit(0.25, 'cm'),
         legend.background = element_rect(fill="transparent"))
 
 ########################################################################
 
 #individual aspects if we just want to look at those
 
 plotheading <- function(n) {
   ggplot(n)+
     geom_line(mapping=aes(x=second, y=angletofish0), color="#307CAB")+
     geom_line(mapping=aes(x=second, y=angletofish1), color="#C1774B")+
     geom_line(mapping=aes(x=second, y=angletofish2), color="#CCB22B")+
     geom_line(mapping=aes(x=second, y=angletofish3), color="#86628C")+
     geom_line(mapping=aes(x=second, y=angletofish4), color="#8FAF61")+
     theme_classic()+
     facet_wrap(~track, ncol=1)
 }
 behaviorlist <- split(bdf, bdf$BehaviorID)
 plotheadinglist <- lapply(behaviorlist, plotheading)
 plotnames <- paste("headingangle", unique(bdf$BehaviorID), sep='-')
 names(plotheadinglist) <- plotnames
 lapply(names(plotheadinglist), 
        function(x) ggsave(filename=paste(x,".png",sep=""),
                           plot=plotheadinglist[[x]], 
                           width = 12,
                           height = 10,
                           unit = "cm",
                           path=here("Figures")))
 #bring back this as a polar histogram?
 
 plotheadingpolar <- function(n) {
   ggplot(n)+
     geom_histogram(mapping=aes(x=angletofish0), fill="#307CAB", alpha = 0.3)+
     geom_histogram(mapping=aes(x=angletofish1), fill="#C1774B", alpha = 0.3)+
     geom_histogram(mapping=aes(x=angletofish2), fill="#CCB22B", alpha = 0.3)+
     geom_histogram(mapping=aes(x=angletofish3), fill="#86628C", alpha = 0.3)+
     geom_histogram(mapping=aes(x=angletofish4), fill="#8FAF61", alpha = 0.3)+
     scale_x_continuous(breaks = seq(0, 360, 60)) +
     theme_classic()+
     coord_polar() +
     xlab("Relative to Fish")+
     ylab("Heading Angle")+
     facet_wrap(~track, ncol=2)
 }
 behaviorlist <- split(bdf, bdf$BehaviorID)
 plotheadingpolarlist <- lapply(behaviorlist, plotheadingpolar)
 plotnames <- paste("headingangle", unique(bdf$BehaviorID), sep='-')
 names(plotheadingpolarlist) <- plotnames
 lapply(names(plotheadingpolarlist), 
        function(x) ggsave(filename=paste(x,".png",sep=""),
                           plot=plotheadingpolarlist[[x]], 
                           path=here("Figures")))
 
 
 plotrueheading <- function(n) {
   ggplot(n)+
     geom_path(mapping=aes(x=second, y=trueheadingangle, group=track, color=as.factor(track)))+
     theme_classic()+
     theme(legend.position="none")
 }
 behaviorlist <- split(bdf, bdf$BehaviorID)
 plottrueheadinglist <- lapply(behaviorlist, plotrueheading)
 plotnames <- paste("trueheadingangle", unique(bdf$BehaviorID), sep='-')
 names(plottrueheadinglist) <- plotnames
 lapply(names(plottrueheadinglist), 
        function(x) ggsave(filename=paste(x,".png",sep=""),
                           plot=plottrueheadinglist[[x]], 
                           width = 10,
                           height = 12,
                           unit = "cm",
                           path=here("Figures")))
 
 ####
 
 #for this function we'll need to caluclate the quantiles, and plot those on top of velocity:
 
 plotvelocity <- function(n) {
   q <- n %>% 
     group_by(BehaviorID) %>%
     summarize(quantile = quantile(velocity_mmpersec, probs=c(0.98)))
   ggplot(n)+
     geom_path(mapping=aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
     theme_classic()+
     geom_hline(q, mapping=aes(yintercept=quantile), color="red")+
     facet_wrap(~track, ncol=1)
 }
 
behaviorlist <- split(bdf, bdf$BehaviorID)
plotvelocitylist <- lapply(behaviorlist, plotvelocity)
plotnames <- paste("chasevelocity_98", unique(bdf$BehaviorID), sep='-')
names(plotvelocitylist) <- plotnames
lapply(names(plotvelocitylist), 
        function(x) ggsave(filename=paste(x,".png",sep=""),
                           plot=plotvelocitylist[[x]],
                           path=here("Figures")))
