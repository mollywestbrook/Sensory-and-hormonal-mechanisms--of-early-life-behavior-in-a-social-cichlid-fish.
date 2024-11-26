## pixel to mm conversion 

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)

tracks_final <- read_csv(here("pilotdataset_organized_tracksonly.csv"))

#we'll need frame 6900:

frameofinterest <- tracks_final %>%
  filter(frame == 6900)

ggplot(frameofinterest, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_point(size=0.8)+
  theme_classic()

#we want the distance between fish 3 and fish 4's noses:

x1 <- frameofinterest[29,5]
x2 <- frameofinterest[37,5]
y1 <- frameofinterest[29,6]
y2 <- frameofinterest[37,6]

sqrt((x2-x1)^2+(y2-y1)^2)

#conversion: 359.11/62
#divide by 5.79

#check

n1 <- x1/5.79
m1 <- y1/5.79
n2 <- x2/5.79 
m2 <- y2/5.79

sqrt((n2-n1)^2+(m2-m1)^2)

#and then down here we're going to standardize the axes real quack:

frameofinterest$x <- frameofinterest$x*5.79
frameofinterest$y <- frameofinterest$y*5.79

ggplot(frameofinterest, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_point(size=0.8)+
  coord_cartesian(xlim=c(0, 1300), ylim=c(0,975))+
  theme_classic()