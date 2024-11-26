### Aggression Clip Analysis for GRAND Talk:

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)

ggplot(featuredf, aes(x=second, y=velocity_mmpersec))+
  geom_line()+
  theme_classic()+
  facet_wrap(~track, ncol=1)

chase <- featuredf %>%
  filter(frame %in% c(330:390))

ggplot(chase, aes(x=second, y=velocity_mmpersec))+
  geom_line()+
  theme_classic()+
  facet_wrap(~track, ncol=1)+
  ylab("Velocity (mm/sec)")

ggplot(chase, aes(x=second, y=acceleration_mmpsps))+
  geom_line()+
  theme_classic()+
  facet_wrap(~track, ncol=1)+
  ylab("Acceleration (mm/sec^2)")

ggplot(chase, aes(x=second, y=truenn, fill=truenndist))+
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