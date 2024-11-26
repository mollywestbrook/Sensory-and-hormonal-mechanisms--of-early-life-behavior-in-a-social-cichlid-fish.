## Theoretical figure

library(here)
library(tidyverse)

#generate some uniform data

dist1 <- rnorm(100, mean=50, sd=5)
dist2 <- rnorm(100, mean=40, sd=5)

df <- data.frame(dist1, dist2)
df$ID <- rownames(df)

ggplot()+
  geom_density(df, mapping=aes(x=dist1), color="#e90b00", fill="#e90b00", alpha=0.2, bw=4)+
  geom_density(df, mapping=aes(x=dist2), color="#002eff", fill="#002eff", alpha=0.2, bw=4)+
  geom_segment(mapping = aes(x = 35, xend = 35, y = 0, yend = 0.047), color="#6f88ff", linewidth=1)+
  geom_segment(mapping = aes(x = 50, xend = 50, y = 0, yend = 0.015), color="#001781", linewidth=1)+
  geom_segment(mapping = aes(x = 43, xend = 43, y = 0, yend = 0.0345), color="#820600", linewidth=1)+
  geom_segment(mapping = aes(x = 60, xend = 60, y = 0, yend = 0.016), color="#e4827d", linewidth=1)+
  theme_classic()+
  theme(text = element_text(size=16))+
  ylab("Proportion of Individuals")+
  xlab("Aggressive Behaviors")