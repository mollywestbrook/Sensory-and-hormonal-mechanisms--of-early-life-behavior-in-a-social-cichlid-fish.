#Model

library(tidyverse)
library(here)
library(data.table)

model <- read_csv(here("modeldataformodel.csv"), col_names = TRUE)

receptors <- model %>%
  filter(Type == "Receptor")

titers <- model %>%
  filter(Type == "Titer")

ggplot(model, aes(x=Type, y=Level, fill=Predisposition))+
  geom_boxplot(coef=0, outlier.shape = NA, position=position_dodge(width=0.8))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=20))+
  scale_fill_manual(values=c("#4d5656", "#95a5a6"))


  

