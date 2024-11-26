#Sex Steroid Dosing Figure

library(tidyverse)
library(here)
library(data.table)

sexingresults <- read_csv(here("sexsteroidsexes.csv"), col_names = TRUE)

sexingsummary <- sexingresults %>%
  group_by(Treatment) %>%
  summarize(male = sum(Sex == 'M')/length(Treatment)*100,
            female = sum(Sex == 'F')/length(Treatment)*100)
sexingsummary <- melt(sexingsummary, id='Treatment')

sexingresults <- sexingresults %>%
  mutate(NumericSex = case_when(Sex == 'M' ~ 100,
                                Sex == 'F' ~ 0))
names(sexingsummary) <- c("Treatment", "Sex", "Ratio")

ggplot(sexingsummary, aes(x=Treatment, y=Ratio))+
  geom_bar(mapping=aes(fill=Sex), stat="identity", width=0.6)+
  geom_jitter(sexingresults, mapping=aes(x=Treatment, y=NumericSex), width=0.3, height=1, alpha=0.5, size=2.5)+
  theme_classic()+
  scale_fill_manual(labels=c('Male', 'Female'),
                    values=c('male' = "#1E88E5", 'female'="#D81B60"))+
  theme(text=element_text(size=15))+
  guides(fill=guide_legend(title="Treatment"))+
  xlab("Treatment")+
  ylab("Sex Ratio")

ggplot(sexingresults, mapping=aes(x=Treatment, y=NumericSex))+
  geom_jitter(width=0.2, height=1, alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values=c("#1E88E5", "#D81B60"))+
  theme(text=element_text(size=15))+
  guides(fill=guide_legend(title="Treatment"))+
  xlab("Treatment")+
  ylab("Sex Ratio")

################################################################################

#Fisher tests

kttable <- read_csv(here("ktcontingencytable.csv"), col_names = TRUE)

kttable <- kttable[-c(1)]
fisher.test(kttable)

# Fisher's Exact Test for Count Data

# data:  kttable
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.2509948 4.5452324
# sample estimates:
# odds ratio 
#   1.059694 

eetable <- read_csv(here("eecontingencytable.csv"), col_names = TRUE)

eetable <- eetable[-c(1)]
fisher.test(eetable)
# 
# Fisher's Exact Test for Count Data
# 
# data:  eetable
# p-value = 0.001089
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     2.359522 1133.565731
# sample estimates:
# odds ratio 
#   22.39046 
