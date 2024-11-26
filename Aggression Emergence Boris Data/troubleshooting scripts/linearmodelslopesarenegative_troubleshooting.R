#Stats troubleshooting

#use in conjunction with sex steroid analysis

#control

hormonedf_DMSO <- hormonedf %>%
  filter(Treatment == "DMSO")

ggplot(hormonedf_DMSO, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(ID)))+
  geom_point(size=2, alpha=0.6)+
  geom_path(size=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  geom_abline(intercept = summary$coefficients[1,1], slope = summary$coefficients[2,1], color="green", linewidth=1)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#004D40", "#D81B60", "#1E88E5"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#EE

hormonedf_EE <- hormonedf %>%
  filter(Treatment == "EE")

ggplot(hormonedf_EE, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(ID)))+
  geom_point(size=2, alpha=0.6)+
  geom_path(size=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  geom_abline(intercept = summary$coefficients[3,1], slope = summary$coefficients[5,1], color="red", linewidth=1)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#e06d97"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")

#KT

hormonedf_KT <- hormonedf %>%
  filter(Treatment == "KT")

ggplot(hormonedf_KT, aes(x=Age, y=Behaviors, color=as.factor(Treatment), group=as.factor(ID)))+
  geom_point(size=2, alpha=0.6)+
  geom_path(size=1, stat="identity", na.rm=TRUE, alpha=0.5)+
  geom_abline(intercept = summary$coefficients[4,1], slope = summary$coefficients[6,1], color="blue", linewidth=1)+
  theme_classic()+
  theme(text = element_text(size = 15), 
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  guides(color=guide_legend(title="Treatment"))+
  scale_color_manual(values=c("#1E88E5"))+
  ylab("Aggressive Behavior Count")+
  xlab("Age (Days Post Fertilization)")
