## Stack Oflow Lapply question

#build a mini bdf

group <- c(rep(1, 30), rep(2, 30))
animal <- c(rep(c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6)), 2))
velocity <- c(rpois(60, 7))
second <- rep(c(1:6), 10)
bdf <- data.frame(group, animal, velocity, second)

bdf_quantiles <- bdf %>%
  group_by(group) %>%
  summarize(quantile = quantile(velocity, probs=c(0.97)))

#works without quantiles:
plotvelocity <- function(n) {
  ggplot(n)+
    geom_path(mapping=aes(x=second, y=velocity, group=group))+
    theme_classic()+
    facet_wrap(~animal, ncol=1)
}
behaviorlist <- split(bdf, bdf$group)
plotvelocitylist <- lapply(behaviorlist, plotvelocity)

#works, but plots all quantiles over the plot
plotvelocity <- function(n, q) {
  ggplot(n)+
    geom_path(mapping=aes(x=second, y=velocity, group=group))+
    theme_classic()+
    geom_hline(q, mapping=aes(yintercept=quantile))+
    facet_wrap(~animal, ncol=1)
}
behaviorlist <- split(bdf, bdf$group)
plotvelocitylist <- lapply(behaviorlist, plotvelocity, bdf_quantiles)