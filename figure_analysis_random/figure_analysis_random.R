library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
load("../analysis_random/analysis_random_window.rda")
load("../analysis_random/analysis_random_simple.rda")

analysis_random_window_summ <- analysis_random_window %>%
  group_by(cut, div, R) %>%
  summarize(
    cor=cor(resilience_true, est, use="complete.obs")
  ) %>%
  mutate(
    R=paste0("R=",R),
    R=factor(R, levels=paste0("R=",(2:7)*2)),
    cut=paste0("Truncation threshold=", cut)
  )

analysis_random_simple_summ <- analysis_random_simple %>%
  group_by(R) %>%
  summarize(
    cor=cor(resilience_true, est, use="complete.obs")
  ) %>%
  mutate(
    R=paste0("R=",R),
    R=factor(R, levels=paste0("R=",(2:7)*2))
  )

analysis_random_window_summ %>%
  filter(cor > 0.6)

g1 <- ggplot(analysis_random_window_summ) +
  geom_point(aes(div, cor)) +
  geom_hline(data=analysis_random_simple_summ, aes(yintercept=cor), lty=2,
             col="red") +
  scale_x_continuous("Number of divisions") +
  facet_grid(cut~R) +
  theme(
    strip.background = element_blank()
  )

ggsave("figure_analysis_random.pdf", g1, width=8, height=6)
