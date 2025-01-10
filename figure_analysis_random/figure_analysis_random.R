library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
load("../analysis_random/analysis_random_window.rda")
load("../analysis_random/analysis_random_simple.rda")
load("../analysis_random/analysis_random_iterative.rda")

analysis_random_all <- bind_rows(
  analysis_random_iterative %>% mutate(method="Iterative"),
  analysis_random_simple %>% mutate(method="Naive"),
  analysis_random_window %>% mutate(method="Window-based")
)

g1 <- ggplot(analysis_random_all) +
  geom_point(aes(resilience_true, est, col=duration)) +
  geom_smooth(aes(resilience_true, est), method="lm", col="red", fill="red") +
  geom_abline(intercept=0, slope=1, lty=2) +
  scale_x_continuous("Intrinsic resilience") +
  scale_y_continuous("Estimated resilience") +
  scale_color_viridis_c("Duration of NPIs") +
  facet_wrap(~method) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )

ggsave("figure_analysis_random.pdf", g1, width=6, height=4)

analysis_random_all_filter <- analysis_random_all %>%
  filter(est > 0, duration < 2)

lapply(split(analysis_random_all_filter, analysis_random_all_filter$method),
       function(x) {
         cor(x$resilience_true, x$est)
       })

ggplot(analysis_random_all_filter) +
  geom_point(aes(resilience_true, est, col=duration)) +
  geom_smooth(aes(resilience_true, est), method="lm", col="red", fill="red") +
  geom_abline(intercept=0, slope=1, lty=2) +
  scale_x_log10("Intrinsic resilience") +
  scale_y_log10("Estimated resilience") +
  scale_color_viridis_c("Duration of NPIs") +
  facet_wrap(~method) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )
