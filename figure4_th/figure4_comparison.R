library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
load("figure4_summ_th.rda")
load("../figure4/figure4_summ.rda")

analysis_all_summ_filter_th2 <- analysis_all_summ_filter_th %>%
  select(key, country, resilience, resilience_lwr, resilience_upr) %>%
  rename(
    resilience_th=resilience, 
    resilience_th_lwr=resilience_lwr, 
    resilience_th_upr=resilience_upr
  )

analysis_all_summ_filter2 <- analysis_all_summ_filter %>%
  select(key, country, resilience, resilience_lwr, resilience_upr)

analysis_comb <- merge(analysis_all_summ_filter2, analysis_all_summ_filter_th2)

g1 <- ggplot(analysis_comb) +
  geom_abline(intercept=0, slope=1, lty=2) +
  geom_point(aes(resilience, resilience_th, col=country, shape=country)) +
  geom_errorbar(aes(x=resilience, ymin=resilience_th_lwr, ymax=resilience_th_upr, col=country), width=0) +
  geom_errorbarh(aes(xmin=resilience_lwr, xmax=resilience_upr, y=resilience_th, col=country), height=0) +
  scale_x_continuous("Original resilience estimates (1/year)") +
  scale_y_continuous("Resilience estimates using higher dimensions (1/year)") +
  scale_color_viridis_d("Country") +
  scale_shape_discrete("Country") +
  theme(
    
  )

ggsave("figure4_comparison.pdf", g1, width=6, height=4)
