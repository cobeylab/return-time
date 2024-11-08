library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_korea/analysis_korea_growth.rda")

summ_korea_resilience <- analysis_korea_growth %>%
  group_by(key, id) %>%
  summarize(
    return_rC=unique(return_rC),
    return_takens=unique(return_takens)
  ) %>%
  gather(
    method, value, -id, -key
  ) %>%
  mutate(
    method=factor(method,
                  levels=c("return_rC", "return_takens"),
                  labels=c("growth rate", "Takens' theorem"))
  )

summ_korea_when <- analysis_korea_growth %>%
  group_by(key, id) %>%
  summarize(
    when_rC=unique(when_rC),
    when_takens=unique(when_takens)
  ) %>%
  gather(
    method, value, -id, -key
  ) %>%
  mutate(
    method=factor(method,
                  levels=c("when_rC", "when_takens"),
                  labels=c("growth rate", "Takens' theorem"))
  )

ggplot(summ_korea_resilience) +
  geom_violin(aes(key, 1/value, fill=method)) +
  scale_y_continuous("Resilience (1/years)") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "top"
  )
