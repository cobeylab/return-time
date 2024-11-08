library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_korea/analysis_korea_growth.rda")

analysis_korea_return <- analysis_korea_growth %>%
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

ggplot(analysis_korea_return) +
  geom_violin(aes(key, value, fill=method)) +
  scale_y_continuous("Return time (years)") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "top"
  )
