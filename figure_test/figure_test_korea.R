library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_korea.R")

g1 <- ggplot(data_korea_scale) +
  geom_line(aes(year+week/52, total)) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2014:2024) +
  scale_y_continuous("Total cases", expand=c(0, 0), limits=c(0, 3.2e3)) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

ggsave("figure_test_korea.pdf", g1, width=10, height=4)
