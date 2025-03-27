library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_nrevss.R")

data_nrevss_resp_comb2 <- data_nrevss_resp_comb %>%
  filter(type %in% c("Adenovirus", "PIV", "CoV", "Human metapneumovirus", "RSV",
                     "Rhinovirus"),
         year >= 2014)

g1 <- ggplot(data_nrevss_resp_comb2) +
  geom_line(aes(year+week/52, tests)) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2014:2024) +
  scale_y_continuous("Total cases", expand=c(0, 0)) +
  facet_wrap(~type) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

ggsave("figure_test_nrevss.pdf", g1, width=10, height=4)
