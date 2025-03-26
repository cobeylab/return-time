library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_canada_resp.R")

data_canada_resp_scaled2 <- data_canada_resp_scaled %>%
  mutate(
    key=factor(key,
               levels=c("AdV", "HMPV", "PIV", "RV/EV", "RSV",
                        "Flu A", "Flu B", "CoV", "noro"),
               labels=c("Adenovirus", "Human metapneumovirus",
                        "Parainfluenza virus", "Rhinovirus/Enterovirus",
                        "RSV", "Flu A", "Flu B", "Human coronavirus",
                        "Norovirus"))
  ) %>% 
  filter(key %in% c("Adenovirus", "Human coronavirus", "Human metapneumovirus",
                    "Parainfluenza virus", "Rhinovirus/Enterovirus",
                    "RSV"))

g1 <- ggplot(data_canada_resp_scaled2) +
  geom_line(aes(year+week/52, tests)) +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Weekly tests", expand=c(0, 0), limits=c(0, 3.9e4)) +
  facet_wrap(~key) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

ggsave("figure_test_canada.pdf", g1, width=10, height=4)
