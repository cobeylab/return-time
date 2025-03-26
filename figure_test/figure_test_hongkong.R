library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_hongkong_resp.R")
source("../script/script_data_hongkong_piv.R")
source("../script/script_data_hongkong_fecal.R")

g1 <- ggplot(data_hongkong_resp_scaled %>% filter(key=="adeno")) +
  geom_line(aes(year+week/52, test)) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2014:2024) +
  scale_y_continuous("Weekly tests", expand=c(0, 0), limits=c(0, 1.1e4)) +
  ggtitle("Respiratory pathogens") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

g2 <- ggplot(data_hongkong_piv_scaled %>% filter(key=="PIV 1")) +
  geom_line(aes(year+week/52, test)) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2014:2024) +
  scale_y_continuous("Weekly tests", expand=c(0, 0), limits=c(0, 1.1e4)) +
  ggtitle("PIV") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

g3 <- ggplot(data_hongkong_fecal_scaled %>% filter(key=="noro")) +
  geom_line(aes(year+as.numeric(month)/12, test)) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2013:2024) +
  scale_y_continuous("Weekly tests", expand=c(0, 0), limits=c(0, 1.9e3)) +
  ggtitle("Gastroenteritis viruses") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

gcomb <- ggarrange(g1, g2, g3, nrow=1)

ggsave("figure_test_hongkong.pdf", gcomb, width=12, height=4)
