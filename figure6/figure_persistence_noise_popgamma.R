library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
source("../R/simulate_stochastic.R")

load("../analysis_factorial/analysis_factorial_noise_pop.rda")
load("../analysis_factorial/analysis_factorial_noise_gamma.rda")

g1 <- ggplot(analysis_factorial_noise_pop) +
  geom_raster(aes(N, 1/delta/365, fill=log10(resilience*365))) +
  scale_x_log10("Population size", expand=c(0, 0)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Resilience (1/years)",
                       breaks=log10(c(0.03, 0.1, 0.3)),
                       labels=c(0.03, 0.1, 0.3)) +
  theme(
    legend.position = "bottom"
  )

g2 <- ggplot(analysis_factorial_noise_pop %>% filter(!extinction)) +
  geom_raster(aes(N, 1/delta/365, fill=log10(amplitude))) +
  scale_x_log10("Population size", expand=c(0, 0)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Epidemic amplitude", option="A",
                       breaks=log10(c(0.01, 0.03, 0.1, 0.3, 1)),
                       labels=c(0.01, 0.03, 0.1, 0.3, 1)) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

g3 <- ggplot(analysis_factorial_noise_gamma) +
  geom_raster(aes(1/gamma, 1/delta/365, fill=log10(resilience*365))) +
  scale_x_log10("Duration of infection (days)", expand=c(0, 0),
                breaks=c(7, 14, 28)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Resilience (1/years)",
                       breaks=log10(c(0.03, 0.1, 0.3)),
                       labels=c(0.03, 0.1, 0.3)) +
  theme(
    legend.position = "bottom"
  )

g4 <- ggplot(analysis_factorial_noise_gamma) +
  geom_raster(aes(1/gamma, 1/delta/365, fill=log10(amplitude))) +
  scale_x_log10("Duration of infection (days)", expand=c(0, 0),
                breaks=c(7, 14, 28)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Epidemic amplitude", option="A",
                       breaks=log10(c(0.01, 0.03, 0.1, 0.3, 1)),
                       labels=c(0.01, 0.03, 0.1, 0.3, 1)) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

gcomb <- ggarrange(g1, g2, g3, g4, nrow=2,
                   labels=c("A", "B", "C", "D"))

ggsave("figure_persistence_noise_popgamma.pdf", gcomb, width=8, height=8)
