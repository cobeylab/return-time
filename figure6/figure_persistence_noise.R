library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
source("../R/simulate_stochastic.R")

load("../analysis_factorial/analysis_factorial_noise.rda")

gamma <- 1/7
N <- 1e8

R01 <- 2
delta1 <- 1/365/40
beta1 <- R01*gamma
S1 <- 1/R01
Istar1 <- delta1 * (1-S1)/(beta1 * S1 + delta1)

R02 <- 2
delta2 <- 1/365
beta2 <- R02*gamma
S2 <- 1/R02
Istar2 <- delta2 * (1-S2)/(beta2 * S2 + delta2)

R03 <- 20
delta3 <- 1/365
beta3 <- R03*gamma
S3 <- 1/R03
Istar3 <- delta3 * (1-S3)/(beta3 * S3 + delta3)

set.seed(101)
ss1 <- simulate_sirs_stochastic(R0=R01, delta=delta1, theta=0, N=N, sigma_env=0, rho_env=0,
                                I0=Istar1, rho=1, k=10,
                                tmax=50) %>%
  filter(time > 25)

ss2 <- simulate_sirs_stochastic(R0=R02, delta=delta2, theta=0, N=N, sigma_env=0, rho_env=0,
                                I0=Istar2, rho=1, k=10,
                                tmax=50) %>%
  filter(time > 25)

ss3 <- simulate_sirs_stochastic(R0=R03, delta=delta3, theta=0, N=N, sigma_env=0, rho_env=0,
                                I0=Istar3, rho=1, k=10,
                                tmax=50) %>%
  filter(time > 25)

g1 <- ggplot(ss1) +
  geom_line(aes(time-25, I/mean(I))) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Relative prevalence") +
  ggtitle("Low resilience") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g2 <- ggplot(ss2) +
  geom_line(aes(time, I/mean(I))) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Relative prevalence") +
  ggtitle("Intermediate resilience") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g3 <- ggplot(ss3) +
  geom_line(aes(time, I/mean(I))) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Relative prevalence") +
  ggtitle("High resilience") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g4 <- ggplot(analysis_factorial_noise) +
  geom_raster(aes(R0, 1/delta/365, fill=log10(resilience*365))) +
  annotate("text", x=R01, y=1/delta1/365, col="white", label="A", size=8, hjust=0, vjust=1) +
  annotate("text", x=R02, y=1/delta2/365, col="white", label="B", size=8, hjust=0, vjust=0) +
  annotate("text", x=R03, y=1/delta3/365, col="white", label="C", size=8, hjust=1, vjust=0) +
  scale_x_continuous("Basic reproduction number", expand=c(0, 0)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Resilience (1/years)",
                       breaks=c(-1, 0, 1),
                       labels=c(0.1, 1, 10),
                       end=0.9) +
  theme(
    legend.position = "bottom"
  )

g5 <- ggplot(analysis_factorial_noise) +
  geom_raster(aes(R0, 1/delta/365, fill=log10(amplitude))) +
  annotate("text", x=R01, y=1/delta1/365, col="white", label="A", size=8, hjust=0, vjust=1) +
  annotate("text", x=R02, y=1/delta2/365, col="white", label="B", size=8, hjust=0, vjust=0) +
  annotate("text", x=R03, y=1/delta3/365, col="white", label="C", size=8, hjust=1, vjust=0) +
  scale_x_continuous("Basic reproduction number", expand=c(0, 0)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Epidemic amplitude", option="A",
                       breaks=log10(c(0.01, 0.03, 0.1, 0.3)),
                       labels=c(0.01, 0.03, 0.1, 0.3)) +
  theme(
    legend.position = "bottom"
  )

g6 <- ggplot(analysis_factorial_noise) +
  geom_point(aes(resilience*365, amplitude)) +
  scale_x_log10("Resilience (1/years)") +
  scale_y_log10("Epidemic amplitude") +
  theme(
    panel.grid = element_blank()
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6,
                   nrow=2,
                   labels=c("A", "B", "C", "D", "E", "F"))

ggsave("figure_persistence_noise.pdf", gcomb, width=12, height=6)

g7 <- ggplot(analysis_factorial_noise) +
  geom_point(aes(period, period_obs)) +
  geom_abline(intercept=0, slope=1, lty=2) +
  scale_x_log10("Predicted periodicity of epidemic cycle (years)") +
  scale_y_log10("Observed periodicity of epidemic cycle (years)") +
  theme(
    panel.grid = element_blank()
  )

ggsave("figure_persistence_noise_period.pdf", g7, width=6, height=4)
