library(rstan)
library(posterior)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
source("../R/eigen_discrete.R")
source("../R/simulate_sirs.R")
source("../R/simulate_bhattacharyya.R")
source("../R/simulate_sirs_stoch.R")

load("../data_simulation/data_simulation_sirs.rda")
load("../data_simulation/data_simulation_bhattacharyya.rda")

load("../stanfit_simulation/stanfit_sirs_sirs.rda")
load("../stanfit_simulation/stanfit_sirs_bhattacharyya_both.rda")
load("../stanfit_simulation/stanfit_sirs_bhattacharyya1.rda")
load("../stanfit_simulation/stanfit_sirs_bhattacharyya2.rda")

ss_sirs <- summary(stanfit_sirs_sirs)
ss_b1 <- summary(stanfit_sirs_bhattacharyya1)
ss_b2 <- summary(stanfit_sirs_bhattacharyya2)
ss_bboth <- summary(stanfit_sirs_bhattacharyya_both)

time <- data_simulation_SIRS$year+data_simulation_SIRS$week/52
time <- time[time > 2014]

npidata <- data.frame(
  time=time,
  npi=sapply(time, npifun_default)
)

fitdata_sirs <- data.frame(
  time=time,
  est=ss_sirs$summary[grepl("C\\[", rownames(ss_sirs$summary)),6],
  lwr=ss_sirs$summary[grepl("C\\[", rownames(ss_sirs$summary)),4],
  upr=ss_sirs$summary[grepl("C\\[", rownames(ss_sirs$summary)),8]
)

fitdata_b1 <- data.frame(
  time=time,
  est=ss_b1$summary[grepl("C\\[", rownames(ss_b1$summary)),6],
  lwr=ss_b1$summary[grepl("C\\[", rownames(ss_b1$summary)),4],
  upr=ss_b1$summary[grepl("C\\[", rownames(ss_b1$summary)),8]
)

fitdata_b2 <- data.frame(
  time=time,
  est=ss_b2$summary[grepl("C\\[", rownames(ss_b2$summary)),6],
  lwr=ss_b2$summary[grepl("C\\[", rownames(ss_b2$summary)),4],
  upr=ss_b2$summary[grepl("C\\[", rownames(ss_b2$summary)),8]
)

fitdata_bboth <- data.frame(
  time=time,
  est=ss_bboth$summary[grepl("C\\[", rownames(ss_bboth$summary)),6],
  lwr=ss_bboth$summary[grepl("C\\[", rownames(ss_bboth$summary)),4],
  upr=ss_bboth$summary[grepl("C\\[", rownames(ss_bboth$summary)),8]
)

npidata_sirs <- data.frame(
  time=time,
  est=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),6],
  lwr=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),4],
  upr=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),8]
)

npidata_b1 <- data.frame(
  time=time,
  est=ss_b1$summary[grepl("npieff\\[", rownames(ss_b1$summary)),6],
  lwr=ss_b1$summary[grepl("npieff\\[", rownames(ss_b1$summary)),4],
  upr=ss_b1$summary[grepl("npieff\\[", rownames(ss_b1$summary)),8]
)

npidata_b2 <- data.frame(
  time=time,
  est=ss_b2$summary[grepl("npieff\\[", rownames(ss_b2$summary)),6],
  lwr=ss_b2$summary[grepl("npieff\\[", rownames(ss_b2$summary)),4],
  upr=ss_b2$summary[grepl("npieff\\[", rownames(ss_b2$summary)),8]
)

npidata_bboth <- data.frame(
  time=time,
  est=ss_bboth$summary[grepl("npieff\\[", rownames(ss_bboth$summary)),6],
  lwr=ss_bboth$summary[grepl("npieff\\[", rownames(ss_bboth$summary)),4],
  upr=ss_bboth$summary[grepl("npieff\\[", rownames(ss_bboth$summary)),8]
)

dd_sirs <- as_draws_matrix(stanfit_sirs_sirs)
dd_b1 <- as_draws_matrix(stanfit_sirs_bhattacharyya1)
dd_b2 <- as_draws_matrix(stanfit_sirs_bhattacharyya2)
dd_bboth <- as_draws_matrix(stanfit_sirs_bhattacharyya_both)

nsamp <- 400

set.seed(101)
eigen_sirs_stan <- sapply(sample(nrow(dd_sirs), nsamp), function(j) {
  eigen_discrete_sirs(
    beta=mean(dd_sirs[j,grepl("beta", colnames(dd_sirs))]),
    omega=c(dd_sirs[j,colnames(dd_sirs)=="omega"]),
    N=1e8,
    mu=1/50/52,
    gamma=1,
    delta=c(dd_sirs[j,colnames(dd_sirs)=="delta"])
  )
})

resilience_sirs <- data.frame(
  value=-log(eigen_sirs_stan)*52
)

eigen_b1_stan <- sapply(sample(nrow(dd_b1), nsamp), function(j) {
  eigen_discrete_sirs(
    beta=mean(dd_b1[j,grepl("beta", colnames(dd_b1))]),
    omega=c(dd_b1[j,colnames(dd_b1)=="omega"]),
    N=1e8,
    mu=1/70/52,
    gamma=1,
    delta=c(dd_b1[j,colnames(dd_b1)=="delta"])
  )
})

resilience_b1 <- data.frame(
  value=-log(eigen_b1_stan)*52
)

eigen_b2_stan <- sapply(sample(nrow(dd_b2), nsamp), function(j) {
  eigen_discrete_sirs(
    beta=mean(dd_b2[j,grepl("beta", colnames(dd_b2))]),
    omega=c(dd_b2[j,colnames(dd_b2)=="omega"]),
    N=1e8,
    mu=1/70/52,
    gamma=1,
    delta=c(dd_b2[j,colnames(dd_b2)=="delta"])
  )
})

resilience_b2 <- data.frame(
  value=-log(eigen_b2_stan)*52
)

eigen_bboth_stan <- sapply(sample(nrow(dd_bboth), nsamp), function(j) {
  eigen_discrete_sirs(
    beta=mean(dd_bboth[j,grepl("beta", colnames(dd_bboth))]),
    omega=c(dd_bboth[j,colnames(dd_bboth)=="omega"]),
    N=1e8,
    mu=1/70/52,
    gamma=1,
    delta=c(dd_bboth[j,colnames(dd_bboth)=="delta"])
  )
})

resilience_bboth <- data.frame(
  value=-log(eigen_bboth_stan)*52
)

ee_sirs_true <- eigen_sirs()
ee_b_true <- eigen_bhattacharyya()

ee_b <- -max(Re(ee_b_true$values[Im(ee_b_true$values)!=0]))

g1 <- ggplot(data_simulation_SIRS) +
  geom_point(aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(data=fitdata_sirs, aes(time, est)) +
  geom_ribbon(data=fitdata_sirs, aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", limits=c(2014, 2024), expand=c(0, 0),
                     breaks=2014:2023) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  ggtitle("SIRS model") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g2 <- ggplot(npidata_sirs) +
  geom_hline(yintercept=1, lty=2) +
  geom_line(aes(time, est)) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=npidata, aes(time, npi), col="#EF6351") +
  scale_x_continuous("Year", limits=c(2020, 2024), expand=c(0, 0),
                     breaks=c(2020:2023)) +
  scale_y_continuous("Relative changes\nin transmission", limits=c(0, 2), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g3 <- ggplot(resilience_sirs) +
  geom_density(aes(value)) +
  geom_vline(xintercept=ee_sirs_true, lty=2)+
  scale_x_continuous("Intrinsic resilience (1/years)", limits=c(0.5, 2.2)) +
  scale_y_continuous("Density")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g4 <- ggplot(data_simulation_bhattacharyya1) +
  geom_point(aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(data=fitdata_b1, aes(time, est)) +
  geom_ribbon(data=fitdata_b1, aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", limits=c(2014, 2024), expand=c(0, 0),
                     breaks=2014:2023) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  ggtitle("Two strain model, strain 1") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g5 <- ggplot(npidata_b1) +
  geom_hline(yintercept=1, lty=2) +
  geom_line(aes(time, est)) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=npidata, aes(time, npi), col="#EF6351") +
  scale_x_continuous("Year", limits=c(2020, 2024), expand=c(0, 0),
                     breaks=c(2020:2023)) +
  scale_y_continuous("Relative changes\nin transmission", limits=c(0, 2), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g6 <- ggplot(resilience_b1) +
  geom_density(aes(value)) +
  geom_vline(xintercept=ee_b, lty=2)+
  scale_x_continuous("Intrinsic resilience (1/years)", limits=c(0.5, 2.2)) +
  scale_y_continuous("Density")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g7 <- ggplot(data_simulation_bhattacharyya2) +
  geom_point(aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(data=fitdata_b2, aes(time, est)) +
  geom_ribbon(data=fitdata_b2, aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", limits=c(2014, 2024), expand=c(0, 0),
                     breaks=2014:2023) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  ggtitle("Two strain model, strain 2") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g8 <- ggplot(npidata_b2) +
  geom_hline(yintercept=1, lty=2) +
  geom_line(aes(time, est)) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=npidata, aes(time, npi), col="#EF6351") +
  scale_x_continuous("Year", limits=c(2020, 2024), expand=c(0, 0),
                     breaks=c(2020:2023)) +
  scale_y_continuous("Relative changes\nin transmission", limits=c(0, 2), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g9 <- ggplot(resilience_b2) +
  geom_density(aes(value)) +
  geom_vline(xintercept=ee_b, lty=2)+
  scale_x_continuous("Intrinsic resilience (1/years)", limits=c(0.5, 2.2)) +
  scale_y_continuous("Density")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g10 <- ggplot(data_simulation_bhattacharyya_both) +
  geom_point(aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(data=fitdata_bboth, aes(time, est)) +
  geom_ribbon(data=fitdata_bboth, aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", limits=c(2014, 2024), expand=c(0, 0),
                     breaks=2014:2023) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  ggtitle("Two strain model, both strains") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g11 <- ggplot(npidata_bboth) +
  geom_hline(yintercept=1, lty=2) +
  geom_line(aes(time, est)) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.4) +
  geom_line(data=npidata, aes(time, npi), col="#EF6351") +
  scale_x_continuous("Year", limits=c(2020, 2024), expand=c(0, 0),
                     breaks=2020:2023) +
  scale_y_continuous("Relative changes\nin transmission", limits=c(0, 2), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g12 <- ggplot(resilience_bboth) +
  geom_density(aes(value)) +
  geom_vline(xintercept=ee_b, lty=2)+
  scale_x_continuous("Intrinsic resilience (1/years)", limits=c(0.5, 2.2)) +
  scale_y_continuous("Density")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12,
          byrow = FALSE,
          ncol=4,
          labels=LETTERS[1:12])

ggsave("figure4.pdf", gcomb, width=16, height=8)
