library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_rsv_pitzer.R")
source("../R/distfun.R")

ss_perturb_sirs <- simulate_sirs(
  theta=0,
  tmax=2040) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_sirs <- simulate_sirs(theta=0,
                                   tmax=2040,npifun=function(x)1) %>%
  filter(time >= 2010)

time <- ss_perturb_sirs$time

mat_perturb_SI_sirs <- matrix(c(ss_perturb_sirs$S, log(ss_perturb_sirs$I)), ncol=2)
mat_unperturb_SI_sirs <- matrix(c(ss_unperturb_sirs$S, log(ss_unperturb_sirs$I)), ncol=2)

dist_SI_sirs <- distfun(mat_perturb_SI_sirs, mat_unperturb_SI_sirs,
                   ss_perturb_sirs$time,
                   2020, "all",
                   scale=FALSE)

g1 <- ggplot(ss_perturb_sirs) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("No seasonal transmission") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_sirs <- eigen_sirs()

g2 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_sirs), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_sirs), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_sirs <- dist_SI_sirs$time[which.max(dist_SI_sirs$dist)]
lfit_SI_sirs <- lm(log(dist)~time, data=filter(dist_SI_sirs, time>=maxt_SI_sirs))

g3 <- ggplot(filter(dist_SI_sirs, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", lwd=1) +
  geom_smooth(aes(time, dist), method="loess", col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI_sirs)[1]) * exp(-ee_sirs*(x-maxt_SI_sirs)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5),
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-3, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ss_perturb_sirs2 <- simulate_sirs(
  theta=0.2,
  tmax=2040) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_sirs2 <- simulate_sirs(theta=0.2,
                                   tmax=2040,npifun=function(x)1) %>%
  filter(time >= 2010)

mat_perturb_SI_sirs2 <- matrix(c(ss_perturb_sirs2$S, log(ss_perturb_sirs2$I)), ncol=2)
mat_unperturb_SI_sirs2 <- matrix(c(ss_unperturb_sirs2$S, log(ss_unperturb_sirs2$I)), ncol=2)

dist_SI_sirs2 <- distfun(mat_perturb_SI_sirs2, mat_unperturb_SI_sirs2,
                        ss_perturb_sirs2$time,
                        2020, "all")

g4 <- ggplot(ss_perturb_sirs2) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs2, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Moderately seasonal transmission") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g5 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_sirs2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_sirs2), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_sirs2 <- dist_SI_sirs2$time[which.max(dist_SI_sirs2$dist)]
lfit_SI_sirs2 <- lm(log(dist)~time, data=filter(dist_SI_sirs2, time>=maxt_SI_sirs2))

g6 <- ggplot(filter(dist_SI_sirs2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", lwd=1) +
  geom_smooth(aes(time, dist), method="loess", col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI_sirs2)[1]) * exp(-ee_sirs*(x-maxt_SI_sirs2)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5),
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-3, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )


ss_perturb_sirs3 <- simulate_sirs(
  theta=0.4,
  tmax=2040) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_sirs3 <- simulate_sirs(theta=0.4,
                                    tmax=2040,npifun=function(x)1) %>%
  filter(time >= 2010)

mat_perturb_SI_sirs3 <- matrix(c(ss_perturb_sirs3$S, log(ss_perturb_sirs3$I)), ncol=2)
mat_unperturb_SI_sirs3 <- matrix(c(ss_unperturb_sirs3$S, log(ss_unperturb_sirs3$I)), ncol=2)

dist_SI_sirs3 <- distfun(mat_perturb_SI_sirs3, mat_unperturb_SI_sirs3,
                         ss_perturb_sirs3$time,
                         2020, "all")

g7 <- ggplot(ss_perturb_sirs3) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs3, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.18)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Strongly seasonal transmission") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g8 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_sirs3), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_sirs3), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_sirs3 <- dist_SI_sirs3$time[which.max(dist_SI_sirs3$dist)]
lfit_SI_sirs3 <- lm(log(dist)~time, data=filter(dist_SI_sirs3, time>=maxt_SI_sirs3))

g9 <- ggplot(filter(dist_SI_sirs3, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", lwd=1) +
  geom_smooth(aes(time, dist), method="loess", col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI_sirs3)[1]) * exp(-ee_sirs*(x-maxt_SI_sirs3)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5),
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-3, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9,
          nrow=3,
          byrow = FALSE,
          labels=LETTERS[1:9])

ggsave("figure2_simple_seas.pdf", gcomb, width=10, height=8)
