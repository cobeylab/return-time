library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/distfun.R")

ss_perturb <- simulate_sirs(tmax=2030) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_sirs(tmax=2030,npifun=function(x)1) %>%
  filter(time >= 2010)

time <- ss_perturb$time

mat_perturb_SI <- matrix(c(ss_perturb$S, log(ss_perturb$I)), ncol=2)
mat_unperturb_SI <- matrix(c(ss_unperturb$S, log(ss_unperturb$I)), ncol=2)

dist_SI <- distfun(mat_perturb_SI, mat_unperturb_SI,
                ss_perturb$time,
                2020, "all")

mat_perturb_RI <- matrix(c(ss_perturb$Reff, log(ss_perturb$I)), ncol=2)
mat_unperturb_RI <- matrix(c(ss_unperturb$Reff, log(ss_unperturb$I)), ncol=2)

dist_RI <- distfun(mat_perturb_RI, mat_unperturb_RI,
                   ss_perturb$time,
                   2020, "all")

mat_perturb_rI <- matrix(c(ss_perturb$r, log(ss_perturb$I)), ncol=2)
mat_unperturb_rI <- matrix(c(ss_unperturb$r, log(ss_unperturb$I)), ncol=2)

dist_rI <- distfun(mat_perturb_rI, mat_unperturb_rI,
                   ss_perturb$time,
                   2020, "all")

g_timeseries <- ggplot(ss_perturb) +
  geom_line(aes(time, I, col=group)) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee <- eigen_sirs()

g1 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI <- dist_SI$time[which.max(dist_SI$dist)]
lfit_SI <- lm(log(dist)~time, data=filter(dist_SI, time>=maxt_SI))

g2 <- ggplot(filter(dist_SI, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_SI, time>=maxt_SI), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_SI)[1]) * exp(-ee*(x-maxt_SI)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g3 <- ggplot(bind_cols(as.data.frame(mat_perturb_RI), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_RI), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Effective reproduction number") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_RI <- dist_RI$time[which.max(dist_RI$dist)]
lfit_RI <- lm(log(dist)~time, data=filter(dist_RI, time>=maxt_RI))

g4 <- ggplot(filter(dist_RI, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_RI, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_RI, time>=maxt_RI), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_RI)[1]) * exp(-ee*(x-maxt_RI)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g5 <- ggplot(bind_cols(as.data.frame(mat_perturb_rI), data.frame(time=time))) +
  geom_path(aes(V1/365, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_rI), data.frame(time=time)),
            aes(V1/365, exp(V2)), col="#EF6351") +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_rI <- dist_rI$time[which.max(dist_rI$dist)]
lfit_rI <- lm(log(dist)~time, data=filter(dist_rI, time>=maxt_rI))

g6 <- ggplot(filter(dist_rI, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_rI, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_rI, time>=maxt_rI), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_rI)[1]) * exp(-ee*(x-maxt_rI)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gtop <- ggarrange(g_timeseries, labels="A")
gbottom <- ggarrange(g1, g2, g3, g4, g5, g6, nrow=2, byrow = FALSE,
                     labels=c("B", "C", "D", "E", "F", "G"))

gcomb <- arrangeGrob(gtop, gbottom, nrow=2, heights=c(1, 2))

ggsave("figure_test_sirs.pdf", gcomb, width=12, height=8)
