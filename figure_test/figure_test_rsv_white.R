library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_rsv_white.R")
source("../R/distfun.R")

ss_perturb <- simulate_rsv_white(tmax=2050,
                                 tmin=1905) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_rsv_white(tmin=1905, tmax=2053,npifun=function(x) 1) %>%
  filter(time >= 2010) %>%
  tail(-365*3)

time <- ss_perturb$time

mat_perturb_SI1 <- matrix(c(ss_perturb$Seff1, log(ss_perturb$prevalence1)), ncol=2)
mat_unperturb_SI1 <- matrix(c(ss_unperturb$Seff1, log(ss_unperturb$prevalence1)), ncol=2)

dist_SI1 <- distfun(mat_perturb_SI1, mat_unperturb_SI1,
                ss_perturb$time,
                2020, "all") %>%
  filter(time>=2020)

mat_perturb_RI1 <- matrix(c(ss_perturb$Reff1, log(ss_perturb$prevalence1)), ncol=2)
mat_unperturb_RI1 <- matrix(c(ss_unperturb$Reff1, log(ss_unperturb$prevalence1)), ncol=2)

dist_RI1 <- distfun(mat_perturb_RI1, mat_unperturb_RI1,
                   ss_perturb$time,
                   2020, "all") %>%
  filter(time>=2020)

mat_perturb_rI1 <- matrix(c(ss_perturb$r1, log(ss_perturb$prevalence1)), ncol=2)
mat_unperturb_rI1 <- matrix(c(ss_unperturb$r1, log(ss_unperturb$prevalence1)), ncol=2)

dist_rI1 <- distfun(mat_perturb_rI1, mat_unperturb_rI1,
                   ss_perturb$time,
                   2020, "all") %>%
  filter(time>=2020)

g_timeseries1 <- ggplot(ss_perturb) +
  geom_line(aes(time, prevalence1, col=group)) +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

eigen_all <- eigen_rsv_white()

ee1 <- -max(Re(eigen_all$values[Im(eigen_all$values)!=0]))

g1 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI1), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI1), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI1 <- dist_SI1$time[which.max(dist_SI1$dist)]
lfit_SI1 <- lm(log(dist)~time, data=filter(dist_SI1, time>=maxt_SI1, time<=maxt_SI1+1/ee1))

g2 <- ggplot(filter(dist_SI1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI1, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_SI1+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI, time>=maxt_SI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI1, time>=maxt_SI1), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_SI1)[1]) * exp(-ee1*(x-maxt_SI1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI1, time>=maxt_SI1, time<=maxt_SI1+1/ee1), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g3 <- ggplot(bind_cols(as.data.frame(mat_perturb_RI1), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_RI1), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Effective reproduction number") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_RI1 <- dist_RI1$time[which.max(dist_RI1$dist)]
lfit_RI1 <- lm(log(dist)~time, data=filter(dist_RI1, time>=maxt_RI1, time<=maxt_RI1+1/ee1))

g4 <- ggplot(filter(dist_RI1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_RI1, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_RI1+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_RI, time>=maxt_RI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_RI1, time>=maxt_RI1), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_RI1)[1]) * exp(-ee1*(x-maxt_RI1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_RI1, time>=maxt_RI1, time<=maxt_RI1+1/ee1), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g5 <- ggplot(bind_cols(as.data.frame(mat_perturb_rI1), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_rI1), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_rI1 <- dist_rI1$time[which.max(dist_rI1$dist)]
lfit_rI1 <- lm(log(dist)~time, data=filter(dist_rI1, time>=maxt_rI1, time<=maxt_rI1+1/ee1))

g6 <- ggplot(filter(dist_rI1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_rI1, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_rI1+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_rI, time>=maxt_rI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_rI1, time>=maxt_rI1), aes(time, dist), col="orange", method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_rI1)[1]) * exp(-ee1*(x-maxt_rI1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_rI1, time>=maxt_rI1, time<=maxt_rI1+1/ee1), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

mat_perturb_SI2 <- matrix(c(ss_perturb$Seff2, log(ss_perturb$prevalence2)), ncol=2)
mat_unperturb_SI2 <- matrix(c(ss_unperturb$Seff2, log(ss_unperturb$prevalence2)), ncol=2)

dist_SI2 <- distfun(mat_perturb_SI2, mat_unperturb_SI2,
                    ss_perturb$time,
                    2020, "all") %>%
  filter(time>=2020)

mat_perturb_RI2 <- matrix(c(ss_perturb$Reff2, log(ss_perturb$prevalence2)), ncol=2)
mat_unperturb_RI2 <- matrix(c(ss_unperturb$Reff2, log(ss_unperturb$prevalence2)), ncol=2)

dist_RI2 <- distfun(mat_perturb_RI2, mat_unperturb_RI2,
                    ss_perturb$time,
                    2020, "all") %>%
  filter(time>=2020)

mat_perturb_rI2 <- matrix(c(ss_perturb$r2, log(ss_perturb$prevalence2)), ncol=2)
mat_unperturb_rI2 <- matrix(c(ss_unperturb$r2, log(ss_unperturb$prevalence2)), ncol=2)

dist_rI2 <- distfun(mat_perturb_rI2, mat_unperturb_rI2,
                    ss_perturb$time,
                    2020, "all") %>%
  filter(time>=2020)

g_timeseries2 <- ggplot(ss_perturb) +
  geom_line(aes(time, prevalence2, col=group)) +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee2 <- -max(Re(eigen_all$values[Im(eigen_all$values)!=0]))

g7 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI2), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI2 <- dist_SI2$time[which.max(dist_SI2$dist)]
lfit_SI2 <- lm(log(dist)~time, data=filter(dist_SI2, time>=maxt_SI2, time<=maxt_SI2+1/ee2))

g8 <- ggplot(filter(dist_SI2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI2, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_SI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI2, time>=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI2, time>=maxt_SI2), aes(time, dist), col="orange", method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_SI2)[1]) * exp(-ee2*(x-maxt_SI2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI2, time>=maxt_SI2, time<=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g9 <- ggplot(bind_cols(as.data.frame(mat_perturb_RI2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_RI2), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Effective reproduction number") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_RI2 <- dist_RI2$time[which.max(dist_RI2$dist)]
lfit_RI2 <- lm(log(dist)~time, data=filter(dist_RI2, time>=maxt_RI2, time<=maxt_RI2+1/ee2))

g10 <- ggplot(filter(dist_RI2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_RI2, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_RI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_RI, time>=maxt_RI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_RI2, time>=maxt_RI2), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_RI2)[1]) * exp(-ee2*(x-maxt_RI2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_RI2, time>=maxt_RI2, time<=maxt_RI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g11 <- ggplot(bind_cols(as.data.frame(mat_perturb_rI2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_rI1), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Epidemic growth rate (1/day)") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_rI2 <- dist_rI2$time[which.max(dist_rI2$dist)]
lfit_rI2 <- lm(log(dist)~time, data=filter(dist_rI2, time>=maxt_rI2, time<=maxt_rI2+1/ee2))

g12 <- ggplot(filter(dist_rI2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_rI2, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_rI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_rI, time>=maxt_rI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_rI2, time>=maxt_rI2), aes(time, dist), col="orange", method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_rI2)[1]) * exp(-ee2*(x-maxt_rI2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_rI2, time>=maxt_rI2, time<=maxt_rI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gtop1 <- ggarrange(g_timeseries1, labels="A")
gbottom1 <- ggarrange(g1, g2, g3, g4, g5, g6, nrow=2, byrow = FALSE,
                     labels=c("B", "C", "D", "E", "F", "G"))

gtop2 <- ggarrange(g_timeseries2, labels="H")
gbottom2 <- ggarrange(g7, g8, g9, g10, g11, g12, nrow=2, byrow = FALSE,
                      labels=c("I", "J", "K", "L", "M", "N"))

gcomb <- arrangeGrob(gtop1, gbottom1, gtop2, gbottom2, nrow=4, heights=c(1, 2, 1, 2))

ggsave("figure_test_rsv_white.pdf", gcomb, width=12, height=16)
