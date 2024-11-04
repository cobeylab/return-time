library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_bhattacharyya.R")
source("../R/distfun.R")
source("../R/takens.R")

ss_perturb <- simulate_bhattacharyya(tmax=2030,
                                 tmin=1905) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_bhattacharyya(tmin=1905, tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010)

time <- ss_perturb$time

mat_perturb_SI1 <- matrix(c(ss_perturb$Seff1, log(ss_perturb$prevalence1)), ncol=2)
mat_unperturb_SI1 <- matrix(c(ss_unperturb$Seff1, log(ss_unperturb$prevalence1)), ncol=2)

dist_SI1 <- distfun(mat_perturb_SI1, mat_unperturb_SI1,
                ss_perturb$time,
                2020, "all") %>%
  filter(time>=2020)

tau <- floor(365/4)

mat_perturb_takens1 <- takens(log(ss_perturb$prevalence1), d=2, tau=tau)
mat_unperturb_takens1 <- takens(log(ss_unperturb$prevalence1), d=2, tau=tau)

time_takens <- tail(ss_perturb$time, -tau)

dist_takens1 <- distfun(mat_perturb_takens1, mat_unperturb_takens1,
                        time_takens,
                   2020, "all") %>%
  filter(time>=2020)

mat_perturb_takens12 <- takens(log(ss_perturb$prevalence1), d=3, tau=tau)
mat_unperturb_takens12 <- takens(log(ss_unperturb$prevalence1), d=3, tau=tau)

time_takens2 <- tail(ss_perturb$time, -tau*2)

dist_takens12 <- distfun(mat_perturb_takens12, mat_unperturb_takens12,
                        time_takens2,
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

eigen_all <- eigen_bhattacharyya()

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
lfit_SI1 <- lm(log(dist)~time, data=filter(dist_SI1, time>=maxt_SI1))

g2 <- ggplot(filter(dist_SI1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI1, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_SI1+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI, time>=maxt_SI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI1, time>=maxt_SI1), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_SI1)[1]) * exp(-ee1*(x-maxt_SI1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI1, time>=maxt_SI1), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g3 <- ggplot(bind_cols(as.data.frame(mat_perturb_takens1), data.frame(time=time_takens))) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_takens1), data.frame(time=time_takens)),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10("I(t)")+
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens1 <- dist_takens1$time[which.max(dist_takens1$dist)]
lfit_takens1 <- lm(log(dist)~time, data=filter(dist_takens1, time>=maxt_takens1))

g4 <- ggplot(filter(dist_takens1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens1, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_takens1+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_takens, time>=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens1, time>=maxt_takens1), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_takens1)[1]) * exp(-ee1*(x-maxt_takens1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_takens1, time>=maxt_takens1), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens12 <- dist_takens12$time[which.max(dist_takens12$dist)]
lfit_takens12 <- lm(log(dist)~time, data=filter(dist_takens12, time>=maxt_takens12))

g6 <- ggplot(filter(dist_takens12, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens12, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_takens12+1/ee1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_takens, time>=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens12, time>=maxt_takens12), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_takens12)[1]) * exp(-ee1*(x-maxt_takens12)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_takens12, time>=maxt_takens12), aes(time, dist), method="lm", col="#224B95") +
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

mat_perturb_takens2 <- takens(log(ss_perturb$prevalence2), d=2, tau=tau)
mat_unperturb_takens2 <- takens(log(ss_unperturb$prevalence2), d=2, tau=tau)

dist_takens2 <- distfun(mat_perturb_takens2, mat_unperturb_takens2,
                    time_takens,
                    2020, "all") %>%
  filter(time>=2020)

mat_perturb_takens22 <- takens(log(ss_perturb$prevalence2), d=3, tau=tau)
mat_unperturb_takens22 <- takens(log(ss_unperturb$prevalence2), d=3, tau=tau)

dist_takens22 <- distfun(mat_perturb_takens22, mat_unperturb_takens22,
                        time_takens2,
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
lfit_SI2 <- lm(log(dist)~time, data=filter(dist_SI2, time>=maxt_SI2))

g8 <- ggplot(filter(dist_SI2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI2, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_SI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI2, time>=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI2, time>=maxt_SI2), aes(time, dist), col="orange", method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_SI2)[1]) * exp(-ee2*(x-maxt_SI2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI2, time>=maxt_SI2), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g9 <- ggplot(bind_cols(as.data.frame(mat_perturb_takens2), data.frame(time=time_takens))) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_takens2), data.frame(time=time_takens)),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10("I(t)")+
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens2 <- dist_takens2$time[which.max(dist_takens2$dist)]
lfit_takens2 <- lm(log(dist)~time, data=filter(dist_takens2, time>=maxt_takens2))

g10 <- ggplot(filter(dist_takens2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens2, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_takens2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_takens, time>=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens2, time>=maxt_takens2), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_takens2)[1]) * exp(-ee2*(x-maxt_takens2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_takens2, time>=maxt_takens2), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens22 <- dist_takens22$time[which.max(dist_takens22$dist)]
lfit_takens22 <- lm(log(dist)~time, data=filter(dist_takens22, time>=maxt_takens22))

g12 <- ggplot(filter(dist_takens22, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens22, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_takens22+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_takens, time>=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens22, time>=maxt_takens22), aes(time, dist), col="orange",
              method="loess") +
  geom_function(fun=function(x) exp(predict(lfit_takens22)[1]) * exp(-ee2*(x-maxt_takens22)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_takens22, time>=maxt_takens22), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gtop1 <- ggarrange(g_timeseries1, labels="A")
gbottom1 <- ggarrange(g1, g2, g3, g4, g3, g6, nrow=2, byrow = FALSE,
                     labels=c("B", "C", "D", "E", "F", "G"))

gtop2 <- ggarrange(g_timeseries2, labels="H")
gbottom2 <- ggarrange(g7, g8, g9, g10, g9, g12, nrow=2, byrow = FALSE,
                      labels=c("I", "J", "K", "L", "M", "N"))

gcomb <- arrangeGrob(gtop1, gbottom1, gtop2, gbottom2, nrow=4, heights=c(1, 2, 1, 2))

ggsave("figure_test_takens_bhattacharyya.pdf", gcomb, width=12, height=16)
