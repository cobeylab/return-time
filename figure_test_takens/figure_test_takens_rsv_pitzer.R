library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_rsv_pitzer.R")
source("../R/distfun.R")
source("../R/takens.R")

ss_perturb <- simulate_rsv_pitzer(tmax=2030) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_rsv_pitzer(tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010)

time <- ss_perturb$time

mat_perturb_SI <- matrix(c(ss_perturb$Seff, log(ss_perturb$prevalence)), ncol=2)
mat_unperturb_SI <- matrix(c(ss_unperturb$Seff, log(ss_unperturb$prevalence)), ncol=2)

dist_SI <- distfun(mat_perturb_SI, mat_unperturb_SI,
                ss_perturb$time,
                2020, "all")

tau <- floor(365/4)

mat_perturb_takens <- takens(log(ss_perturb$prevalence), d=2, tau=tau)
mat_unperturb_takens <- takens(log(ss_unperturb$prevalence), d=2, tau=tau)

time_takens <- tail(ss_perturb$time, -tau)

dist_takens <- distfun(mat_perturb_takens, mat_unperturb_takens,
                       tail(ss_perturb$time,-tau),
                       2020, "all")

mat_perturb_takens2 <- takens(log(ss_perturb$prevalence), d=3, tau=tau)
mat_unperturb_takens2 <- takens(log(ss_unperturb$prevalence), d=3, tau=tau)

time_takens2 <- tail(ss_perturb$time, -tau*2)

dist_takens2 <- distfun(mat_perturb_takens2, mat_unperturb_takens2,
                        tail(ss_perturb$time,-tau*2),
                        2020, "all")

g_timeseries <- ggplot(ss_perturb) +
  geom_line(aes(time, prevalence, col=group)) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

eigen_all <- eigen_rsv_pitzer()

ee <- -Re(eigen_all$values[Im(eigen_all$values)!=0])[1]

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
lfit_SI <- lm(log(dist)~time, data=filter(dist_SI, time>=maxt_SI, time<=maxt_SI+1/ee))

g2 <- ggplot(filter(dist_SI, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_SI+1/ee, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_SI, time>=maxt_SI, time<=maxt_SI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI, time>=maxt_SI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_SI, time>=maxt_SI), aes(time, dist), col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI)[1]) * exp(-ee*(x-maxt_SI)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g3 <- ggplot(bind_cols(as.data.frame(mat_perturb_takens), data.frame(time=time_takens))) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_takens), data.frame(time=time_takens)),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10("I(t)")+
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens <- dist_takens$time[which.max(dist_takens$dist)]
lfit_takens <- lm(log(dist)~time, data=filter(dist_takens, time>=maxt_takens, time<=maxt_takens+1/ee))

g4 <- ggplot(filter(dist_takens, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_takens+1/ee, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_takens, time>=maxt_takens, time<=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens, time>=maxt_takens+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens, time>=maxt_takens), aes(time, dist), col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_takens)[1]) * exp(-ee*(x-maxt_takens)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 40)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_takens2 <- dist_takens2$time[which.max(dist_takens2$dist)]
lfit_takens2 <- lm(log(dist)~time, data=filter(dist_takens2, time>=maxt_takens2, time<=maxt_takens2+1/ee))

g6 <- ggplot(filter(dist_takens2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_takens2, lty=2, col="#224B95") +
  geom_vline(xintercept=maxt_takens2+1/ee, lty=2, col="#224B95") +
  geom_smooth(data=filter(dist_takens2, time>=maxt_takens2, time<=maxt_takens2+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens2, time>=maxt_takens2+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_smooth(data=filter(dist_takens2, time>=maxt_takens2), aes(time, dist), col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_takens2)[1]) * exp(-ee*(x-maxt_takens2)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 40)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gtop <- ggarrange(g_timeseries, labels="A")
gbottom <- ggarrange(g1, g2, g3, g4, g3, g6, nrow=2, byrow = FALSE,
                     labels=c("B", "C", "D", "E", "F", "G"))

gcomb <- arrangeGrob(gtop, gbottom, nrow=2, heights=c(1, 2))

ggsave("figure_test_takens_rsv_pitzer.pdf", gcomb, width=12, height=8)
