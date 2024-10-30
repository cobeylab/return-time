library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_bhattacharyya.R")
source("../R/distfun.R")

ss_perturb <- simulate_bhattacharyya(tmax=2030,
                                 tmin=1905) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_bhattacharyya(tmin=1905, tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010)

time <- ss_perturb$time

mat_perturb_all <- as.matrix(ss_perturb[,2:9])
mat_unperturb_all <- as.matrix(ss_unperturb[,2:9])

dist_all <- distfun(mat_perturb_all, mat_unperturb_all,
                   ss_perturb$time,
                   2020, "all") %>%
  filter(time>=2020)

mat_perturb_SI <- matrix(c(ss_perturb$Seff1, log(ss_perturb$prevalence1),
                           ss_perturb$Seff2, log(ss_perturb$prevalence2)), ncol=4)

mat_unperturb_SI <- matrix(c(ss_unperturb$Seff1, log(ss_unperturb$prevalence1),
                             ss_unperturb$Seff2, log(ss_unperturb$prevalence2)), ncol=4)

dist_SI <- distfun(mat_perturb_SI, mat_unperturb_SI,
                ss_perturb$time,
                2020, "all") %>%
  filter(time>=2020)

eigen_all_bh <- eigen_bhattacharyya()

ee <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

maxt_all <- dist_all$time[which.max(dist_all$dist)]
lfit_all <- lm(log(dist)~time, data=filter(dist_all, time>=maxt_all))

g1 <- ggplot(filter(dist_all, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_all, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_SI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI2, time>=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_all)[1]) * exp(-ee*(x-maxt_all)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_all, time>=maxt_all), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 20)) +
  ggtitle("All 8 dimensions") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI <- dist_SI$time[which.max(dist_SI$dist)]
lfit_SI <- lm(log(dist)~time, data=filter(dist_SI, time>=maxt_SI))

g2 <- ggplot(filter(dist_SI, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_SI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI2, time>=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_SI)[1]) * exp(-ee*(x-maxt_SI)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI, time>=maxt_SI), aes(time, dist), method="lm", col="#224B95") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  ggtitle("4 dimensions, both susceptibles and infected") +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, nrow=1)

ggsave("figure_test_bhattacharyya_4d.pdf", gcomb, width=8, height=4)
