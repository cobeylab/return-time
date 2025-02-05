library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_bhattacharyya.R")
source("../R/distfun.R")

ss_perturb_bh <- simulate_bhattacharyya(b01=2*52,
                                        b02=4*52,
                                        theta1=0.1,
                                        theta2=0.1,
                                        phi1=0,
                                        phi2=0,
                                        epsilon12=0.9,
                                        epsilon21=0.5,
                                        gamma1=52,
                                        gamma2=52,
                                        rho1=1,
                                        rho2=1,
                                        q1=1,
                                        q2=1,
                                        mu=1/70,
                                        tmin=1900,
                                        tmax=2030,
                                        npifun=npifun_default) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_bh <- simulate_bhattacharyya(b01=2*52,
                                          b02=4*52,
                                          theta1=0.1,
                                          theta2=0.1,
                                          phi1=0,
                                          phi2=0,
                                          epsilon12=0.9,
                                          epsilon21=0.5,
                                          gamma1=52,
                                          gamma2=52,
                                          rho1=1,
                                          rho2=1,
                                          q1=1,
                                          q2=1,
                                          mu=1/70,
                                          tmin=1900,
                                          tmax=2030,
                                          npifun=function(x) 1) %>%
  filter(time >= 2010)

time <- ss_perturb_bh$time

mat_perturb_SI_bh1 <- matrix(c(ss_perturb_bh$Seff1, log(ss_perturb_bh$prevalence1)), ncol=2)
mat_unperturb_SI_bh1 <- matrix(c(ss_unperturb_bh$Seff1, log(ss_unperturb_bh$prevalence1)), ncol=2)

dist_SI_bh1 <- distfun(mat_perturb_SI_bh1, mat_unperturb_SI_bh1,
                       ss_perturb_bh$time,
                       2020, "all") %>%
  filter(time>=2020)

g1 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence1, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence1), col="#EF6351", lty=2) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, strain 1") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

eigen_all_bh <- eigen_bhattacharyya(b01=2*52,
                                    b02=4*52,
                                    theta1=0,
                                    theta2=0,
                                    phi1=0,
                                    phi2=0,
                                    epsilon12=0.9,
                                    epsilon21=0.5,
                                    gamma1=52,
                                    gamma2=52,
                                    rho1=1,
                                    rho2=1,
                                    q1=1,
                                    q2=1,
                                    mu=1/70)

ee_bh1 <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

g2 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_bh1), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_bh1), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_bh1 <- dist_SI_bh1$time[which.max(dist_SI_bh1$dist)]
lfit_SI_bh1 <- lm(log(dist)~time, data=filter(dist_SI_bh1, time>=maxt_SI_bh1))

g3 <- ggplot(filter(dist_SI_bh1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI_bh1, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI, time>=maxt_SI+1/ee), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_SI_bh1)[1]) * exp(-ee_bh1*(x-maxt_SI_bh1)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI_bh1, time>=maxt_SI_bh1), aes(time, dist), method="loess", col="#224B95") +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

mat_perturb_SI_bh2 <- matrix(c(ss_perturb_bh$Seff2, log(ss_perturb_bh$prevalence2)), ncol=2)
mat_unperturb_SI_bh2 <- matrix(c(ss_unperturb_bh$Seff2, log(ss_unperturb_bh$prevalence2)), ncol=2)

dist_SI_bh2 <- distfun(mat_perturb_SI_bh2, mat_unperturb_SI_bh2,
                       ss_perturb_bh$time,
                       2020, "all") %>%
  filter(time>=2020)

g4 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence2, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence2), col="#EF6351", lty=2) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, strain 2") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_bh2 <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

g5 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_bh2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_bh2), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_bh2 <- dist_SI_bh2$time[which.max(dist_SI_bh2$dist)]
lfit_SI_bh2 <- lm(log(dist)~time, data=filter(dist_SI_bh2, time>=maxt_SI_bh2))

g6 <- ggplot(filter(dist_SI_bh2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI_bh2, lty=2, col="#224B95") +
  # geom_vline(xintercept=maxt_SI2+1/ee2, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI2, time>=maxt_SI2+1/ee2), aes(time, dist), method="lm", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_SI_bh2)[1]) * exp(-ee_bh2*(x-maxt_SI_bh2)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI_bh1, time>=maxt_SI_bh1), aes(time, dist), method="loess", col="#224B95") +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

mat_perturb_SI_bh3 <- matrix(c(ss_perturb_bh$Sall, log(ss_perturb_bh$prevalence1+ss_perturb_bh$prevalence2)), ncol=2)
mat_unperturb_SI_bh3 <- matrix(c(ss_unperturb_bh$Sall, log(ss_unperturb_bh$prevalence1+ss_unperturb_bh$prevalence2)), ncol=2)

dist_SI_bh3 <- distfun(mat_perturb_SI_bh3, mat_unperturb_SI_bh3,
                       ss_perturb_bh$time,
                       2020, "all") %>%
  filter(time>=2020)

g7 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence1+prevalence2, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence1+prevalence2), col="#EF6351", lty=2) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, both strains") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g8 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_bh3), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_bh3), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_bh3 <- dist_SI_bh3$time[which.max(dist_SI_bh3$dist)]
lfit_SI_bh3 <- lm(log(dist)~time, data=filter(dist_SI_bh3, time>=maxt_SI_bh3))

g9 <- ggplot(filter(dist_SI_bh3, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_vline(xintercept=maxt_SI_bh3, lty=2, col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_SI_bh3)[1]) * exp(-ee_bh1*(x-maxt_SI_bh3)),
                lwd=0.7, lty=3) +
  geom_smooth(data=filter(dist_SI_bh3, time>=maxt_SI_bh3), aes(time, dist), method="loess", col="#224B95") +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9,
                   nrow=3,
                   byrow = FALSE,
                   labels=LETTERS[1:9])

ggsave("figure2_multi.pdf", gcomb, width=15, height=8)
