library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_rsv_pitzer.R")
source("../R/distfun.R")

ss_perturb_sirs <- simulate_sirs(
  tmax=2040) %>%
  filter(time >= 2010, time < 2029.5) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_sirs <- simulate_sirs(tmax=2040,npifun=function(x)1) %>%
  filter(time >= 2010, time < 2029.5)

time <- ss_perturb_sirs$time

mat_perturb_SI_sirs <- matrix(c(ss_perturb_sirs$S, log(ss_perturb_sirs$I)), ncol=2)
mat_unperturb_SI_sirs <- matrix(c(ss_unperturb_sirs$S, log(ss_unperturb_sirs$I)), ncol=2)

dist_SI_sirs <- distfun(mat_perturb_SI_sirs, mat_unperturb_SI_sirs,
                   ss_perturb_sirs$time,
                   2020, "all")

g1 <- ggplot(ss_perturb_sirs) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("SIRS model, same attractor") +
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
  # geom_vline(xintercept=maxt_SI_sirs, lty=2, col="#224B95") +
  # geom_smooth(data=filter(dist_SI_sirs, time>=maxt_SI_sirs), aes(time, dist), method="lm", col="#224B95") +
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
  npifun=function(t) {
    if (t >= 2020 & t < 2020.5) {
      npi <- 0.5
    } else if (t >= 2020.5) {
      npi <- 0.9
    } else {
      npi <- 1
    }
    
    return(npi)
  },
  tmax=2040) %>%
  filter(time >= 2010, time < 2029.5) %>%
  mutate(
    group=time < 2020
  )

mat_perturb_SI_sirs2 <- matrix(c(ss_perturb_sirs2$S, log(ss_perturb_sirs2$I)), ncol=2)

dist_SI_sirs2 <- distfun(mat_perturb_SI_sirs2, mat_unperturb_SI_sirs,
                        ss_perturb_sirs$time,
                        2020, "all")

g4 <- ggplot(ss_perturb_sirs2) +
  annotate("rect", xmin=2020, xmax=Inf, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("SIRS model, different attractor") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

g5 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_sirs2), data.frame(time=time))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_sirs), data.frame(time=time)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_sirs2 <- eigen_sirs(b_1=2*(365/7+1/50)*0.9)

maxt_SI_sirs2 <- dist_SI_sirs2$time[which.max(dist_SI_sirs2$dist)]
lfit_SI_sirs2 <- lm(log(dist)~time, data=filter(dist_SI_sirs2, time>=maxt_SI_sirs2, time < maxt_SI_sirs2+1/ee_sirs2))

g6 <- ggplot(filter(dist_SI_sirs2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", lwd=1) +
  geom_smooth(aes(time, dist), method="loess", col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI_sirs2)[1]) * exp(-ee_sirs2*(x-maxt_SI_sirs2)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5),
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-3, 10)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ss_perturb_rsv_pitzer <- simulate_rsv_pitzer(tmax=2030) %>%
  filter(time >= 2010, time < 2029.5) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_rsv_pitzer <- simulate_rsv_pitzer(tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010, time < 2029.5)

time_rsv_pitzer <- ss_perturb_rsv_pitzer$time

mat_perturb_SI_rsv_pitzer <- matrix(c(ss_perturb_rsv_pitzer$Seff, log(ss_perturb_rsv_pitzer$prevalence)), ncol=2)
mat_unperturb_SI_rsv_pitzer <- matrix(c(ss_unperturb_rsv_pitzer$Seff, log(ss_unperturb_rsv_pitzer$prevalence)), ncol=2)

dist_SI_rsv_pitzer <- distfun(mat_perturb_SI_rsv_pitzer, mat_unperturb_SI_rsv_pitzer,
                              ss_perturb_rsv_pitzer$time,
                              2020, "all")

g7 <- ggplot(ss_perturb_rsv_pitzer) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence, col=group)) +
  geom_line(data=ss_unperturb_rsv_pitzer, aes(time, prevalence), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Stage-structured model") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

eigen_all_rsv_pitzer <- eigen_rsv_pitzer()

ee_rsv_pitzer <- -Re(eigen_all_rsv_pitzer$values[Im(eigen_all_rsv_pitzer$values)!=0])[1]

g8 <- ggplot(bind_cols(as.data.frame(mat_perturb_SI_rsv_pitzer), data.frame(time=time_rsv_pitzer))) +
  geom_path(aes(V1, exp(V2)), col="#224B95") +
  geom_path(data=bind_cols(as.data.frame(mat_unperturb_SI_rsv_pitzer), data.frame(time=time_rsv_pitzer)),
            aes(V1, exp(V2)), col="#EF6351") +
  scale_x_continuous("Susceptibles") +
  scale_y_log10("Infected") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

maxt_SI_rsv_pitzer <- dist_SI_rsv_pitzer$time[which.max(dist_SI_rsv_pitzer$dist)]
lfit_SI_rsv_pitzer <- lm(log(dist)~time, data=filter(dist_SI_rsv_pitzer, time>=maxt_SI_rsv_pitzer, time<=maxt_SI_rsv_pitzer+1/ee_rsv_pitzer))

g9 <- ggplot(filter(dist_SI_rsv_pitzer, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", lwd=1) +
  geom_smooth(aes(time, dist), method="loess", col="orange") +
  # geom_smooth(data=filter(dist_SI_rsv_pitzer, time>=maxt_SI_rsv_pitzer), 
  #             aes(time, dist), col="orange") +
  geom_function(fun=function(x) exp(predict(lfit_SI_rsv_pitzer)[1]) * exp(-ee_rsv_pitzer*(x-maxt_SI_rsv_pitzer)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5),
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-2, 20)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9,
                   byrow = FALSE,
                   labels=LETTERS[1:9])

ggsave("figure2_simple.pdf", gcomb, width=9, height=8)
