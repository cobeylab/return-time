library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_rsv_pitzer.R")
source("../R/simulate_bhattacharyya.R")
source("../R/distfun.R")
source("../R/takens.R")

ss_perturb_sirs <- simulate_sirs(tmax=2030) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_sirs <- simulate_sirs(tmax=2030,npifun=function(x)1) %>%
  filter(time >= 2010)

logI_pre_sirs <- log(filter(ss_perturb_sirs, time < 2020)$I)

acfout_sirs <- acf(logI_pre_sirs, lag.max=(365*2), plot=FALSE)

tau_sirs <- which(head(acfout_sirs$acf, -1) > 0 & tail(acfout_sirs$acf, -1) < 0)[1]

n.fnn_sirs <- fnn(logI_pre_sirs, dmax=4, tau=tau_sirs, R_tol=10)

d_sirs <- which(n.fnn_sirs==0)[1]+1

takens_perturb_sirs <- takens(log(ss_perturb_sirs$I), d=d_sirs, tau=tau_sirs)
takens_unperturb_sirs <- takens(logI_pre_sirs, d=d_sirs, tau=tau_sirs)

dist_sirs <- sapply(1:nrow(takens_perturb_sirs), function(i) {
  dd <- sqrt(colSums((takens_perturb_sirs[i,] - t(takens_unperturb_sirs))^2))
  
  min(dd[dd>0])
})

distdata_sirs <- data.frame(
  time=tail(ss_perturb_sirs$time, -tau_sirs*(d_sirs-1)),
  dist=dist_sirs
)

g1 <- ggplot(ss_perturb_sirs) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, I, col=group)) +
  geom_line(data=ss_unperturb_sirs, aes(time, I), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.13)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("SIRS model") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_sirs <- eigen_sirs()

g2 <- ggplot(as.data.frame(takens_perturb_sirs)) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=as.data.frame(takens_unperturb_sirs),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10(expression(I(t))) +
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

lfit_sirs <- lm(log(dist)~time, data=filter(distdata_sirs, time>=2022))

g3 <- ggplot(filter(distdata_sirs, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1)  +
  geom_smooth(aes(time, dist), method="loess", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_sirs)[1]) * exp(-ee_sirs*(x-2022)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 30)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ss_perturb_rsv_pitzer <- simulate_rsv_pitzer(tmax=2030) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_rsv_pitzer <- simulate_rsv_pitzer(tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010)

logI_pre_rsv_pitzer <- log(filter(ss_perturb_rsv_pitzer, time < 2020)$prevalence)

acfout_rsv_pitzer <- acf(logI_pre_rsv_pitzer, lag.max=(365*2), plot=FALSE)

tau_rsv_pitzer <- which(head(acfout_rsv_pitzer$acf, -1) > 0 & tail(acfout_rsv_pitzer$acf, -1) < 0)[1]

n.fnn_rsv_pitzer <- fnn(logI_pre_rsv_pitzer, dmax=4, tau=tau_rsv_pitzer, R_tol=10)

d_rsv_pitzer <- which(n.fnn_rsv_pitzer==0)[1]+1

takens_perturb_rsv_pitzer <- takens(log(ss_perturb_rsv_pitzer$prevalence), d=d_rsv_pitzer, tau=tau_rsv_pitzer)
takens_unperturb_rsv_pitzer <- takens(logI_pre_rsv_pitzer, d=d_rsv_pitzer, tau=tau_rsv_pitzer)

dist_rsv_pitzer <- sapply(1:nrow(takens_perturb_rsv_pitzer), function(i) {
  dd <- sqrt(colSums((takens_perturb_rsv_pitzer[i,] - t(takens_unperturb_rsv_pitzer))^2))
  
  min(dd[dd>0])
})

distdata_rsv_pitzer <- data.frame(
  time=tail(ss_perturb_rsv_pitzer$time, -tau_rsv_pitzer*(d_rsv_pitzer-1)),
  dist=dist_rsv_pitzer
)

g4 <- ggplot(ss_perturb_rsv_pitzer) +
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

g5 <- ggplot(as.data.frame(takens_perturb_rsv_pitzer)) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=as.data.frame(takens_unperturb_rsv_pitzer),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10(expression(I(t))) +
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

lfit_rsv_pitzer <- lm(log(dist)~time, data=filter(distdata_rsv_pitzer, time>=2022))

g6 <- ggplot(filter(distdata_rsv_pitzer, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1)  +
  geom_smooth(aes(time, dist), method="loess", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_rsv_pitzer)[1]) * exp(-ee_rsv_pitzer*(x-2022)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 30)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ss_perturb_bh <- simulate_bhattacharyya(tmax=2030,
                                        tmin=1905) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb_bh <- simulate_bhattacharyya(tmin=1905, tmax=2030,npifun=function(x) 1) %>%
  filter(time >= 2010)

logI_pre_bh1 <- log(filter(ss_perturb_bh, time < 2020)$prevalence1)

acfout_bh1 <- acf(logI_pre_bh1, lag.max=(365*2), plot=FALSE)

tau_bh1 <- which(head(acfout_bh1$acf, -1) > 0 & tail(acfout_bh1$acf, -1) < 0)[1]

n.fnn_bh1 <- fnn(logI_pre_bh1, dmax=4, tau=tau_bh1, R_tol=10)

d_bh1 <- which(n.fnn_bh1==0)[1]+1

takens_perturb_bh1 <- takens(log(ss_perturb_bh$prevalence1), d=d_bh1, tau=tau_bh1)
takens_unperturb_bh1 <- takens(logI_pre_bh1, d=d_bh1, tau=tau_bh1)

dist_bh1 <- sapply(1:nrow(takens_perturb_bh1), function(i) {
  dd <- sqrt(colSums((takens_perturb_bh1[i,] - t(takens_unperturb_bh1))^2))
  
  min(dd[dd>0])
})

distdata_bh1 <- data.frame(
  time=tail(ss_perturb_bh$time, -tau_bh1*(d_bh1-1)),
  dist=dist_bh1
)

g7 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence1, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence1), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, strain 1") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

eigen_all_bh <- eigen_bhattacharyya()

ee_bh1 <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

g8 <- ggplot(as.data.frame(takens_perturb_bh1)) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=as.data.frame(takens_unperturb_bh1),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10(expression(I(t))) +
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

lfit_bh1 <- lm(log(dist)~time, data=filter(distdata_bh1, time>=2022))

g9 <- ggplot(filter(distdata_bh1, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1)  +
  geom_smooth(aes(time, dist), method="loess", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_bh1)[1]) * exp(-ee_bh1*(x-2022)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 30)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

logI_pre_bh2 <- log(filter(ss_perturb_bh, time < 2020)$prevalence2)

acfout_bh2 <- acf(logI_pre_bh2, lag.max=(365*2), plot=FALSE)

tau_bh2 <- which(head(acfout_bh2$acf, -1) > 0 & tail(acfout_bh2$acf, -1) < 0)[1]

n.fnn_bh2 <- fnn(logI_pre_bh2, dmax=4, tau=tau_bh2, R_tol=10)

d_bh2 <- which(n.fnn_bh2==0)[1]+1

takens_perturb_bh2 <- takens(log(ss_perturb_bh$prevalence2), d=d_bh2, tau=tau_bh2)
takens_unperturb_bh2 <- takens(logI_pre_bh2, d=d_bh2, tau=tau_bh2)

dist_bh2 <- sapply(1:nrow(takens_perturb_bh2), function(i) {
  dd <- sqrt(colSums((takens_perturb_bh2[i,] - t(takens_unperturb_bh2))^2))
  
  min(dd[dd>0])
})

distdata_bh2 <- data.frame(
  time=tail(ss_perturb_bh$time, -tau_bh2*(d_bh2-1)),
  dist=dist_bh2
)

g10 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence2, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence2), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, strain 2") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_bh2 <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

g11 <- ggplot(as.data.frame(takens_perturb_bh2)) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=as.data.frame(takens_unperturb_bh2),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10(expression(I(t))) +
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

lfit_bh2 <- lm(log(dist)~time, data=filter(distdata_bh2, time>=2022))

g12 <- ggplot(filter(distdata_bh2, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1)  +
  geom_smooth(aes(time, dist), method="loess", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_bh2)[1]) * exp(-ee_bh2*(x-2022)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 30)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )


logI_pre_bh_both <- log(filter(ss_perturb_bh, time < 2020)$prevalence1+filter(ss_perturb_bh, time < 2020)$prevalence2)

acfout_bh_both <- acf(logI_pre_bh_both, lag.max=(365*2), plot=FALSE)

tau_bh_both <- which(head(acfout_bh_both$acf, -1) > 0 & tail(acfout_bh_both$acf, -1) < 0)[1]

n.fnn_bh_both <- fnn(logI_pre_bh_both, dmax=4, tau=tau_bh_both, R_tol=10)

d_bh_both <- which(n.fnn_bh_both==0)[1]+1

takens_perturb_bh_both <- takens(log(ss_perturb_bh$prevalence1+ss_perturb_bh$prevalence2), d=d_bh_both, tau=tau_bh_both)
takens_unperturb_bh_both <- takens(logI_pre_bh_both, d=d_bh_both, tau=tau_bh_both)

dist_bh_both <- sapply(1:nrow(takens_perturb_bh_both), function(i) {
  dd <- sqrt(colSums((takens_perturb_bh_both[i,] - t(takens_unperturb_bh_both))^2))
  
  min(dd[dd>0])
})

distdata_bh_both <- data.frame(
  time=tail(ss_perturb_bh$time, -tau_bh_both*(d_bh_both-1)),
  dist=dist_bh_both
)

g13 <- ggplot(ss_perturb_bh) +
  annotate("rect", xmin=2020, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="gray80") +
  geom_line(aes(time, prevalence2, col=group)) +
  geom_line(data=ss_unperturb_bh, aes(time, prevalence2), lty=2, col="#EF6351") +
  scale_x_continuous("Year", expand=c(0, 0), limits=c(NA, 2029.5)) +
  scale_y_continuous("Infected", expand=c(0, 0), limits=c(0, 0.2)) +
  scale_color_manual(values=c("#224B95", "#EF6351")) +
  ggtitle("Two strain model, both strains") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ee_bh_both <- -max(Re(eigen_all_bh$values[Im(eigen_all_bh$values)!=0]))

g14 <- ggplot(as.data.frame(takens_perturb_bh_both)) +
  geom_path(aes(exp(V1), exp(V2)), col="#224B95") +
  geom_path(data=as.data.frame(takens_unperturb_bh_both),
            aes(exp(V1), exp(V2)), col="#EF6351") +
  scale_x_log10(expression(I(t))) +
  scale_y_log10(expression(I(t-tau))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

lfit_bh_both <- lm(log(dist)~time, data=filter(distdata_bh_both, time>=2022))

g15 <- ggplot(filter(distdata_bh_both, time >= 2020)) +
  geom_line(aes(time, dist), col="#224B95", alpha=0.2, lwd=1)  +
  geom_smooth(aes(time, dist), method="loess", col="#224B95") +
  geom_function(fun=function(x) exp(predict(lfit_bh_both)[1]) * exp(-ee_bh_both*(x-2022)),
                lwd=0.7, lty=3) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2)) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e-5, 30)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12,
                   g13, g14, g15,
                   nrow=3,
                   byrow = FALSE,
                   labels=LETTERS[1:15])

ggsave("figure2_takens_nn.pdf", gcomb, width=15, height=8)
