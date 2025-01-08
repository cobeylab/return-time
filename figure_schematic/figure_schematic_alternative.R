library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_sirs_stoch.R")
source("../R/npi_random.R")
source("../R/takens.R")

npifun_random <- npifun_random_generate(duration=4,
                                        npimin=0.5,
                                        seed=11)

npidata <- data.frame(
  time=seq(2015, 2030, by=1/52/7),
  npi=sapply(seq(2015, 2030, by=1/52/7), npifun_random)
)

plot(npidata, type="l")

ss <- simulate_SIRS_stoch(
  R0=2,
  theta=0.1,
  delta=1/52/7/1,
  npifun=npifun_random)

ee <- eigen_sirs(b_1=2*(365/7+1/50),
                 mu=1/50,
                 gamma=365/7,
                 delta=1/1)

ss_weekly <- ss %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases)
  ) %>%
  mutate(
    time=year+week/52
  ) %>%
  filter(
    year > 2014, year < 2030
  )

plot(ss_weekly$cases, type="l", log="y")

logcases_pre <- log(filter(ss_weekly, year < 2020)$cases+1)

acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)

tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]

n.fnn <- fnn(logcases_pre, dmax=25, tau=tau, R_tol=10)

d <- which(n.fnn==0)[1]+1

takens_perturb <- takens(log(ss_weekly$cases+1), d=d, tau=tau)
takens_unperturb <- takens(logcases_pre, d=d, tau=tau)

takens_data <- as.data.frame(takens_perturb)
takens_data$time <- tail(ss_weekly$year+ss_weekly$week/52, -tau*(d-1))

dist <- sapply(1:nrow(takens_perturb), function(i) {
  dd <- sqrt(colSums((takens_perturb[i,] - t(takens_unperturb))^2))
  
  min(dd[dd>0])
})

plot(dist, log="y")

time <- tail(ss_weekly$year+ss_weekly$week/52, -tau*(d-1))

distdata_takens <- data.frame(
  time=time,
  dist=dist
)

g1 <- ggplot(npidata) +
  annotate("rect", xmin=2020, xmax=2024, ymin=-Inf, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept=2015:2030, lty=3, col="gray70", lwd=0.5) +
  geom_line(aes(time, npi*100)) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_continuous("% transmission reduction", limits=c(0, 105), expand=c(0, 0)) +
  ggtitle("A. Non-pharmaceutical interventions") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g2 <- ggplot(ss_weekly) +
  annotate("rect", xmin=2020, xmax=2024, ymin=-Inf, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept=2015:2030, lty=3, col="gray70", lwd=0.5) +
  # geom_line(data=ss, aes(time, S/1e8*10000+10000), col="orange", lty=2) +
  geom_line(aes(time, cases)) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_continuous("Cases", expand=c(0, 0)) +
  ggtitle("B. Observed case time series") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g3 <- ggplot(takens_data)  +
  annotate("rect", xmin=2020, xmax=2024, ymin=-Inf, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept = 2015:2030, lty=3, col="gray") +
  geom_line(aes(time, V1, col=time)) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_continuous(expression(X(t)), expand=c(0, 0)) +
  ggtitle("C. Delayed copies of logged cases") +
  scale_color_viridis_c() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

g4 <- ggplot(takens_data) +
  annotate("rect", xmin=2020+tau/52, xmax=2024+tau/52, ymin=-Inf, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept = 2015:2030, lty=3, col="gray") +
  geom_line(aes(time, V2, col=time)) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_continuous(expression(X(t-tau)), expand=c(0, 0)) +
  scale_color_viridis_c() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

g5 <- ggplot(takens_data) +
  annotate("rect", xmin=2020+2*tau/52, xmax=2024+2*tau/52, ymin=-Inf, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept = 2015:2030, lty=3, col="gray") +
  geom_line(aes(time, V3, col=time)) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_continuous(expression(X(t-2~tau)), expand=c(0, 0)) +
  scale_color_viridis_c() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

g6 <- ggplot(takens_data) +
  geom_path(aes(V1, V2, col=time)) +
  scale_x_continuous(expression(X(t))) +
  scale_y_continuous(expression(X(t-tau))) +
  scale_color_viridis_c() +
  ggtitle("D. Delayed embedding") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

g7 <- ggplot(distdata_takens) +
  annotate("rect", xmin=2020, xmax=2024, ymin=0.08, ymax=Inf, fill="gray", alpha=0.4) +
  geom_vline(xintercept = 2015:2030, lty=3, col="gray") +
  geom_point(aes(time, dist), shape=1, size=0.5) +
  # geom_smooth(data=filter(distdata_takens, time>time[dist==max(dist)]), aes(time, dist), col="#224B95", fill="#224B95",
  #             method="loess") +
  # geom_function(fun=function(x) 10*exp(-ee*(x-2022)), lty=2, col="#EF6351", lwd=0.7) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_log10("Distance from attractor", limits=c(0.08, 14), expand=c(0, 0)) +
  ggtitle("E. Nearest neighbor distance from attractor") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

plot(distdata_takens, log="y")
curve(10*exp(-ee*(x-2023)), add=TRUE)

gcomb1 <- ggarrange(g1, g2)

gcomb2 <- ggarrange(g3, g4, g5, ncol=1)

gcomb3 <- ggarrange(g6)

gcomb4 <- ggarrange(g7)

lay <- rbind(
  c(1, 2, 3),
  c(1, 2, 4)
)

gfinal <- arrangeGrob(gcomb1, gcomb2, gcomb3, gcomb4, 
                      layout_matrix = lay)

ggsave("figure_schematic_alternative.pdf", gfinal, width=12, height=4)
