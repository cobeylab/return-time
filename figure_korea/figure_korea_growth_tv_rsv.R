library(posterior)
library(rstan)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_korea.R")
source("../R/distfun.R")

load("../stanfit_korea/stanfit_korea_growth_tv_rsv.rda")

data_fit <- data_korea_scaled %>%
  filter(key=="RSV")

dd <- as_draws(stanfit_korea_growth_tv_rsv)

ss <- summary(stanfit_korea_growth_tv_rsv)

fitdata <- data.frame(
  year=data_fit$year,
  week=data_fit$week,
  est=ss$summary[grepl("C\\[", rownames(ss$summary)),6],
  lwr=ss$summary[grepl("C\\[", rownames(ss$summary)),4],
  upr=ss$summary[grepl("C\\[", rownames(ss$summary)),8]
)

rdata <- data.frame(
  year=data_fit$year,
  week=data_fit$week,
  est=c(ss$summary[grepl("r\\[", rownames(ss$summary)),6], NA),
  lwr=c(ss$summary[grepl("r\\[", rownames(ss$summary)),4], NA),
  upr=c(ss$summary[grepl("r\\[", rownames(ss$summary)),8], NA)
)

phasedata <- data.frame(
  year=data_fit$year,
  week=data_fit$week,
  c=ss$summary[grepl("C\\[", rownames(ss$summary)),6],
  r=c(ss$summary[grepl("r\\[", rownames(ss$summary)),6], NA)
)

g1 <- ggplot(fitdata) +
  geom_vline(xintercept=2016:2024, lty=2, col="gray") +
  geom_point(data=data_fit, aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(aes(year+week/52, est)) +
  geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2016:2024) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g2 <- ggplot(rdata) +
  geom_vline(xintercept=2016:2024, lty=2, col="gray") +
  geom_line(aes(year+week/52, est)) +
  geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2016:2024) +
  scale_y_continuous("Growth rate (1/week)")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

g3 <- ggplot(phasedata) +
  geom_path(aes(r, c, col=year+week/52)) +
  scale_x_continuous("Growth rate (1/week)") +
  scale_y_log10("Predicted cases") +
  scale_color_viridis_c("Year") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

gcomb <- ggarrange(g1, g2, g3)
ggsave("fig1.pdf", gcomb, width=8, height=12)

phasedata_mean <- phasedata %>%
  group_by(week) %>%
  mutate(
    r=mean(r[year<2020]),
    c=exp(mean(log(c[year<2020])))
  )

mp <- as.matrix(phasedata)[,3:4]
mp[,1] <- log(mp[,1])
mup <- as.matrix(phasedata_mean)[,3:4]
mup[,1] <- log(mup[,1])

dist <- distfun(
  mat_perturb = mp,
  mat_unperturb = mup,
  time = phasedata$year+phasedata$week/52,
  tbreak = 2020+1/52,
  out="all"
)

g4 <- ggplot(dist %>% filter(time >= 2020)) +
  geom_line(aes(time, dist)) +
  geom_smooth(aes(time, dist)) +
  scale_x_continuous("Year") +
  scale_y_log10("Distance")

ggsave("fig2.pdf", g4, width=8, height=6)

phasedata_mean2 <- phasedata %>%
  mutate(even=year %% 2) %>%
  group_by(even, week) %>%
  mutate(
    r=mean(r[year<2020]),
    c=exp(mean(log(c[year<2020])))
  )

mup2 <- as.matrix(phasedata_mean2)[,3:4]
mup2[,1] <- log(mup2[,1])

dist2 <- distfun(
  mat_perturb = mp,
  mat_unperturb = mup2,
  time = phasedata$year+phasedata$week/52,
  tbreak = 2020+1/52,
  out="all"
)

g5 <- ggplot(fitdata) +
  geom_vline(xintercept=2016:2024, lty=2, col="gray") +
  # geom_point(data=data_fit, aes(year+week/52, cases), col="#EF6351", shape=1) +
  geom_line(data=phasedata_mean2, aes(year+week/52, c), col="blue", lwd=2, alpha=0.4) +
  geom_line(aes(year+week/52, est)) +
  geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), alpha=0.4) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2016:2024) +
  scale_y_continuous("Cases", limits=c(0, NA), expand=c(0, 0))  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

maxt <- dist2$time[which.max(dist2$dist)]

lfit <- lm(log(dist)~time, data=filter(dist2, time>=maxt))

-1/coef(lfit)[2]

-1/confint(lfit)[2,]

g6 <- ggplot(dist2 %>% filter(time >= 2020)) +
  geom_vline(xintercept=maxt, lty=2) +
  geom_line(aes(time, dist)) +
  geom_smooth(aes(time, dist)) +
  geom_smooth(data=dist2 %>% filter(time >= maxt), aes(time, dist), method="lm",
              col="red") +
  scale_x_continuous("Year") +
  scale_y_log10("Distance")  +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth=1)
  )

gcomb2 <- ggarrange(g5, g6)

ggsave("fig3.pdf", gcomb2, width=8, height=6)
