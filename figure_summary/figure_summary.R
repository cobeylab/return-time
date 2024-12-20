library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(ggrepel)
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
load("../figure_takens/figure_takens_acf_fnn_pred.rda")

R0vec <- seq(1.1, 3,
             length.out=101)

deltavec <- exp(seq(log(1/2.5), log(2), length.out=101))

paramdata <- expand.grid(R0vec, deltavec)

summdata <- apply(paramdata, 1, function(x) {
  R0 <- x[[1]]
  delta <- x[[2]]
  mu <- 1/80
  gamma <- 365/7
  b_1 <- R0*(gamma+mu)
  
  S <- 1/R0
  
  I <- (delta * (1-S) + mu - mu * S)/(delta + b_1 * S)
  
  R <- (1-S-I)
  
  resilience <- ((delta + mu + b_1 * I)/2)
  
  data.frame(
    R0=R0,
    delta=delta,
    resilience=resilience,
    replenish=(delta*R+mu)/S
  )
}, simplify=FALSE) %>%
  bind_rows

resopt <- function(delta, resilience) {
  R0 <- 3
  delta <- delta
  mu <- 1/80
  gamma <- 365/7
  b_1 <- R0*(gamma+mu)
  
  S <- 1/R0
  
  I <- (delta * (1-S) + mu - mu * S)/(delta + b_1 * S)
  
  R <- (1-S-I)
  
  resilience - (delta + mu + b_1 * I)/2
}

analysis_all_summ_mean <- analysis_all_summ %>%
  group_by(key) %>%
  summarize(
    resilience=mean(resilience)
  ) %>%
  arrange(resilience) %>%
  group_by(key) %>%
  mutate(
    delta=uniroot(resopt, interval=c(0.01, 3),
                  resilience=resilience)$root,
    R0=3
  )

g1 <- ggplot(summdata) +
  geom_raster(aes(R0, 1/delta, fill=resilience)) +
  geom_contour(aes(R0, 1/delta, z=resilience),
               breaks=analysis_all_summ_mean$resilience,
               col="white") +
  geom_text_repel(data=analysis_all_summ_mean, aes(3.04, 1/delta, label=key,
                                                   col=resilience),
                  direction = "y",
                  segment.color = NA,
                  hjust=0) +
  scale_x_continuous("Basic reproduction number", expand=c(0, 0),
                     breaks=c(1.5, 2, 2.5, 3),
                     limits = c(NA,4)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Resilience\n(1/years)", limits=range(summdata$resilience)) +
  scale_color_viridis_c("Resilience\n(1/years)", limits=range(summdata$resilience)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

g2 <- ggplot(summdata) +
  geom_raster(aes(R0, 1/delta, fill=replenish*100)) +
  geom_contour(aes(R0, 1/delta, z=resilience),
               breaks=analysis_all_summ_mean$resilience,
               col="white") +
  geom_text_repel(data=analysis_all_summ_mean, aes(3.04, 1/delta, label=key),
                  direction = "y",
                  segment.color = NA,
                  hjust=0) +
  scale_x_continuous("Basic reproduction number", expand=c(0, 0),
                     breaks=c(1.5, 2, 2.5, 3),
                     limits = c(NA,4)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0)) +
  scale_fill_viridis_c("Susceptible replenishment rate\n(%/years)",
                       option="A") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

gcomb <- ggarrange(g1, g2, nrow=1,
                   labels=c("A", "B"))

ggsave("figure_summary.pdf", gcomb, width=12, height=5)
