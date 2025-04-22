library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(ggrepel)
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_seir.R")
load("../figure4/figure4_summ.rda")

R0vec <- seq(1.1, 6,
             length.out=101)

deltavec <- exp(seq(log(1/4), log(2), length.out=101))

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
  R0 <- 6
  delta <- delta
  mu <- 1/80
  gamma <- 365/7
  b_1 <- R0*(gamma+mu)
  
  S <- 1/R0
  
  I <- (delta * (1-S) + mu - mu * S)/(delta + b_1 * S)
  
  R <- (1-S-I)
  
  resilience - (delta + mu + b_1 * I)/2
}

analysis_all_summ_mean <- analysis_all_summ_filter %>%
  group_by(key) %>%
  summarize(
    resilience=mean(resilience)
  ) %>%
  arrange(resilience) %>%
  group_by(key) %>%
  mutate(
    delta=uniroot(resopt, interval=c(0.01, 20),
                  resilience=resilience)$root,
    R0=10
  )

# measles_resilience <- -eigen_seir()

g1 <- ggplot(summdata) +
  geom_raster(aes(R0, 1/delta, fill=log10(resilience))) +
  geom_contour(aes(R0, 1/delta, z=log10(resilience)),
               breaks=log10(analysis_all_summ_mean$resilience),
               col="white") +
  geom_text_repel(data=analysis_all_summ_mean, aes(6.04, 1/delta, label=key,
                                                   col=log10(resilience)),
                  direction = "y",
                  segment.color = NA,
                  hjust=0) +
  scale_x_continuous(expression("Basic reproduction number,"~R[0]~"                                                     "), expand=c(0, 0),
                     breaks=c(2, 4, 6),
                     limits = c(NA,9)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0),
                breaks=c(0.5, 1, 2, 4)) +
  scale_fill_viridis_c("Resilience\n(1/years)", limits=log10(range(summdata$resilience)),
                       breaks=log10(c(0.25, 0.5, 1, 2, 4)),
                       labels=c(0.25, 0.5, 1, 2, 4)) +
  scale_color_viridis_c("Resilience\n(1/years)", limits=log10(range(summdata$resilience)),
                        breaks=log10(c(0.25, 0.5, 1, 2, 4)),
                        labels=c(0.25, 0.5, 1, 2, 4)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

g2 <- ggplot(summdata) +
  geom_raster(aes(R0, 1/delta, fill=log10(replenish*100))) +
  geom_contour(aes(R0, 1/delta, z=log10(resilience)),
               breaks=log10(analysis_all_summ_mean$resilience),
               col="white") +
  geom_text_repel(data=analysis_all_summ_mean, aes(6.04, 1/delta, label=key),
                  direction = "y",
                  segment.color = NA,
                  hjust=0) +
  scale_x_continuous(expression("Basic reproduction number,"~R[0]~"                                                     "), expand=c(0, 0),
                     breaks=c(2, 4, 6),
                     limits = c(NA,9)) +
  scale_y_log10("Duration of immunity (years)", expand=c(0, 0),
                breaks=c(0.5, 1, 2, 4)) +
  scale_fill_viridis_c("Susceptible replenishment rate\n(%/years)",
                       option="A",
                       breaks=log10(c(6.25, 25, 100, 400)),
                       labels=c(6.25, 25, 100, 400)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

replenish_range <- lapply(split(analysis_all_summ_mean, analysis_all_summ_mean$key)[-7], function(x) {
  summdata %>%
    filter(
      resilience < x$resilience * 1.01 & resilience > x$resilience * 0.99
    ) %>%
    mutate(
      key=x$key[1]
    )
}) %>%
  bind_rows

range(filter(replenish_range, !key %in% c("Bocavirus", "Norovirus"))$replenish)

1/50 * 17

gcomb <- ggarrange(g1, g2, nrow=1,
                   labels=c("A", "B"))

ggsave("figure5.pdf", gcomb, width=12, height=6)
