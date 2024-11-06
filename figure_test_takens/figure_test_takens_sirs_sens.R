library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/distfun.R")
source("../R/takens.R")

ss_perturb <- simulate_sirs(tmax=2030) %>%
  filter(time >= 2010) %>%
  mutate(
    group=time < 2020
  )

ss_unperturb <- simulate_sirs(tmax=2030,npifun=function(x)1) %>%
  filter(time >= 2010)

time <- ss_perturb$time

delayvec <- c(7, 14, 7*4, 7*13)

reslist <- vector('list', length(delayvec))

for (i in 1:length(delayvec)) {
  delay <- delayvec[i]
  maxdimension <- pmin(floor(365/delay), 8)
  
  out <- lapply(2:maxdimension, function(j){
    mat_perturb_takens <- takens(log(ss_perturb$I), d=j, tau=delay)
    mat_unperturb_takens <- takens(log(ss_unperturb$I), d=j, tau=delay)
    
    time_takens <- tail(ss_perturb$time, -delay*(j-1))
    
    dist_takens <- distfun(mat_perturb_takens, mat_unperturb_takens,
                           time_takens,
                           2020,
                           out="all")
    
    
    
    maxt <- dist_takens$time[which.max(dist_takens$dist)]
    lfit <- lm(log(dist)~time, data=filter(dist_takens, time>=maxt))
    
    data.frame(
      time=time_takens,
      dist=dist_takens$dist,
      delay=delay,
      dimension=j,
      returntime=1/-coef(lfit)[[2]]
    )    
  }) %>%
    bind_rows
  
  
  reslist[[i]] <- out
}

resdata <- reslist %>%
  bind_rows %>%
  mutate(
    delay=factor(delay,
                 levels=delayvec,
                 labels=paste0(delayvec, "-day lag")),
    dimension=paste0(dimension, " dimensions")
  )

mat_perturb_SI <- matrix(c(ss_perturb$S, log(ss_perturb$I)), ncol=2)
mat_unperturb_SI <- matrix(c(ss_unperturb$S, log(ss_unperturb$I)), ncol=2)

dist_SI <- distfun(mat_perturb_SI, mat_unperturb_SI,
                   ss_perturb$time,
                   2020, "all")

ee <- eigen_sirs()

maxt_SI <- dist_SI$time[which.max(dist_SI$dist)]
lfit_SI <- lm(log(dist)~time, data=filter(dist_SI, time>=maxt_SI))

g1 <- ggplot(resdata) +
  geom_line(data=filter(dist_SI, time>=2020), aes(time, dist), col="#224B95", alpha=0.2, lwd=1) +
  geom_function(fun=function(x) exp(predict(lfit_SI)[1]) * exp(-ee*(x-maxt_SI)),
                lwd=0.7, lty=2) +
  geom_line(aes(time, dist)) +
  scale_x_continuous("Year", breaks=seq(2020, 2030, 2),) +
  scale_y_log10("Distance from attractor", expand=c(0, 0)) +
  coord_cartesian(xlim=c(2020, 2030), ylim=c(1e-5, 10)) +
  facet_grid(delay~dimension) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank()
  )

resdata_summ <- resdata %>%
  group_by(delay, dimension) %>%
  summarize(
    returntime=unique(returntime)
  ) %>%
  bind_rows(
    data.frame(
      delay="SI",
      dimension="2 dimensions",
      returntime=-1/coef(lfit_SI)[[2]]
    )
  ) %>%
  mutate(
    delay=factor(delay,
                 levels=c("SI", "7-day lag", "14-day lag",
                          "28-day lag", "91-day lag"))
  )

g2 <- ggplot(resdata_summ) +
  geom_point(aes(delay, returntime, col=dimension)) +
  geom_hline(yintercept=1/ee, lty=2) +
  scale_x_discrete("Delays") +
  scale_y_continuous("Return time (years)") +
  scale_color_discrete("Dimensions")

gcomb <- ggarrange(g1, g2, nrow=2, heights=c(4, 1))

ggsave("figure_test_takens_sirs_sens.pdf", gcomb, width=15, height=10)
