library(dplyr)
source("../R/simulate_stochastic.R")

R0vec <- seq(2, 20, length.out=21)
deltavec <- exp(seq(log(1/365/40), log(1/365), length.out=21))

gamma <- 1/7
N <- 1e8

paramdata <- expand.grid(R0vec, deltavec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  print(i)
  pp <- paramdata[i,]
  
  R0 <- pp[[1]]
  delta <- pp[[2]]
  
  beta <- R0*gamma
  S <- 1/R0
  Istar <- delta * (1-S)/(beta * S + delta)
  
  resilience <- (delta+beta*Istar)/2
  
  imaginary <- sqrt(4*(delta+beta*S)*beta*Istar-(delta+beta*Istar)^2)/2
  period <- 2 * pi/imaginary
  
  ss_det <- simulate_sirs_deterministic(R0=R0, delta=delta, theta=0.04, N=N, sigma_env=0, rho_env=0,
                           I0=Istar, rho=1, k=10,
                           tmax=100)
  
  ss <- simulate_sirs_stochastic(R0=R0, delta=delta, theta=0.04, N=N, sigma_env=0, rho_env=0,
                                 I0=Istar, rho=1, k=10,
                                 tmax=100)
  
  ss_diff <- (ss$I-ss_det$I)/ss_det$I
  amplitude <- (max(ss_diff[ss$time>25]) - min(ss_diff[ss$time>25]))/2
  
  reslist[[i]] <- data.frame(
    R0=R0,
    delta=delta,
    resilience=resilience,
    period=period/365,
    amplitude=amplitude
  )
}

analysis_factorial_noise_cycle <- reslist %>%
  bind_rows

save("analysis_factorial_noise_cycle", 
     file="analysis_factorial_noise_cycle.rda")
