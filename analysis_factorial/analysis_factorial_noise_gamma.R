library(dplyr)
source("../R/simulate_stochastic.R")

R0 <- 2
deltavec <- exp(seq(log(1/365/40), log(1/365), length.out=21))
gammavec <- exp(seq(log(1/7/4), log(1/7), length.out=21))

N <- 1e8

paramdata <- expand.grid(gammavec, deltavec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  print(i)
  pp <- paramdata[i,]
  
  gamma <- pp[[1]]
  delta <- pp[[2]]
  
  # approximate
  # period <- 2 * pi/(sqrt(gamma*delta*(R0-1))) # in days
  
  beta <- R0*gamma
  S <- 1/R0
  Istar <- delta * (1-S)/(beta * S + delta)
  
  resilience <- (delta+beta*Istar)/2
  
  imaginary <- sqrt(4*(delta+beta*S)*beta*Istar-(delta+beta*Istar)^2)/2
  period <- 2 * pi/imaginary
  
  ss_det <- simulate_sirs_deterministic(R0=R0, gamma=gamma, delta=delta, theta=0, N=N, sigma_env=0, rho_env=0,
                                        I0=Istar, rho=1, k=10,
                                        tmax=400)
  
  S0 <- tail(ss_det$S,1)/N
  I0 <- tail(ss_det$I,1)/N
  
  ss <- simulate_sirs_stochastic(R0=R0, gamma=gamma, delta=delta, theta=0, N=N, sigma_env=0, rho_env=0,
                                 S0=S0, I0=I0, rho=1, k=10,
                                 tmax=100)
  
  amplitude <- (max(filter(ss, time>25)$I)-min(filter(ss, time>25)$I))/mean(filter(ss, time>25)$I)/2
  
  spec <- spectrum(log(filter(ss, time > 25)$I), plot=FALSE)
  
  spx <- spec$freq*365
  spy <- rev(2*spec$spec)
  
  px <- rev(1/spx)
  
  reslist[[i]] <- data.frame(
    R0=R0,
    delta=delta,
    gamma=gamma,
    resilience=resilience,
    period=period/365,
    period_obs=px[which.max(spy)],
    amplitude=amplitude
  )
}

analysis_factorial_noise_gamma <- reslist %>%
  bind_rows

save("analysis_factorial_noise_gamma", 
     file="analysis_factorial_noise_gamma.rda")
