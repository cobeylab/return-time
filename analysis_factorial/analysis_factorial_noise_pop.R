library(dplyr)
source("../R/simulate_stochastic.R")

R0 <- 2
deltavec <- exp(seq(log(1/365/40), log(1/365), length.out=21))
Nvec <- exp(seq(log(1e6), log(1e8), length.out=21))

gamma <- 1/7

paramdata <- expand.grid(Nvec, deltavec)

reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  print(i)
  pp <- paramdata[i,]
  
  N <- pp[[1]]
  delta <- pp[[2]]
  
  beta <- R0*gamma
  S <- 1/R0
  Istar <- delta * (1-S)/(beta * S + delta)
  
  resilience <- (delta+beta*Istar)/2
  
  imaginary <- sqrt(4*(delta+beta*S)*beta*Istar-(delta+beta*Istar)^2)/2
  period <- 2 * pi/imaginary
  
  ss_det <- simulate_sirs_deterministic(R0=R0, delta=delta, theta=0, N=N, sigma_env=0, rho_env=0,
                                        I0=Istar, rho=1, k=10,
                                        tmax=400)
  
  S0 <- tail(ss_det$S,1)/N
  I0 <- tail(ss_det$I,1)/N
  
  ss <- simulate_sirs_stochastic(R0=R0, delta=delta, theta=0, N=round(N), sigma_env=0, rho_env=0,
                                 S0=S0, I0=I0, rho=1, k=10,
                                 tmax=100)
  
  extinction <- any(ss$I==0)
  
  amplitude <- (max(filter(ss, time>25)$I)-min(filter(ss, time>25)$I))/mean(filter(ss, time>25)$I)/2
  
  reslist[[i]] <- data.frame(
    R0=R0,
    delta=delta,
    N=N,
    resilience=resilience,
    period=period/365,
    amplitude=amplitude,
    extinction=extinction
  )
}

analysis_factorial_noise_pop <- reslist %>%
  bind_rows

save("analysis_factorial_noise_pop", 
     file="analysis_factorial_noise_pop.rda")
