library(dplyr)
library(deSolve)
library(bbmle)
source("../R/simulate_sir.R")
load("../data_simulation/data_simulation_sirs.rda")
 
nllfun <- function(log.b_1=8.0842508,
                   logit.theta=-2.1653668,
                   logit.phi=0.3229489,
                   logit.reduction=1.1,
                   logit.S0=-3.7458692,
                   logit.I0=-7.6473678,
                   logit.rho=0.6929699,
                   log.size=-0.1767736) {
  npifun_fit <- function(t) {
    if (t >= 2020 & t < 2020.5) {
      npi <- 1-plogis(logit.reduction)
    } else {
      npi <- 1
    }
    
    return(npi)
  }
  
  S0 <- plogis(logit.S0)
  I0 <- plogis(logit.I0)
  rho <- plogis(logit.rho)
  size <- exp(log.size)
  
  ss <- try(simulate_sir(
    b_1=exp(log.b_1),
    theta=plogis(logit.theta),
    phi=plogis(logit.phi)-0.5,
    mu=1/50,
    gamma=365/7,
    npifun=npifun_fit,
    yini=c(S=S0, I=I0, R=1-S0-I0, C=0) ,
    tmin=head(data_simulation_sirs$time,1),
    tmax=tail(data_simulation_sirs$time,1)
  ))
  
  if (inherits(ss, "try-error")) {
    return(NA)
  }
  
  cases_fit <- c(NA, diff(ss$C) * 1e6)
  
  plot(data_simulation_sirs$cases, log="y")
  lines(cases_fit*rho, col=2)
  
  nll <- -sum(dnbinom(data_simulation_sirs$cases, mu=cases_fit*rho, size=size,log=TRUE),
       na.rm=TRUE)
  
  print(nll)
  
  nll
}

mfit <- mle2(nllfun,
             method="Nelder-Mead")


