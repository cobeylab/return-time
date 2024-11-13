library(dplyr)
library(deSolve)
library(bbmle)
source("../R/simulate_sirs.R")
load("../data_simulation/data_simulation_sirs.rda")
 
nllfun <- function(log.b_1=4.59612970,
                   logit.theta=-1.33201283,
                   logit.phi=-0.11612715,
                   logit.reduction=-0.25463844,
                   log.delta=0.02786531,
                   logit.S0=0.99288075,
                   logit.I0=-5.95103070,
                   logit.rho=-3.94463813,
                   log.size=1.78602304) {
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
  delta <- exp(log.delta)
  
  ss <- try(simulate_sirs(
    b_1=exp(log.b_1),
    theta=plogis(logit.theta),
    phi=plogis(logit.phi)-0.5,
    mu=1/50,
    gamma=365/7,
    delta=delta,
    npifun=npifun_fit,
    yini=c(S=S0, I=I0, R=1-S0-I0, C=0) ,
    tmin=head(data_simulation_sirs$time,1),
    tmax=tail(data_simulation_sirs$time,1)
  ))
  
  if (inherits(ss, "try-error")) {
    return(NA)
  }
  
  cases_fit <- c(NA, diff(ss$C) * 1e6)
  
  plot(data_simulation_sirs$cases)
  lines(cases_fit*rho, col=2)
  
  nll <- -sum(dnbinom(data_simulation_sirs$cases, mu=cases_fit*rho, size=size,log=TRUE),
       na.rm=TRUE)
  
  if(nll==0) {
    return(NA)
  }
  
  print(nll)
  
  nll
}

fit_sirs_sirs <- mle2(nllfun, method="BFGS",
                      skip.hessian = TRUE)

save(fit_sirs_sirs, file="fit_sirs_sirs.rda")
