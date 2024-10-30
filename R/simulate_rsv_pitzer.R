npifun_default <- function(t) {
  if (t >= 2020 & t < 2020.5) {
    npi <- 0.5
  } else {
    npi <- 1
  }
  
  return(npi)
}

model_rsv_pitzer <- function(t, y, par) {
  with(as.list(c(y, par)), {
    beta <- b_1 * (1 + theta * cos(2 * pi * (t-phi))) * npifun(t)
    
    lambda <- beta * (I1 + rho1 * I2 + rho2 * (I3 + I4))
    
    dM <- mu - (omega+mu) * M
    dS0 <- omega * M - (lambda + mu) * S0
    dI1 <- lambda * S0 - (gamma1 + mu) * I1
    dS1 <- gamma1 * I1 - (sigma1 * lambda + mu) * S1
    dI2 <- sigma1 * lambda * S1 - (gamma2 + mu) * I2
    dS2 <- gamma2 * I2 - (sigma2 * lambda + mu) * S2
    dI3 <- sigma2 * lambda * S2 - (gamma3 + mu) * I3
    dS3 <- gamma3 * I3 - (sigma3 * lambda + mu) * S3 + gamma3 * I4
    dI4 <- sigma3 * lambda * S3 - (gamma3 + mu) * I4
    
    R0 <- beta/(gamma1+mu)
    
    trans_eff <- (I1 + rho1 * I2 + rho2 * (I3 + I4))/(I1+I2+I3+I4)
    
    sigma_D2 <- 1/(gamma2+mu)/(1/(gamma1+mu))
    sigma_D3 <- 1/(gamma3+mu)/(1/(gamma1+mu))
    
    dur_eff <- (I1+I2+I3+I4)/(I1+I2/sigma_D2+I3/sigma_D3+I4/sigma_D3)
    
    Seff <- S0 + sigma1 * S1 + sigma2 * S2 + sigma3 * S3
    
    Reff <- R0 * trans_eff * dur_eff * Seff
    
    list(c(dM, dS0, dI1, dS1, dI2, dS2, dI3, dS3, dI4),
         beta=beta,
         Reff=Reff,
         Seff=Seff,
         prevalence=I1+I2+I3+I4)
  })
}

simulate_rsv_pitzer <- function(R0=9,
                                theta=0.2,
                                phi=-0.1,
                                omega=1/16/7*365,
                                gamma1=1/10*365,
                                gamma2=1/7*365,
                                gamma3=1/5*365,
                                sigma1=0.76,
                                sigma2=0.6,
                                sigma3=0.4,
                                rho1=0.75,
                                rho2=0.51,
                                mu=1/80,
                                npifun=npifun_default,
                                yini,
                                tmin=1900,
                                tmax=2030) {
  b_1 <- R0 * (gamma1+mu)
  
  parms <- c(b_1=b_1,
             theta=theta,
             phi=phi,
             omega=omega,
             gamma1=gamma1,
             gamma2=gamma2,
             gamma3=gamma3,
             sigma1=sigma1,
             sigma2=sigma2,
             sigma3=sigma3,
             rho1=rho1,
             rho2=rho2,
             mu=mu,
             npifun=npifun)
  
  if (missing(yini)) {
    yini <- c(M=0, S0=1/R0-1e-6, I1=1e-6, S1=1-1/R0, I2=0, S2=0, I3=0, S3=0, I4=0)
  }
  
  times <- seq(tmin, tmax, by=1/365)
  
  out <- as.data.frame(deSolve::ode(yini, times, model_rsv_pitzer, parms))
  
  out$r <- c(NA, diff(out$prevalence))/out$prevalence
  
  out
}

eigen_rsv_pitzer <- function(R0=9,
                             omega=1/16/7*365,
                             gamma1=1/10*365,
                             gamma2=1/7*365,
                             gamma3=1/5*365,
                             sigma1=0.76,
                             sigma2=0.6,
                             sigma3=0.4,
                             rho1=0.75,
                             rho2=0.51,
                             mu=1/80) {
  b_1 <- R0 * (gamma1+mu)
  
  parms <- c(b_1=b_1,
             theta=0,
             phi=0,
             omega=omega,
             gamma1=gamma1,
             gamma2=gamma2,
             gamma3=gamma3,
             sigma1=sigma1,
             sigma2=sigma2,
             sigma3=sigma3,
             rho1=rho1,
             rho2=rho2,
             mu=mu,
             npifun=function(x) 1)
  
  yini <- c(M=0, S0=1/R0-1e-6, I1=1e-6, S1=1-1/R0, I2=0, S2=0, I3=0, S3=0, I4=0)
  
  times <- seq(1900, 2100, by=1/365)
  
  out <- as.data.frame(deSolve::ode(yini, times, model_rsv_pitzer, parms))
  
  equi <- unlist(tail(out,1)[2:10])
  
  function_deriv_only <- function(y) {
    model_rsv_pitzer(0, y, parms)[[1]]
  }
  
  jac <- numDeriv::jacobian(function_deriv_only, equi)
  
  ee <- eigen(jac)
  
  ee
}
