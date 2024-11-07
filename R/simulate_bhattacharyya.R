npifun_default <- function(t) {
  if (t >= 2020 & t < 2020.5) {
    npi <- 0.5
  } else {
    npi <- 1
  }
  
  return(npi)
}

model_bhattacharyya <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    beta1 <- b01 * (1 + theta1 * sin(2*pi*(t-phi1))) * q1 * npifun(t)
    beta2 <- b02 * (1 + theta2 * sin(2*pi*(t-phi2))) * q2 * npifun(t)
    
    prevalence1 <- I1 + J1
    prevalence2 <- I2 + J2
    
    lambda1 <- beta1 * prevalence1
    lambda2 <- beta2 * prevalence2
    
    dS <- mu - mu * S - (lambda1 + lambda2) * S + rho1 * R1 + rho2 * R2
    dI1 <- lambda1 * S - (gamma1 + mu) * I1
    dI2 <- lambda2 * S - (gamma2 + mu) * I2
    dR1 <- gamma1 * I1 - lambda2 * epsilon21 * R1 - rho1 * R1 + rho2 * R - mu * R1
    dR2 <- gamma2 * I2 - lambda1 * epsilon12 * R2 - rho2 * R2 + rho1 * R - mu * R2
    dJ1 <- lambda1 * epsilon12 * R2 - (gamma1 + mu) * J1
    dJ2 <- lambda2 * epsilon21 * R1 - (gamma2 + mu) * J2
    dR <- gamma1 * J1 + gamma2 * J2 - rho1 * R - rho2 * R - mu * R
    
    Seff1 <- S + R2 * epsilon12
    Seff2 <- S + R1 * epsilon21
    Sall <- S + R2 + R1
    
    list(c(dS, dI1, dI2, dR1, dR2, dJ1, dJ2, dR),
         prevalence1=prevalence1,
         prevalence2=prevalence2,
         Seff1=Seff1,
         Seff2=Seff2,
         Sall=Sall,
         Reff1=beta1/(gamma1+mu)*Seff1,
         Reff2=beta2/(gamma2+mu)*Seff2)
  })
}

simulate_bhattacharyya <- function(b01=3.4*52,
                                   b02=3.9*52,
                                   theta1=0.4,
                                   theta2=0.3,
                                   phi1=0.005*7/365,
                                   phi2=4.99*7/365,
                                   epsilon12=0.92,
                                   epsilon21=0.45,
                                   gamma1=365/10,
                                   gamma2=365/10,
                                   rho1=1,
                                   rho2=1,
                                   q1=0.5,
                                   q2=0.5,
                                   mu=1/70,
                                   yini,
                                   tmin=1900,
                                   tmax=2030,
                                   npifun=npifun_default) {
  parms <- c(b01=b01, b02=b02, theta1=theta1, theta2=theta2, phi1=phi1, phi2=phi2, epsilon12=epsilon12, epsilon21=epsilon21,
             gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
             npifun=npifun)
  
  if (missing(yini)) {
    yini <- c(S=1-2e-6, I1=1e-6, I2=1e-6, R1=0, R2=0, J1=0, J2=0, R=0)
  }
  
  times <- seq(tmin, tmax, by=1/365)
  
  out <- as.data.frame(ode(yini, times, model_bhattacharyya, parms))
  
  out$r1 <- c(NA, diff(out$prevalence1))/out$prevalence1
  out$r2 <- c(NA, diff(out$prevalence2))/out$prevalence2
  
  out
}

eigen_bhattacharyya <- function(b01=3.4*52,
                                b02=3.9*52,
                                theta1=0.4,
                                theta2=0.3,
                                phi1=0.005*7/365,
                                phi2=4.99*7/365,
                                epsilon12=0.92,
                                epsilon21=0.45,
                                gamma1=365/10,
                                gamma2=365/10,
                                rho1=1,
                                rho2=1,
                                q1=0.5,
                                q2=0.5,
                                mu=1/70) {
  parms <- c(b01=b01, b02=b02, theta1=0, theta2=0, phi1=0, phi2=0, epsilon12=epsilon12, epsilon21=epsilon21,
             gamma1=gamma1, gamma2=gamma2, rho1=rho1, rho2=rho2, q1=q1, q2=q2, mu=mu,
             npifun=function(x) 1)
  
  yini <- c(S=1-2e-6, I1=1e-6, I2=1e-6, R1=0, R2=0, J1=0, J2=0, R=0)
  
  times <- seq(1900, 2000, by=1/365)
  
  out <- as.data.frame(deSolve::ode(yini, times, model_bhattacharyya, parms))
  
  equi <- unlist(tail(out,1)[2:9])
  
  function_deriv_only <- function(y) {
    model_bhattacharyya(0, y, parms)[[1]]
  }
  
  jac <- numDeriv::jacobian(function_deriv_only, equi)
  
  ee <- eigen(jac)
  
  ee  
}
