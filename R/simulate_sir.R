npifun_default <- function(t) {
  if (t >= 2020 & t < 2020.5) {
    npi <- 0.5
  } else {
    npi <- 1
  }
  
  return(npi)
}

model_sir <- function(t, y, par) {
  with(as.list(c(y, par)), {
    beta <- b_1 * (1 + theta * cos(2 * pi * (t-phi))) * npifun(t)
    
    lambda <- beta * I
    
    dS <- mu - (lambda + mu) * S
    dI <- lambda * S - (gamma + mu) * I
    dR <- gamma * I - mu * R
    
    list(c(dS, dI, dR),
         beta=beta,
         Reff=beta/(gamma+mu)*S,
         r=beta*S-(gamma+mu))
  })
}

simulate_sir <- function(b_1=17*(365/14+1/50),
                         theta=0.2,
                         phi=0,
                         mu=1/50,
                         gamma=365/14,
                         npifun=npifun_default,
                         yini,
                         tmin=1900,
                         tmax=2100) {
  parms <- c(b_1=b_1,
             theta=theta,
             phi=phi,
             mu=mu,
             gamma=gamma,
             npifun=npifun)
  
  R0 <- b_1/(gamma+mu)
  
  if (missing(yini)) {
    yini <- c(S=1/R0, I=1e-6, R=1-1/R0-1e-6)
  }
  
  times <- seq(tmin, tmax, by=1/365)
  
  out <- as.data.frame(deSolve::ode(yini, times, model_sir, parms))
  
  out
}
