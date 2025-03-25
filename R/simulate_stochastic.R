simulate_sir_stochastic <- function(R0=17,
                                    gamma=1/7,
                                    mu=1/365/50,
                                    theta=0.2,
                                    phi=0,
                                    omega=0.01,
                                    N=1e7,
                                    sigma_env=0.02,
                                    rho_env=0,
                                    S0,
                                    I0=1e-5,
                                    rho=0.1,
                                    k=100,
                                    tmax=50,
                                    seed) {
  if (!missing(seed)) {
    set.seed(seed)
  }
  time <- (1:(tmax*365))/365
  Svec <- Ivec <- Rvec <- Nvec <- Cvec <- Cexpectedvec <- betavec <- noisevec <- betavec_realized <- rep(0, tmax*365)
  
  if (!missing(S0)) {
    Svec[1] <- round(1/(R0) * N)
  } else {
    Svec[1] <- round(S0 * N)
  }
  
  Ivec[1] <- round(I0 * N)
  Rvec[1] <- round(N-Svec[1]-Ivec[1])
  
  Nvec[1] <- Svec[1] + Ivec[1] + Rvec[1]
  
  betavec[1] <- R0 * (1 + theta * cos(2 * pi * (1/365+phi))) * (1-exp(-(gamma+mu)))
  noisevec[1] <- rnorm(1, 0, sigma_env)
  betavec_realized[1] <- betavec[1] * exp(noisevec[1])
  
  for (i in 2:(tmax*365)) {
    beta <- betavec[i] <- R0 * (1 + theta * cos(2 * pi * (i/365+phi))) * (1-exp(-(gamma+mu)))
    
    noisevec[i] <- rnorm(1, mean = rho_env * noisevec[i-1], sd=sqrt(1-rho_env^2) * sigma_env)
    
    betavec_realized[i] <- beta * exp(noisevec[i])
    
    lambda_expected <- beta * (Ivec[i-1]+omega)/Nvec[i-1]
    
    lambda <- betavec_realized[i] * (Ivec[i-1]+omega)/Nvec[i-1]
    
    birth <- rpois(1, mu * Nvec[i-1])
    
    Sout <- rbinom(1, Svec[i-1], 1-exp(-(lambda+mu)))
    StoI <- rbinom(1, Sout, lambda/(lambda+mu))
    
    Cexpectedvec[i] <- Svec[i-1] * (1-exp(-(lambda_expected+mu))) * lambda_expected/(lambda_expected+mu)
    
    Iout <- rbinom(1, Ivec[i-1], 1-exp(-(gamma+mu)))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    
    Rout <- rbinom(1, Rvec[i-1], 1-exp(-(mu)))
    
    Svec[i] <- Svec[i-1] + birth - Sout
    Ivec[i] <- Ivec[i-1] - Iout + StoI
    Rvec[i] <- Rvec[i-1] - Rout + ItoR  
    Nvec[i] <- Svec[i] + Ivec[i] + Rvec[i] 
    Cvec[i] <- StoI
  }
  
  rvec <- c(diff(Cvec)/head(Cvec, -1), NA)
  rexpectedvec <- c((tail(Cexpectedvec,-1)-head(Cvec,-1))/head(Cvec, -1), NA)
  
  cases <- rnbinom(length(Cvec), mu=Cvec*rho, size=100)
  
  data.frame(
    time=time,
    S=Svec,
    I=Ivec,
    R=Rvec,
    N=Nvec,
    C=Cvec,
    Cexpected=Cexpectedvec,
    cases=cases,
    beta=betavec,
    noise=noisevec,
    beta_realized=betavec_realized,
    r=rvec,
    rexpected=rexpectedvec
  )
}

simulate_seir_stochastic <- function(R0=17,
                                     sigma=1/8,
                                     gamma=1/5,
                                     mu=1/365/50,
                                     theta=0.2,
                                     omega=0.01,
                                     N=1e7,
                                     sigma_env=0.02,
                                     rho_env=0,
                                     I0=1e-5,
                                     rho=0.1,
                                     k=100,
                                     tmax=50,
                                     seed) {
  if (!missing(seed)) {
    set.seed(seed)
  }
  time <- (1:(tmax*365))/365
  Svec <- Evec <- Ivec <- Rvec <- Nvec <- Cvec <- Cexpectedvec <- betavec <- noisevec <- betavec_realized <- rep(0, tmax*365)
  
  Svec[1] <- round(1/(R0) * N)
  Ivec[1] <- round(I0 * N)
  Rvec[1] <- N-Svec[1]-Ivec[1]
  
  Nvec[1] <- Svec[1] + Evec[1] + Ivec[1] + Rvec[1]
  
  betavec[1] <- R0 * (1 + theta * cos(2 * pi * 1/365)) * (1-exp(-(gamma+mu))) * (sigma+mu)/sigma
  noisevec[1] <- rnorm(1, 0, sigma_env)
  betavec_realized[1] <- betavec[1] * exp(noisevec[1])
  
  for (i in 2:(tmax*365)) {
    beta <- betavec[i] <- R0 * (1 + theta * cos(2 * pi * i/365)) * (1-exp(-(gamma+mu))) * (sigma+mu)/sigma
    
    noisevec[i] <- rnorm(1, mean = rho_env * noisevec[i-1], sd=sqrt(1-rho_env^2) * sigma_env)
    
    betavec_realized[i] <- beta * exp(noisevec[i])
    
    lambda_expected <- beta * (Ivec[i-1]+omega)/Nvec[i-1]
    
    lambda <- betavec_realized[i] * (Ivec[i-1]+omega)/Nvec[i-1]
    
    birth <- rpois(1, mu * Nvec[i-1])
    
    Sout <- rbinom(1, Svec[i-1], 1-exp(-(lambda+mu)))
    StoE <- rbinom(1, Sout, lambda/(lambda+mu))
    
    Cexpectedvec[i] <- Svec[i-1] * (1-exp(-(lambda_expected+mu))) * lambda_expected/(lambda_expected+mu)
    
    Eout <- rbinom(1, Evec[i-1], 1-exp(-(sigma+mu)))
    EtoI <- rbinom(1, Eout, sigma/(sigma+mu))
    
    Iout <- rbinom(1, Ivec[i-1], 1-exp(-(gamma+mu)))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    
    Rout <- rbinom(1, Rvec[i-1], 1-exp(-(mu)))
    
    Svec[i] <- Svec[i-1] + birth - Sout
    Evec[i] <- Evec[i-1] - Eout + StoE
    Ivec[i] <- Ivec[i-1] - Iout + EtoI
    Rvec[i] <- Rvec[i-1] - Rout + ItoR  
    Nvec[i] <- Svec[i] + Evec[i] + Ivec[i] + Rvec[i] 
    Cvec[i] <- StoE
  }
  
  rvec <- c(diff(Cvec)/head(Cvec, -1), NA)
  rexpectedvec <- c((tail(Cexpectedvec,-1)-head(Cvec,-1))/head(Cvec, -1), NA)
  
  cases <- rnbinom(length(Cvec), mu=Cvec*rho, size=100)
  
  data.frame(
    time=time,
    S=Svec,
    E=Evec,
    I=Ivec,
    R=Rvec,
    N=Nvec,
    C=Cvec,
    Cexpected=Cexpectedvec,
    cases=cases,
    beta=betavec,
    noise=noisevec,
    beta_realized=betavec_realized,
    r=rvec,
    rexpected=rexpectedvec
  )
}

simulate_sirs_demo_stochastic <- function(R0=2,
                                     gamma=1/7,
                                     delta=1/365,
                                     mu=1/365/80,
                                     theta=0.1,
                                     N=1e8,
                                     sigma_env=0.02,
                                     rho_env=0.99,
                                     I0=1e-5,
                                     rho=0.1,
                                     k=100,
                                     tmax=50) {
  time <- (1:(tmax*365))/365
  Svec <- Ivec <- Rvec <- Nvec <- Cvec <- betavec <- noisevec <- betavec_realized <- rep(0, tmax*365)
  
  Svec[1] <- round(1/R0 * N)
  Ivec[1] <- round(I0 * N)
  Rvec[1] <- N-Svec[1]-Ivec[1]
  
  Nvec[1] <- Svec[1] + Ivec[1] + Rvec[1]
  
  betavec[1] <- R0 * (1 + theta * cos(2 * pi * 1/365)) * (1-exp(-(gamma+mu)))
  noisevec[1] <- rnorm(1, 0, sigma_env)
  betavec_realized[1] <- betavec[1] * exp(noisevec[1])
  
  for (i in 2:(tmax*365)) {
    beta <- betavec[i] <- R0 * (1 + theta * cos(2 * pi * i/365)) * (1-exp(-(gamma+mu)))
    
    noisevec[i] <- rnorm(1, mean = rho_env * noisevec[i-1], sd=sqrt(1-rho_env^2) * sigma_env)
    
    betavec_realized[i] <- beta * exp(noisevec[i])
    
    lambda <- betavec_realized[i] * Ivec[i-1]/Nvec[i-1]
    
    birth <- rpois(1, mu * Nvec[i-1])
    
    Sout <- rbinom(1, Svec[i-1], 1-exp(-(lambda+mu)))
    StoI <- rbinom(1, Sout, lambda/(lambda+mu))
    
    Iout <- rbinom(1, Ivec[i-1], 1-exp(-(gamma+mu)))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    
    Rout <- rbinom(1, Rvec[i-1], 1-exp(-(delta+mu)))
    RtoS <- rbinom(1, Rout, delta/(delta+mu))
    
    Svec[i] <- Svec[i-1] + birth - Sout + RtoS
    Ivec[i] <- Ivec[i-1] - Iout + StoI
    Rvec[i] <- Rvec[i-1] - Rout + ItoR  
    Nvec[i] <- Svec[i] + Ivec[i] + Rvec[i] 
    Cvec[i] <- StoI
  }
  
  rvec <- c(diff(Cvec)/head(Cvec, -1), NA)
  
  cases <- rnbinom(length(Cvec), mu=Cvec*rho, size=100)
  
  data.frame(
    time=time,
    S=Svec,
    I=Ivec,
    R=Rvec,
    N=Nvec,
    C=Cvec,
    cases=cases,
    beta=betavec,
    noise=noisevec,
    beta_realized=betavec_realized,
    r=rvec
  )
}

simulate_sirs_stochastic <- function(R0=2,
                                     gamma=1/7,
                                     delta=1/365,
                                     theta=0.1,
                                     N=1e8,
                                     sigma_env=0.02,
                                     rho_env=0.99,
                                     S0=1/R0,
                                     I0=1e-5,
                                     rho=0.1,
                                     k=100,
                                     tmax=50) {
  time <- (1:(tmax*365))/365
  Svec <- Ivec <- Rvec <- Nvec <- Cvec <- betavec <- noisevec <- betavec_realized <- rep(0, tmax*365)
  
  Svec[1] <- round(S0 * N)
  Ivec[1] <- round(I0 * N)
  Rvec[1] <- N-Svec[1]-Ivec[1]
  
  Nvec[1] <- Svec[1] + Ivec[1] + Rvec[1]
  
  betavec[1] <- R0 * (1 + theta * cos(2 * pi * 1/365)) * (1-exp(-(gamma)))
  noisevec[1] <- rnorm(1, 0, sigma_env)
  betavec_realized[1] <- betavec[1] * exp(noisevec[1])
  
  for (i in 2:(tmax*365)) {
    beta <- betavec[i] <- R0 * (1 + theta * cos(2 * pi * i/365)) * (1-exp(-(gamma)))
    
    noisevec[i] <- rnorm(1, mean = rho_env * noisevec[i-1], sd=sqrt(1-rho_env^2) * sigma_env)
    
    betavec_realized[i] <- beta * exp(noisevec[i])
    
    lambda <- betavec_realized[i] * Ivec[i-1]/Nvec[i-1]
    
    Sout <- rbinom(1, Svec[i-1], 1-exp(-(lambda)))
    
    Iout <- rbinom(1, Ivec[i-1], 1-exp(-(gamma)))
    
    Rout <- rbinom(1, Rvec[i-1], 1-exp(-(delta)))
    
    Svec[i] <- Svec[i-1] - Sout + Rout
    Ivec[i] <- Ivec[i-1] - Iout + Sout
    Rvec[i] <- Rvec[i-1] - Rout + Iout
    Nvec[i] <- Svec[i] + Ivec[i] + Rvec[i] 
    Cvec[i] <- Sout
  }
  
  rvec <- c(diff(Cvec)/head(Cvec, -1), NA)
  
  cases <- rnbinom(length(Cvec), mu=Cvec*rho, size=100)
  
  data.frame(
    time=time,
    S=Svec,
    I=Ivec,
    R=Rvec,
    N=Nvec,
    C=Cvec,
    cases=cases,
    beta=betavec,
    noise=noisevec,
    beta_realized=betavec_realized,
    r=rvec
  )
}

simulate_sirs_deterministic <- function(R0=2,
                                     gamma=1/7,
                                     delta=1/365,
                                     theta=0.1,
                                     N=1e8,
                                     sigma_env=0.02,
                                     rho_env=0.99,
                                     S0=1/R0,
                                     I0=1e-5,
                                     rho=0.1,
                                     k=100,
                                     tmax=50) {
  time <- (1:(tmax*365))/365
  Svec <- Ivec <- Rvec <- Nvec <- Cvec <- betavec <- noisevec <- betavec_realized <- rep(0, tmax*365)
  
  Svec[1] <- S0 * N
  Ivec[1] <- I0 * N
  Rvec[1] <- N-Svec[1]-Ivec[1]
  
  Nvec[1] <- Svec[1] + Ivec[1] + Rvec[1]
  
  betavec[1] <- R0 * (1 + theta * cos(2 * pi * 1/365)) * (1-exp(-(gamma)))
  noisevec[1] <- rnorm(1, 0, sigma_env)
  betavec_realized[1] <- betavec[1] * exp(noisevec[1])
  
  for (i in 2:(tmax*365)) {
    beta <- betavec[i] <- R0 * (1 + theta * cos(2 * pi * i/365)) * (1-exp(-(gamma)))
    
    noisevec[i] <- 0
    
    betavec_realized[i] <- beta * exp(noisevec[i])
    
    lambda <- betavec_realized[i] * Ivec[i-1]/Nvec[i-1]
    
    Sout <- Svec[i-1] * (1-exp(-(lambda)))
    
    Iout <- Ivec[i-1] * (1-exp(-(gamma)))
    
    Rout <- Rvec[i-1] * (1-exp(-(delta)))
    
    Svec[i] <- Svec[i-1] - Sout + Rout
    Ivec[i] <- Ivec[i-1] - Iout + Sout
    Rvec[i] <- Rvec[i-1] - Rout + Iout
    Nvec[i] <- Svec[i] + Ivec[i] + Rvec[i] 
    Cvec[i] <- Sout
  }
  
  rvec <- c(diff(Cvec)/head(Cvec, -1), NA)
  
  cases <- rnbinom(length(Cvec), mu=Cvec*rho, size=100)
  
  data.frame(
    time=time,
    S=Svec,
    I=Ivec,
    R=Rvec,
    N=Nvec,
    C=Cvec,
    cases=cases,
    beta=betavec,
    noise=noisevec,
    beta_realized=betavec_realized,
    r=rvec
  )
}
