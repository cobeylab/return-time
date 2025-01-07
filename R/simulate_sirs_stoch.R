npifun_default <- function(t) {
  if (t >= 2020 & t < 2020.5) {
    npi <- 0.5
  } else {
    npi <- 1
  }
  
  return(npi)
}

simulate_SIRS_stoch <- function(R0=3,
                                omega=0,
                                theta=0.1,
                                npifun=npifun_default,
                                gamma=1/7,
                                delta=1/52/7/2,
                                mu=1/364/50,
                                pop=1e8,
                                S0=1/3,
                                I0=1e-6,
                                rho=0.002,
                                theta_obs=1000,
                                tstart=1900,
                                tend=2030,
                                seed=101) {
  set.seed(seed)
  time <- seq(tstart, tend, by=1/52/7)
  tmax <- length(time)
  
  year <- rep(tstart:tend, each=52*7)[1:tmax]
  week <- rep(rep(1:52, each=7), ceiling(tend-tstart+1))[1:tmax]
  
  S <- I <- R <- C <- Reff <- recruitment <- rep(0, tmax)
  
  S[1] <- round(S0 * pop)
  I[1] <- round(I0 * pop)
  R[1] <- round((1-S0-I0) * pop)
  
  for (i in 2:tmax) {
    R0eff <- R0 * (1 + theta * cos(2*pi*time[i])) * npifun(time[i])
    
    Reff[i] <- R0eff * S[i-1]/pop
    foi <- R0eff * (I[i-1]+omega)/pop * (1-exp(-gamma))
    
    Sout <- rbinom(1, S[i-1], 1-exp(-(foi+mu)))
    StoI <- rbinom(1, Sout, foi/(foi+mu))
    Iout <- rbinom(1, I[i-1], 1-exp(-(gamma+mu)))
    ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
    Rout <- rbinom(1, R[i-1], 1-exp(-(delta+mu)))
    RtoS <- rbinom(1, Rout, delta/(delta+mu))
    
    birth <- rpois(1, mu * pop)
    
    recruitment[i] <- birth + RtoS
    
    S[i] <- S[i-1] - Sout + birth + RtoS
    I[i] <- I[i-1] + StoI - Iout
    R[i] <- R[i-1] + ItoR - Rout
    C[i] <- StoI
  }
  
  cases <- emdbook::rbetabinom(tmax, rho, C, theta=theta_obs)
  
  # plot(time, cases, type="l", xlim=c(2014, 2030))
  
  data.frame(
    time=time,
    year=year,
    week=week,
    recruitment=recruitment,
    S=S,
    I=I,
    C=C,
    cases=cases,
    model="SIRS"
  )
}
