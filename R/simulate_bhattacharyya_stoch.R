npifun_default <- function(t) {
  if (t >= 2020.25 & t < 2020.5) {
    npi <- 0.6
  } else if (t >= 2020.5 & t < 2021) {
    npi <- 0.9
  } else if (t >= 2021 & t < 2021.25) {
    npi <- 0.8
  } else {
    npi <- 1
  }
  
  return(npi)
}

simulate_bhattacharyya <- function(b01=3.4*52/364,
                                   b02=3.9*52/364,
                                   theta1=0.4,
                                   theta2=0.3,
                                   phi1=0.005*7/364,
                                   phi2=4.99*7/364,
                                   epsilon12=0.92,
                                   epsilon21=0.45,
                                   gamma1=1/10,
                                   gamma2=1/10,
                                   rho1=1/364,
                                   rho2=1/364,
                                   q1=0.5,
                                   q2=0.5,
                                   mu=1/70/364,
                                   pop=1e8,
                                   I0=1e-6,
                                   rho=0.02,
                                   theta_obs=500,
                                   npifun=npifun_default,
                                   tstart=1900,
                                   tend=2030,
                                   seed=101) {
  set.seed(seed)
  time <- seq(tstart, tend, by=1/52/7)
  tmax <- length(time)
  
  year <- rep(tstart:tend, each=52*7)[1:tmax]
  week <- rep(rep(1:52, each=7), ceiling(tend-tstart+1))[1:tmax]
  
  S <- I1 <- I2 <- R1 <- R2 <- J1 <- J2 <- R <- C1 <- C2 <- recruitment1 <- recruitment2 <- rep(0, tmax)
  
  S[1] <- round((1-2*I0)*pop)
  I1[1] <- round(I0*pop)
  I2[1] <- round(I0*pop)
  
  for (i in 2:tmax) {
    beta1 <- b01 * (1 + theta1 * sin(2*pi*(time[i]-phi1))) * q1 * npifun(time[i])
    beta2 <- b02 * (1 + theta2 * sin(2*pi*(time[i]-phi2))) * q2 * npifun(time[i])
    
    foi1 <- beta1 * (I1[i-1] + J1[i-1])/pop
    foi2 <- beta2 * (I2[i-1] + J2[i-1])/pop
    
    birth <- rpois(1, mu * pop)
    
    Sout <- rbinom(1, S[i-1], (1-exp(-(foi1+foi2+mu))))
    StoI1 <- rbinom(1, Sout, foi1/(foi1+foi2+mu))
    StoI2 <- rbinom(1, Sout, foi2/(foi1+foi2+mu))
    I1out <- rbinom(1, I1[i-1], (1-exp(-(gamma1+mu))))
    I1toR1 <- rbinom(1, I1out, gamma1/(gamma1+mu))
    I2out <- rbinom(1, I2[i-1], (1-exp(-(gamma2+mu))))
    I2toR2 <- rbinom(1, I2out, gamma2/(gamma2+mu))
    R1out <- rbinom(1, R1[i-1], (1-exp(-(rho1+mu+epsilon21*foi2))))
    R1toS <- rbinom(1, R1out, rho1/(rho1+mu+epsilon21*foi2))
    R1toJ2 <- rbinom(1, R1out, epsilon21*foi2/(rho1+mu+epsilon21*foi2))
    R2out <- rbinom(1, R2[i-1], (1-exp(-(rho2+mu+epsilon12*foi1))))
    R2toS <- rbinom(1, R2out, rho2/(rho2+mu+epsilon12*foi1))
    R2toJ1 <- rbinom(1, R2out, epsilon12*foi1/(rho2+mu+epsilon12*foi1))
    J1out <- rbinom(1, J1[i-1], (1-exp(-(gamma1+mu))))
    J1toR <- rbinom(1, J1out, gamma1/(gamma1+mu))
    J2out <- rbinom(1, J2[i-1], (1-exp(-(gamma2+mu))))
    J2toR <- rbinom(1, J2out, gamma2/(gamma2+mu))
    Rout <- rbinom(1, R[i-1], (1-exp(-(rho1+rho2+mu))))
    RtoR1 <- rbinom(1, Rout, rho2/(rho1+rho2+mu))
    RtoR2 <- rbinom(1, Rout, rho1/(rho1+rho2+mu))
    
    S[i] <- S[i-1] - Sout + birth + R1toS + R2toS
    I1[i] <- I1[i-1] + StoI1 - I1out
    I2[i] <- I2[i-1] + StoI2 - I2out
    R1[i] <- R1[i-1] + I1toR1 - R1out + RtoR1
    R2[i] <- R2[i-1] + I2toR2 - R2out + RtoR2
    J1[i] <- J1[i-1] + R2toJ1 - J1out
    J2[i] <- J2[i-1] + R1toJ2 - J2out
    R[i] <- R[i-1] + J1toR + J2toR - Rout
    
    C1[i] <- StoI1 + R2toJ1
    C2[i] <- StoI2 + R1toJ2
    
    recruitment1[i] <- birth + R1toS + epsilon21*RtoR2
    recruitment2[i] <- birth + R2toS + epsilon12*RtoR1
  }
  
  cases1 <- emdbook::rbetabinom(tmax, rho, C1, theta=theta_obs)
  cases2 <- emdbook::rbetabinom(tmax, rho, C2, theta=theta_obs)
  
  data.frame(
    time=time,
    year=year,
    week=week,
    recruitment1=recruitment1,
    recruitment2=recruitment2,
    S=S,
    I1=I1,
    I2=I2,
    R1=R1,
    R2=R2,
    J1=J1,
    J2=J2,
    R=R,
    cases1=cases1,
    cases2=cases2,
    S1=S+epsilon21*R2,
    S2=S+epsilon12*R1
  )
}
