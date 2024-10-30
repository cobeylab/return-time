npifun_default <- function(t) {
  if (t >= 2020 & t < 2020.5) {
    npi <- 0.5
  } else {
    npi <- 1
  }
  
  return(npi)
}

model_rsv_white <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		beta_1 <- b_1 * (1 + theta * cos(2 * pi * (t-phi))) * npifun(t)
		beta_2 <- alpha * beta_1 * npifun(t)
		
		lambda_1 <- beta_1 * (Ip1 + eta * (Is_he1 + Is_ho1 + It1))
		lambda_2 <- beta_2 * (Ip2 + eta * (Is_he2 + Is_ho2 + It2))
		
		dS <- mu - (lambda_1 + lambda_2 + mu) * S + delta * (R1 + R2)
		dIp1 <- lambda_1 * S - (gamma + mu) * Ip1
		dIp2 <- lambda_2 * S - (gamma + mu) * Ip2
		dR1 <- gamma * Ip1 - (sigma_ho * lambda_1 + sigma_he * lambda_2 + delta  + mu) * R1 + gamma * Is_ho1 + delta * R12
		dR2 <- gamma * Ip2 - (sigma_ho * lambda_2 + sigma_he * lambda_1 + delta  + mu) * R2 + gamma * Is_ho2 + delta * R12
		dIs_he1 <- sigma_he * lambda_1 * R2 - (gamma + mu) * Is_he1
		dIs_he2 <- sigma_he * lambda_2 * R1 - (gamma + mu) * Is_he2
		dIs_ho1 <- sigma_ho * lambda_1 * R1 - (gamma + mu) * Is_ho1
		dIs_ho2 <- sigma_ho * lambda_2 * R2 - (gamma + mu) * Is_ho2
		dR12 <- gamma * (Is_he1 + Is_he2 + It1 + It2) - (sigma_he * sigma_ho * (lambda_1 + lambda_2) + 2 * delta + mu) * R12
		dIt1 <- sigma_he * sigma_ho * lambda_1 * R12 - (gamma + mu) * It1
		dIt2 <- sigma_he * sigma_ho * lambda_2 * R12 - (gamma + mu) * It2
		
		prevalence1 <- Ip1 + Is_he1 + Is_ho1 + It1
		prevalence2 <- Ip2 + Is_he2 + Is_ho2 + It2
		
		R01 <- beta_1/(gamma+mu)
		R02 <- alpha * R01
		
		trans_eff1 <- (Ip1 + eta * (Is_he1 + Is_ho1 + It1))/prevalence1
		trans_eff2 <- (Ip2 + eta * (Is_he2 + Is_ho2 + It2))/prevalence2
		
		Seff1 <- S + sigma_ho * R1 + sigma_he * R2 + sigma_ho * sigma_he * R12
		Seff2 <- S + sigma_ho * R2 + sigma_he * R1 + sigma_ho * sigma_he * R12
		
		Reff1 <- R01 * trans_eff1 * Seff1
		Reff2 <- R02 * trans_eff2 * Seff2
		
		list(c(dS, dIp1, dIp2, dR1, dR2, dIs_he1, dIs_he2, dIs_ho1, dIs_ho2, dR12, dIt1, dIt2),
		     beta1=beta_1,
		     beta2=beta_2,
		     prevalence1=prevalence1, prevalence2=prevalence2,
		     Seff1=Seff1,
		     Seff2=Seff2,
		     Reff1=Reff1,
		     Reff2=Reff2)
	})
}

simulate_rsv_white <- function(alpha=0.9159,
                               eta=0.4126,
                               sigma_ho=0.3569,
                               sigma_he=0.8426,
                               gamma=40.56,
                               delta=0.51,
                               b_1=99.51,
                               mu=0.012,
                               phi=0.97,
                               theta=0.347,
                               npifun=npifun_default,
                               yini,
                               tmin=1900,
                               tmax=2030) {
	parms <- c(alpha=alpha,
						 eta=eta,
						 sigma_ho=sigma_ho,
						 sigma_he=sigma_he,
						 gamma=gamma,
						 delta=delta,
						 b_1=b_1,
						 mu=mu,
						 phi=phi,
						 theta=theta,
						 npifun=npifun)
	
	if (missing(yini)) {
		yini <- c(S=1-2e-6, Ip1=1e-6, Ip2=1e-6,
							R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)
	}
	
	times <- seq(tmin, tmax, by=1/365)
	
	out <- as.data.frame(deSolve::ode(yini, times, model_rsv_white, parms))
	
	out$r1 <- c(NA, diff(out$prevalence1))/out$prevalence1
	out$r2 <- c(NA, diff(out$prevalence2))/out$prevalence2
	
	out
}

eigen_rsv_white <- function(alpha=0.9159,
                            eta=0.4126,
                            sigma_ho=0.3569,
                            sigma_he=0.8426,
                            gamma=40.56,
                            delta=0.51,
                            b_1=99.51,
                            mu=0.012) {
  parms <- c(alpha=alpha,
             eta=eta,
             sigma_ho=sigma_ho,
             sigma_he=sigma_he,
             gamma=gamma,
             delta=delta,
             b_1=b_1,
             mu=mu,
             phi=0,
             theta=0,
             npifun=function(x) 1)
  
  yini <- c(S=1-2e-6, Ip1=1e-6, Ip2=1e-6,
            R1=0, R2=0, Is_he1=0, Is_he2=0, Is_ho1=0, Is_ho2=0, R12=0, It1=0, It2=0)
  
  times <- seq(1900, 2000, by=1/365)
  
  out <- as.data.frame(deSolve::ode(yini, times, model_rsv_white, parms))
  
  equi <- unlist(tail(out,1)[2:13])
  
  function_deriv_only <- function(y) {
    model_rsv_white(0, y, parms)[[1]]
  }
  
  jac <- numDeriv::jacobian(function_deriv_only, equi)
  
  ee <- eigen(jac)
  
  ee  
}
