npifun_random_generate <- function(duration=2,
                                   sigma_npi=0.02,
                                   threshold=0.1,
                                   npimin=0.5,
                                   npiend=1,
                                   seed) {
  if (!missing(seed)) set.seed(seed)
    
  Nnpi <- floor(duration * 364)
  
  npi <- rev(1 + cumsum(rnorm(Nnpi, 0, sigma_npi)))
  
  while (mean(npi > 1) > threshold) {
    npi <- rev(1 + cumsum(rnorm(Nnpi, 0, sigma_npi)))
  }
  
  npi[npi < 0] <- 0
  npi[npi > 1] <- 1

  npimin_emp <- min(npi)
  
  npi <- (npi-npimin_emp)/(1-npimin_emp) * (npiend-npimin) + npimin
  
  time <- 2020+1:floor(duration * 364)/364-1/364
  
  npifun_random <- function(t) {
    if (t < 2020) {
      return(1)
    } else if (t >= 2020+floor(duration * 364)/364) {
      return(npiend)
    } else {
      return(npi[match(round(t, 3), round(time,3))])
    }  
  }
  
  npifun_random
}
