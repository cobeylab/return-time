eigen_seir <- function(R0=17,
                       mu=1/50,
                       sigma=1/8*365,
                       gamma=1/5*365) {
  beta <- R0/(sigma/(sigma+mu)) * (gamma + mu)
  
  # mu - beta * S * I - mu * S
  # beta * S * I - (sigma+mu) * E
  # sigma * E - (gamma+mu) * I
  # R = 1 - S - E - I
  
  # beta * S * I = (sigma+mu) * E
  # E = (gamma+mu) * I/sigma
  # (sigma+mu) * (gamma+mu)/sigma
  # S =  (sigma+mu) * (gamma+mu)/sigma/beta = 1/R0
  
  # mu - mu * S = beta * S * I
  
  S <- 1/R0
  
  I <- (mu - mu * S)/(beta * S)
  
  E <- (gamma+mu)*I/sigma
  
  mat <- matrix(
    c(-beta*I-mu, 0, -beta*S,
      beta*I,-(sigma+mu),beta*S,
      0,sigma,-(gamma+mu)),
    3, 3,
    byrow=TRUE
  )
  
  ee <- eigen(mat)
  
  max(Re(ee$values))
}
