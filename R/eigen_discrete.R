eigen_discrete <- function(beta,
         omega,
         N,
         mu,
         gamma) {
  FS <- expression(S + mu * N - (1-exp(-(beta*(I+omega)/N+mu))) * S)
  FI <- expression(I + beta*(I+omega)/N/(beta*(I+omega)/N + mu) * (1-exp(-(beta*(I+omega)/N+mu))) * S
                   - (1-exp(-(gamma+mu)))*I)
  
  R0 <- beta/(1-exp(-(gamma+mu)))
  Sstar <- 1/R0 * N
  Istar <- mu * (1-Sstar/N)/(beta*Sstar/N) * N
  
  oo <- optim(c(x=log(Sstar), y=log(Istar)),
        function(par) {
          x <- par[1]
          y <- par[2]
          (eval(FS, list(S=exp(x), I=exp(y)))-exp(x))^2 +
            (eval(FI, list(S=exp(x), I=exp(y)))-exp(y))^2
          })
  
  Sstar_disc <- exp(oo$par[1])
  Istar_disc <- exp(oo$par[2])
  
  Amat <- matrix(c(
    eval(D(FS, "S"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FS, "I"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FI, "S"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FI, "I"), list(S=Sstar_disc, I=Istar_disc))
  ), 2, 2)
  
  ee <- eigen(Amat)
  
  max(Re(ee$values))
}

## beta * S * I - gamma * I
## delta * (N - S - I) + mu * N - beta * S * I/N - mu * S = 0
## delta * N - delta * S - delta * I + mu * N - beta * S * I/N - mu * S = 0
## (delta + mu) * (N-S) = (delta + beta * S/N) * I

eigen_discrete_sirs <- function(beta,
                                omega,
                                N,
                                mu,
                                gamma,
                                delta) {
  FS <- expression(S + mu * N - (1-exp(-(beta*(I+omega)/N+mu))) * S 
                   + delta/(delta + mu) * (1-exp(-(delta+mu))) * (N-S-I))
  FI <- expression(I + beta*(I+omega)/N/(beta*(I+omega)/N + mu) * (1-exp(-(beta*(I+omega)/N+mu))) * S
                   - (1-exp(-(gamma+mu)))*I)
  
  R0 <- beta/(1-exp(-(gamma+mu)))
  Sstar <- 1/R0 * N
  Istar <- (delta + mu) * (N-Sstar)/(delta + beta * Sstar/N)
  
  oo <- optim(c(x=log(Sstar), y=log(Istar)),
              function(par) {
                x <- par[1]
                y <- par[2]
                (eval(FS, list(S=exp(x), I=exp(y)))-exp(x))^2 +
                  (eval(FI, list(S=exp(x), I=exp(y)))-exp(y))^2
              },
              control=list(maxit=1e5))
  
  oo2 <- optim(c(x=oo$par[1], y=oo$par[2]),
               function(par) {
                 x <- par[1]
                 y <- par[2]
                 (eval(FS, list(S=exp(x), I=exp(y)))-exp(x))^2 +
                   (eval(FI, list(S=exp(x), I=exp(y)))-exp(y))^2
               },
               control=list(maxit=1e5))
  
  Sstar_disc <- exp(oo2$par[1])
  Istar_disc <- exp(oo2$par[2])
  
  Amat <- matrix(c(
    eval(D(FS, "S"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FS, "I"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FI, "S"), list(S=Sstar_disc, I=Istar_disc)),
    eval(D(FI, "I"), list(S=Sstar_disc, I=Istar_disc))
  ), 2, 2)
  
  ee <- eigen(Amat)
  
  max(Re(ee$values))
}
