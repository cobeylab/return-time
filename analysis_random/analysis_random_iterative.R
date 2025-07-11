library(tidyr)
library(dplyr)
source("../R/simulate_sirs.R")
source("../R/simulate_sirs_stoch.R")
source("../R/npi_random.R")
source("../R/takens.R")
source("../R/lm_iterative.R")

nsim <- 500

reslist <- vector('list', nsim)

for (i in 1:nsim) {
  set.seed(i)
  print(i)
  extinction <- TRUE
  
  while (extinction) {
    duration <- runif(1, 0.5, 3.5)
    npimin <- runif(1, 0.5, 0.7)
    npiend <- runif(1, 0.8, 1)
    R0 <- runif(1, 1.5, 3)
    immunity <- runif(1, 0.5, 2)
    
    npifun_random <- npifun_random_generate(duration=duration,
                                            npimin=npimin,
                                            npiend=npiend,
                                            seed=i)
    
    npidata <- data.frame(
      time=seq(2015, 2030, by=1/52/7),
      npi=sapply(seq(2015, 2030, by=1/52/7), npifun_random)
    )
    
    ss <- simulate_SIRS_stoch(
      R0=R0,
      S0=1/R0,
      theta=0.1,
      delta=1/52/7/immunity,
      tstart=1990,
      tend=2025,
      npifun=npifun_random,
      seed=i)
    
    extinction <- all(tail(ss$I, 10)==0)
  }
  
  ee <- eigen_sirs(b_1=R0*(365/7+1/50),
                   mu=1/50,
                   gamma=365/7,
                   delta=1/immunity)
  
  ss_weekly <- ss %>%
    group_by(year, week) %>%
    summarize(
      cases=sum(cases)
    ) %>%
    mutate(
      time=year+week/52
    ) %>%
    filter(
      year >= 2014, year < 2025
    ) %>%
    arrange(year, week)
  
  # plot(ss_weekly$cases, type="l", log="y")
  
  logcases_pre <- log(filter(ss_weekly, year < 2020)$cases+1)
  
  acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)
  
  tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]
  
  n.fnn <- fnn(logcases_pre, dmax=25, tau=tau, R_tol=10)
  
  d <- which(n.fnn==0)[1]+1
  
  takens_perturb <- takens(log(ss_weekly$cases+1), d=d, tau=tau)
  takens_unperturb <- takens(logcases_pre, d=d, tau=tau)
  
  takens_data <- as.data.frame(takens_perturb)
  takens_data$time <- tail(ss_weekly$year+ss_weekly$week/52, -tau*(d-1))
  
  dist <- sapply(1:nrow(takens_perturb), function(i) {
    dd <- sqrt(colSums((takens_perturb[i,] - t(takens_unperturb))^2))
    
    min(dd[dd>0])
  })
  
  # plot(dist, log="y")
  
  time <- tail(ss_weekly$year+ss_weekly$week/52, -tau*(d-1))
  
  distdata_takens <- data.frame(
    time=time,
    dist=dist
  )
  
  lfit <- lm_iterative(dist, time, iter=1)
  
  if(lfit[[2]] < 0) {
    reslist[[i]] <- data.frame(
      duration=duration,
      npimin=npimin,
      npiend=npiend,
      R0=R0,
      immunity=immunity,
      resilience_true=ee,
      est=lfit$resilience[1],
      lwr=lfit$resilience_lwr[1],
      upr=lfit$resilience_upr[1],
      reldist=mean(tail(dist, 52))/mean(dist[time < 2020])
    )
  } else {
    if (mean(tail(dist, 52))/mean(dist[time < 2020]) > exp(2)) {
      # plot(time, dist, log="y")
      # curve(exp(-ee*(x-2024)), add=TRUE)
      
      lfit <- lm_iterative(dist, time, iter=1)
      
      reslist[[i]] <- data.frame(
        duration=duration,
        npimin=npimin,
        npiend=npiend,
        R0=R0,
        immunity=immunity,
        resilience_true=ee,
        est=lfit$resilience[1],
        lwr=lfit$resilience_lwr[1],
        upr=lfit$resilience_upr[1],
        reldist=mean(tail(dist, 52))/mean(dist[time < 2020])
      )
    } else {
      lfit <- lm_iterative(dist, time, iter=20)
      
      reslist[[i]] <- data.frame(
        duration=duration,
        npimin=npimin,
        npiend=npiend,
        R0=R0,
        immunity=immunity,
        resilience_true=ee,
        est=lfit$resilience[20],
        lwr=lfit$resilience_lwr[20],
        upr=lfit$resilience_upr[20],
        reldist=mean(tail(dist, 52))/mean(dist[time < 2020])
      )
    }
  }  
}

analysis_random_iterative <- reslist %>%
  bind_rows

save("analysis_random_iterative",
     file="analysis_random_iterative.rda")

plot(analysis_random_iterative$resilience_true, 
     analysis_random_iterative$est)
