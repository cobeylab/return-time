library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
source("../script/script_data_korea.R")
source("../script/script_data_korea_int.R")
source("../R/distfun.R")
source("../R/takens.R")

name <- c("adeno", "boca", "hcov", "hmpv", "noro", "piv", "rhino", "rsv")
realname <- c("Adenovirus", "Bocavirus", "Human coronavirus", "Human metapneumovirus",
              "Norovirus", "Parainfluenza virus", "Rhinovirus", "RSV")

reslist <- vector('list', length(name))

for (i in 1:length(name)) {
  print(i)
  if (name[i] != "noro") {
    truedata <- data_korea_ari_scaled
  } else {
    truedata <- data_korea_int_scaled
  }
  
  tmp <- truedata_filter <- truedata %>%
    filter(key==realname[i]) %>%
    arrange(
      year, week
    )
  
  logcases_pre <- log(filter(truedata_filter, year < 2020)$cases+1)
  
  acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)
  
  tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]
  
  n.fnn <- fnn(logcases_pre, dmax=15, tau=tau, R_tol=10)
  
  d <- which(n.fnn==0)[1]+1
  
  takens_perturb <- takens(log(tmp$cases+1), d=d, tau=tau)
  takens_unperturb <- takens(logcases_pre, d=d, tau=tau)
  
  dist <- sapply(1:nrow(takens_perturb), function(i) {
    dd <- sqrt(colSums((takens_perturb[i,] - t(takens_unperturb))^2))
    
    min(dd[dd>0])
  })
  
  time <- tail(tmp$year+tmp$week/52, -tau*(d-1))
  
  distdata_takens <- data.frame(
    time=time,
    dist=dist
  )
  
  tmp$tau <- tau
  tmp$dist_takens <- c(rep(NA, tau*(d-1)), distdata_takens$dist)
  tmp$d <- d+1
  
  reslist[[i]] <- tmp
}

analysis_korea_takens_acf_fnn <- reslist %>%
  bind_rows

save("analysis_korea_takens_acf_fnn", file="analysis_korea_takens_acf_fnn.rda")
