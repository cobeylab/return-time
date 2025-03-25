library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
source("../script/script_data_nrevss.R")
source("../R/distfun.R")
source("../R/takens.R")

data_nrevss_resp_comb_proxy2 <- data_nrevss_resp_comb_proxy %>%
  filter(type %in% c("Adenovirus", "PIV", "CoV", "Human metapneumovirus", "RSV",
                     "Rhinovirus"),
         year >= 2014)

unique(data_nrevss_resp_comb_proxy2$type)

realname <- c("Rhinovirus", "Human metapneumovirus", "PIV", "Adenovirus", 
              "CoV", "RSV")
reslist <- vector('list', length(realname))

Rvec <- c(5, 10, 15)

for (i in 1:length(realname)) {
  print(i)
  truedata <- data_nrevss_resp_comb_proxy2
  
  tmp <- truedata_filter <- truedata %>%
    filter(type==realname[i]) %>%
    arrange(
      year, week
    )
  
  logcases_pre <- log((filter(truedata_filter, year < 2020)$positive+1)/filter(truedata_filter, year < 2020)$tests*filter(truedata_filter, year < 2020)$ili)
  
  acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)
  
  tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]
  
  Rlist <- vector('list', length(Rvec))
  
  for (j in 1:length(Rvec)) {
    n.fnn <- fnn(logcases_pre, dmax=15, tau=tau, R_tol=Rvec[j])
    
    d <- which(n.fnn==0)[1]+1
    
    takens_perturb <- takens(log((tmp$positive+1)/tmp$tests*tmp$ili), d=d, tau=tau)
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
    
    tmp2 <- tmp
    
    tmp2$tau <- tau
    tmp2$dist_takens <- c(rep(NA, tau*(d-1)), distdata_takens$dist)
    tmp2$d <- d+1
    tmp2$R <- Rvec[j]
    
    Rlist[[j]] <- tmp2
  }
  
  reslist[[i]] <- Rlist %>%
    bind_rows
}

analysis_nrevss_takens_acf_fnn_R <- reslist %>%
  bind_rows

save("analysis_nrevss_takens_acf_fnn_R", file="analysis_nrevss_takens_acf_fnn_R.rda")
