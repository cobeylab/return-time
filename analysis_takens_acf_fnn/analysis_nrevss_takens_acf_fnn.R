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
  
  n.fnn <- fnn(logcases_pre, dmax=25, tau=tau, R_tol=10)
  
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
  
  maxt_takens <- distdata_takens$time[which.max(distdata_takens$dist)]
  
  lfit_takens <- lm(log(dist)~time, data=filter(distdata_takens, time >= maxt_takens))
  
  # pp_takens <- predict(lfit_takens, newdata=data.frame(time=seq(maxt_takens, 2040, by=0.1)))
  
  tmp$tau <- tau
  tmp$maxt_takens <- maxt_takens
  tmp$dist_takens <- c(rep(NA, tau*(d-1)), distdata_takens$dist)
  tmp$d <- d+1
  tmp$intercept_takens <- coef(lfit_takens)[[1]]
  tmp$resilience <- coef(lfit_takens)[[2]]
  tmp$resilience_lwr <- confint(lfit_takens)[2,1]
  tmp$resilience_upr <- confint(lfit_takens)[2,2]
  
  reslist[[i]] <- tmp
}

analysis_nrevss_takens_acf_fnn <- reslist %>%
  bind_rows

save("analysis_nrevss_takens_acf_fnn", file="analysis_nrevss_takens_acf_fnn.rda")
