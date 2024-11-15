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
  
  truedata_filter <- truedata %>%
    filter(key==realname[i])
  
  logcases_pre <- log(filter(truedata_filter, year < 2020)$cases+1)
  
  acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)
  
  tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]
  
  tmp <- truedata_filter %>%
    group_by(week) %>%
    mutate(
      meanlogcase=mean(log(cases+1)[year < 2020])
    )
  
  takens_perturb <- takens(log(tmp$cases+1), d=2, tau=tau)
  takens_unperturb <- takens(tmp$meanlogcase, d=2, tau=tau)
  
  distdata_takens <- distfun(
    mat_perturb = takens_perturb,
    mat_unperturb = takens_unperturb,
    time=tail(tmp$year+tmp$week/52, -tau),
    tbreak=2020,
    out="all"
  )
  
  maxt_takens <- distdata_takens$time[which.max(distdata_takens$dist)]
  
  lfit_takens <- lm(log(dist)~time, data=filter(distdata_takens, time >= maxt_takens))
  
  # pp_takens <- predict(lfit_takens, newdata=data.frame(time=seq(maxt_takens, 2040, by=0.1)))
  
  tmp$tau <- tau
  tmp$maxt_takens <- maxt_takens
  tmp$dist_takens <- c(rep(NA, tau), distdata_takens$dist)
  tmp$intercept_takens <- coef(lfit_takens)[[1]]
  tmp$resilience <- coef(lfit_takens)[[2]]
  tmp$resilience_lwr <- confint(lfit_takens)[2,1]
  tmp$resilience_upr <- confint(lfit_takens)[2,2]
  
  tmp$key <- realname[i]
  
  reslist[[i]] <- tmp
}

analysis_korea_takens_acf_2d <- reslist %>%
  bind_rows

save("analysis_korea_takens_acf_2d", file="analysis_korea_takens_acf_2d.rda")
