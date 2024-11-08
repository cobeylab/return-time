library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
source("../script/script_data_hongkong_piv.R")
source("../R/distfun.R")
source("../R/takens.R")
load("../stanfit_hongkong/stanfit_hongkong_growth_tv_piv.rda")

nsamp <- 200

truedata <- data_hongkong_piv_scaled %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases)
  )

truedata_filter <- truedata

ss <- stanfit_hongkong_growth_tv_piv

dd <- as_draws_df(ss)

npost <- nrow(dd)

whichsamp <- sample(npost, nsamp)

reslist <- vector('list', nsamp)
reslist_pred  <- vector('list', nsamp)

for (j in 1:nsamp) {
  dd_samp <- dd[whichsamp[j],]
  
  r <- c(NA, unname(unlist(dd_samp[,grepl("r\\[", colnames(dd_samp))])))
  
  C <- unname(unlist(dd_samp[,grepl("C\\[", colnames(dd_samp))]))
  
  tmp <- data.frame(
    year=truedata_filter$year,
    week=truedata_filter$week,
    r=r,
    logC=log(C)
  ) %>%
    group_by(week) %>%
    mutate(
      meanr=mean(r[year < 2020], na.rm=TRUE),
      meanlogC=mean(logC[year < 2020], na.rm=TRUE)
    )
  
  distdata_rC <- distfun(
    mat_perturb = as.matrix(tmp[,c("r", "logC")]),
    mat_unperturb = as.matrix(tmp[,c("meanr", "meanlogC")]),
    time=tmp$year+tmp$week/52,
    tbreak=2020,
    out="all"
  )
  
  maxt_rC <- distdata_rC$time[which.max(distdata_rC$dist)]
  
  lfit_rC <- lm(log(dist)~time, data=filter(distdata_rC, time >= maxt_rC))
  
  pp_rC <- predict(lfit_rC, newdata=data.frame(time=seq(maxt_rC, 2040, by=0.1)))
  
  takens_perturb <- takens(tmp$logC, d=2, tau=4)
  takens_unperturb <- takens(tmp$meanlogC, d=2, tau=4)
  
  distdata_takens <- distfun(
    mat_perturb = takens_perturb,
    mat_unperturb = takens_unperturb,
    time=tail(tmp$year+tmp$week/52, -4),
    tbreak=2020,
    out="all"
  )
  
  maxt_takens <- distdata_takens$time[which.max(distdata_takens$dist)]
  
  lfit_takens <- lm(log(dist)~time, data=filter(distdata_takens, time >= maxt_takens))
  
  pp_takens <- predict(lfit_takens, newdata=data.frame(time=seq(maxt_takens, 2040, by=0.1)))
  
  reslist_pred[[j]] <- bind_rows(
    data.frame(
      time=seq(maxt_rC, 2040, by=0.1),
      pred=exp(pp_rC),
      method="growth rate"
    ),
    data.frame(
      time=seq(maxt_takens, 2040, by=0.1),
      pred=exp(pp_takens),
      method="Takens' theorem"
    )
  ) %>%
    mutate(
      key="PIV",
      id=j
    )
  
  tmp$maxt_rC <- maxt_rC
  tmp$dist_rC <- distdata_rC$dist
  tmp$intercept_rC <-  coef(lfit_rC)[[1]]
  tmp$return_rC <- -1/coef(lfit_rC)[[2]]
  
  tmp$maxt_takens <- maxt_takens
  tmp$dist_takens <- c(rep(0, 4), distdata_takens$dist)
  tmp$intercept_takens <- coef(lfit_takens)[[1]]
  tmp$return_takens <- -1/coef(lfit_takens)[[2]]
  
  tmp$key <- "PIV"
  tmp$id <- j
  tmp$whichsamp <- whichsamp[j]
  
  reslist[[j]] <- tmp
}

analysis_hongkong_piv_growth <- reslist %>%
  bind_rows

analysis_hongkong_piv_growth_pred <- reslist_pred %>%
  bind_rows

save("analysis_hongkong_piv_growth", "analysis_hongkong_piv_growth_pred", file="analysis_hongkong_piv_growth.rda")
