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

nsamp <- 200

estlist <- vector('list', length(name))
predlist <- vector('list', length(name))

set.seed(101)
for (i in 1:length(name)) {
  print(i)
  load(paste0("../stanfit_korea/stanfit_korea_growth_tv_", name[i], ".rda"))
  
  if (name[i] != "noro") {
    truedata <- data_korea_ari_scaled
  } else {
    truedata <- data_korea_int_scaled
  }
  
  truedata_filter <- truedata %>%
    filter(key==realname[i])
  
  ss <- get(paste0("stanfit_korea_growth_tv_", name[i]))
  
  dd <- as_draws_df(ss)
  
  npost <- nrow(dd)
  
  whichsamp <- sample(npost, nsamp)
  
  reslist <- vector('list', nsamp)
  reslist_pred  <- vector('list', nsamp)
  
  for (j in 1:nsamp) {
    dd_samp <- dd[whichsamp[j],]
    
    r <- c(NA, unname(unlist(dd_samp[,grepl("r\\[", colnames(dd_samp))])))
    
    C <- unname(unlist(dd_samp[,grepl("C\\[", colnames(dd_samp))]))
    
    if (name[i] == "noro") {
      tmp <- data.frame(
        year=truedata_filter$year,
        week=truedata_filter$week,
        r=r,
        logC=log(C)
      ) %>%
        mutate(even=year %% 2) %>%
        group_by(even, week) %>%
        mutate(
          meanr=mean(r[year < 2020], na.rm=TRUE),
          meanlogC=mean(logC[year < 2020], na.rm=TRUE)
        )
      
    } else {
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
    }
    
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
        key=realname[i],
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
    
    tmp$key <- realname[i]
    tmp$id <- j
    tmp$whichsamp <- whichsamp[j]
    
    reslist[[j]] <- tmp
  }
  
  resdf <- reslist %>%
    bind_rows
  
  estlist[[i]] <- resdf
  
  predlist[[i]] <- reslist_pred %>%
    bind_rows
}

analysis_korea_growth <- estlist %>%
  bind_rows

analysis_korea_growth_pred <- predlist %>%
  bind_rows

save("analysis_korea_growth", "analysis_korea_growth_pred", file="analysis_korea_growth.rda")
