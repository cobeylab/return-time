library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
source("../script/script_data_korea.R")
source("../script/script_data_korea_int.R")
source("../R/distfun.R")

name <- c("adeno", "boca", "hcov", "hmpv", "noro", "piv", "rhino", "rsv")
realname <- c("Adenovirus", "Bocavirus", "Human coronavirus", "Human metapneumovirus",
              "Norovirus", "Parainfluenza virus", "Rhinovirus", "RSV")

nsamp <- 200

estlist <- vector('list', length(name))

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
    
    distdata <- distfun(
      mat_perturb = as.matrix(tmp[,c("r", "logC")]),
      mat_unperturb = as.matrix(tmp[,c("meanr", "meanlogC")]),
      time=tmp$year+tmp$week/52,
      tbreak=2020,
      out="all"
    )
    
    maxt <- distdata$time[which.max(distdata$dist)]
    
    lfit <- lm(log(dist)~time, data=filter(distdata, time >= maxt))
    
    tmp$maxt <- maxt
    tmp$dist <- distdata$dist
    tmp$return <- -1/coef(lfit)[[2]]
    tmp$key <- realname[i]
    tmp$id <- j
    tmp$whichsamp <- whichsamp[j]
    
    reslist[[j]] <- tmp
  }
  
  resdf <- reslist %>%
    bind_rows
  
  estlist[[i]] <- resdf
}

analysis_korea_growth <- estlist %>%
  bind_rows

save("analysis_korea_growth", file="analysis_korea_growth.rda")
