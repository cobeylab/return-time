library(dplyr)
library(rstan)
load("../data_simulation/data_simulation_bhattacharyya.rda")

data_all <- data_simulation_bhattacharyya1 %>%
  filter(year >= 2014)

data_fit <- data_all %>%
  filter(year < 2024)

model <- stan_model("../stanmodel/sirs_npi.stan")

npiwidth <- 1
npistart <- which(data_fit$year==2020 & data_fit$week==13)
npiend <- nrow(data_fit)
npiwhich <- ceiling((1:(npiend-npistart+1))/npiwidth)
Nnpi <- max(npiwhich)

standata <- list(
  N=nrow(data_fit),
  Npred=nrow(data_all)-nrow(data_fit),
  Nnpi=Nnpi,
  npistart=npistart,
  npiend=npiend,
  npiwhich=npiwhich,
  week=data_all$week,
  cases=data_fit$cases,
  mu=1/70/52,
  pop=1e8,
  gamma=1
)

stanfit_sirs_bhattacharyya1 <- sampling(model,
                             data = standata,
                             seed=103,
                             chain=4,
                             cores=4,
                             iter=4000,
                             control=list(
                               adapt_delta=0.9,
                               max_treedepth=15
                             ))

check_hmc_diagnostics(stanfit_sirs_bhattacharyya1)
get_num_divergent(stanfit_sirs_bhattacharyya1)
get_num_max_treedepth(stanfit_sirs_bhattacharyya1)
get_low_bfmi_chains(stanfit_sirs_bhattacharyya1)
get_bfmi(stanfit_sirs_bhattacharyya1)

save("stanfit_sirs_bhattacharyya1", file="stanfit_sirs_bhattacharyya1.rda")

ss <- summary(stanfit_sirs_bhattacharyya1)

plot(data_fit$rel_recruitment)
lines((1/50/52*1e8)/ss$summary[grepl("S\\[", rownames(ss$summary)),6], col=2)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.012
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 256.7

plot(data_all$cases)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("npieff\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("beta\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),4], type="l")
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),8], type="l")

plot(ss$summary[grepl("S\\[", rownames(ss$summary)),6])
