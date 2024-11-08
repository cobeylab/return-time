library(dplyr)
library(rstan)
source("../script/script_data_hongkong_resp.R")

data_fit <- data_hongkong_resp_scaled %>%
  filter(key=="rvev")

model <- stan_model("../stanmodel/growth_tv_poisson.stan")

standata <- list(
  N=nrow(data_fit),
  cases=data_fit$cases
)

stanfit_hongkong_growth_tv_rvev <- sampling(model,
                                         data = standata,
                                         seed=101,
                                         chain=4,
                                         cores=4,
                                         iter=4000,
                                         control=list(
                                           max_treedepth=12
                                         ))

check_hmc_diagnostics(stanfit_hongkong_growth_tv_rvev)
get_num_divergent(stanfit_hongkong_growth_tv_rvev)
get_num_max_treedepth(stanfit_hongkong_growth_tv_rvev)
get_low_bfmi_chains(stanfit_hongkong_growth_tv_rvev)
get_bfmi(stanfit_hongkong_growth_tv_rvev)

save("stanfit_hongkong_growth_tv_rvev", file="stanfit_hongkong_growth_tv_rvev.rda")

ss <- summary(stanfit_hongkong_growth_tv_rvev)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.003
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 774

plot(data_fit$cases)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("r\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),8])
