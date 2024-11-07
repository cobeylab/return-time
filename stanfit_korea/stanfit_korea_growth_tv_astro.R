library(dplyr)
library(rstan)
source("../script/script_data_korea_int.R")

data_fit <- data_korea_int_scaled %>%
  filter(key=="Astrovirus")

model <- stan_model("../stanmodel/growth_tv.stan")

standata <- list(
  N=nrow(data_fit),
  cases=data_fit$cases
)

stanfit_korea_growth_tv_astro <- sampling(model,
                                         data = standata,
                                         seed=101,
                                         chain=4,
                                         cores=4,
                                         iter=8000)

check_hmc_diagnostics(stanfit_korea_growth_tv_astro)
get_num_divergent(stanfit_korea_growth_tv_astro)
get_num_max_treedepth(stanfit_korea_growth_tv_astro)
get_low_bfmi_chains(stanfit_korea_growth_tv_astro)
get_bfmi(stanfit_korea_growth_tv_astro)

save("stanfit_korea_growth_tv_astro", file="stanfit_korea_growth_tv_astro.rda")

ss <- summary(stanfit_korea_growth_tv_astro)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.009509
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 410.344

plot(data_fit$cases)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("r\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),8])
