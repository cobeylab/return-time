library(dplyr)
library(rstan)
source("../script/script_data_korea.R")

data_fit <- data_korea_ari_scaled %>%
  filter(key=="Adenovirus")

model <- stan_model("../stanmodel/growth_tv.stan")

standata <- list(
  N=nrow(data_fit),
  cases=data_fit$cases
)

stanfit_korea_growth_tv_adeno <- sampling(model,
                                         data = standata,
                                         seed=101,
                                         chain=4,
                                         cores=4,
                                         iter=4000)

check_hmc_diagnostics(stanfit_korea_growth_tv_adeno)
get_num_divergent(stanfit_korea_growth_tv_adeno)
get_num_max_treedepth(stanfit_korea_growth_tv_adeno)
get_low_bfmi_chains(stanfit_korea_growth_tv_adeno)
get_bfmi(stanfit_korea_growth_tv_adeno)

save("stanfit_korea_growth_tv_adeno", file="stanfit_korea_growth_tv_adeno.rda")

ss <- summary(stanfit_korea_growth_tv_adeno)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.002193
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 991.3812

plot(data_fit$cases)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("r\\[", rownames(ss$summary)),6], type="l")
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("r\\[", rownames(ss$summary)),8])

ss$summary[grepl("r_sd", rownames(ss$summary)),]
ss$summary[grepl("phi", rownames(ss$summary)),]

plot(density(rnbinom(10000, mu=400, size=172)))
hist(rpois(10000, 400))
