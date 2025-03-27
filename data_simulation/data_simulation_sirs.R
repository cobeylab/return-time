library(dplyr)

source("../R/simulate_sirs_stoch.R")

npifun_default <- function(t) {
  if (t >= 2020.25 & t < 2020.5) {
    npi <- 0.6
  } else if (t >= 2020.5 & t < 2021) {
    npi <- 0.9
  } else if (t >= 2021 & t < 2021.25) {
    npi <- 0.8
  } else {
    npi <- 1
  }
  
  return(npi)
}

simulation_SIRS <- simulate_SIRS()

data_simulation_SIRS <- simulation_SIRS %>%
  filter(year >= 2010) %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases),
    model=unique(model),
    S=S[1],
    recruitment=sum(recruitment)
  ) %>%
  mutate(
    rel_recruitment=recruitment/S
  )

save("simulation_SIRS",
     "data_simulation_SIRS",
     file="data_simulation_sirs.rda")
