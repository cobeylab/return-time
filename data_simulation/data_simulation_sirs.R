library(dplyr)

source("../R/simulate_sirs_stoch.R")

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
