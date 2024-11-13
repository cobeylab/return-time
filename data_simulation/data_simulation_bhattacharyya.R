library(dplyr)

source("../R/simulate_bhattacharyya_stoch.R")

simulation_bhattacharyya <- simulate_bhattacharyya()

data_simulation_bhattacharyya1 <- simulation_bhattacharyya %>%
  filter(year >= 2010) %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases1),
    S=S1[1],
    recruitment=sum(recruitment1)
  ) %>%
  mutate(
    rel_recruitment=recruitment/S
  )

data_simulation_bhattacharyya2 <- simulation_bhattacharyya %>%
  filter(year >= 2010) %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases2),
    S=S2[1],
    recruitment=sum(recruitment2)
  ) %>%
  mutate(
    rel_recruitment=recruitment/S
  )

data_simulation_bhattacharyya_both <- simulation_bhattacharyya %>%
  filter(year >= 2010) %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases1+cases2)
  )

save(simulation_bhattacharyya,
     data_simulation_bhattacharyya1,
     data_simulation_bhattacharyya2,
     data_simulation_bhattacharyya_both,
     file="data_simulation_bhattacharyya.rda")
