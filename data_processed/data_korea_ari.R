library(readxl)
library(dplyr)
library(tidyr)

korea_raw <- read_xlsx("../data_korea/korea-ari-2024-11-20.xlsx", skip=28) %>%
  setNames(c("year", "week", "total", 
             "Adenovirus", "Bocavirus", "Parainfluenza virus",
             "RSV", "Rhinovirus", "Human metapneumovirus", "Human coronavirus")) %>%
  select(-total)

data_korea_ari <- korea_raw %>%
  filter(!is.na(week)) %>%
  mutate(
    time=((1:n())+22)/52.25+2015
  ) %>%
  gather(key, value, -year, -week, -time) %>%
  mutate(
    year=as.numeric(as.character(year)),
    week=as.numeric(as.character(week)),
    value=as.numeric(gsub(",", "", as.character(value)))
  ) %>%
  rename(cases=value) %>%
  filter(!is.na(week))

write.csv(data_korea_ari, file="data_korea_ari.csv")
