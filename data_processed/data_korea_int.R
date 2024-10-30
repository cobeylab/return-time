library(readxl)
library(dplyr)
library(tidyr)

korea_raw <- read_xlsx("../data_korea/korea-intestinal-2024-09-04.xlsx", skip=7) %>%
  setNames(c("year", "week", "total",
             "Salmonella", "Vibrio", "ETEC", "EIEC", "EPEC",
             "Campylobacter", "Clostridium perfringens",
             "Staphylococcus aureus", "Bacillus cereus",
             "Yersinia enterocolitica", "Listeria monocytogenes",
             "Group A Rotavirus", "Astrovirus", "Adenovirus",
             "Norovirus", "Sapovirus", "Entamoeba histolytica", 
             "Giardia duodenalis", "Cryptosporidium parvum",
             "Cyclospora")) %>%
  select(-total)

data_korea_int <- korea_raw %>%
  filter(
    !is.na(year),
    Salmonella != "-",
    Salmonella != "집계 중",
    !is.na(Salmonella)
  ) %>%
  mutate(
    year=as.numeric(year),
    week=as.numeric(week)
  ) %>%
  gather(key, value, -year, -week) %>%
  mutate(
    value=as.numeric(value)
  )

write.csv(data_korea_int, file="data_korea_int.csv")
