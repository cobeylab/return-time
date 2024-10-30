data_canada_resp <- read.csv("../data_processed/data_canada_resp.csv") %>%
  mutate(
    date=as.Date(date)
  ) %>%
  dplyr::mutate(
    week=ifelse(week==53, 52, week),
    year=lubridate::epiyear(date)
  ) %>%
  dplyr::group_by(province, type, year, week) %>%
  dplyr::summarize(
    tests=round(mean(tests)),
    positive=round(mean(positive))
  )

data_canada_resp_scaled <- data_canada_resp %>% 
  filter(province=="Canada", type != "RV") %>%
  group_by(type) %>%
  arrange(type, year, week) %>%
  mutate(
    tests_ra=zoo::rollapply(tests, 52*2, mean, partial=TRUE),
    scale=tail(tests_ra,1)/tests_ra,
    cases=round(tests*positive/100*scale)
  ) %>%
  rename(
    key=type
  )
