data_korea_int <- read.csv("../data_processed/data_korea_int.csv") %>%
  dplyr::mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  dplyr::group_by(key, year, week) %>%
  dplyr::summarize(
    cases=round(mean(value))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    key, year, week
  ) %>%
  dplyr::group_by(key) %>%
  dplyr::mutate(
    time=1:n()
  )

data_korea_scale_int <- data_korea_int %>%
  group_by(year, week) %>%
  summarize(
    total=sum(cases)
  )

roll_total_int <- zoo::rollapply(filter(data_korea_scale_int, year<2020)$total,
                             width=2*52, mean, partial=TRUE)

scale_korea_int <- c(
  roll_total_int/tail(roll_total_int, 1),
  rep(1, nrow(data_korea_scale_int)-length(roll_total_int))
)

data_korea_int_scaled <- data_korea_int %>%
  group_by(key) %>%
  arrange(key, year, week) %>%
  mutate(
    cases=round(cases/scale_korea_int)
  )
