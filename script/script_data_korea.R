data_korea_ari <- read.csv("../data_processed/data_korea_ari.csv") %>%
  dplyr::filter(year+week/52 <= 2024 + 26/52) %>%
  dplyr::mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  dplyr::group_by(key, year, week) %>%
  dplyr::summarize(
    cases=round(mean(cases))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    key, year, week
  ) %>%
  dplyr::group_by(key) %>%
  dplyr::mutate(
    time=1:n()
  )

data_korea_scale <- data_korea_ari %>%
  group_by(year, week) %>%
  summarize(
    total=sum(cases)
  )

roll_total <- zoo::rollapply(filter(data_korea_scale, year < 2020)$total,
                             width=2*52, mean, partial=TRUE)

scale_korea <- c(
  roll_total/tail(roll_total, 1),
  rep(1, nrow(data_korea_scale)-length(roll_total))
)

data_korea_ari_scaled <- data_korea_ari %>%
  group_by(key) %>%
  arrange(key, year, week) %>%
  mutate(
    cases=round(cases/scale_korea)
  )
