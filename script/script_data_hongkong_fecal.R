data_hongkong_fecal <- read.csv("../data_processed/data_hongkong_fecal.csv") %>%
  dplyr::select(-X) %>%
  tidyr::gather(key, value, -year, -month, -test) %>%
  dplyr::rename(
    cases=value
  ) %>%
  mutate(
    month=factor(month,
                 levels=month.name)
  )

data_hongkong_scale <- data_hongkong_fecal %>%
  filter(key=="noro")

roll_total <- zoo::rollapply(filter(data_hongkong_scale, year < 2020)$test,
                             width=2*52, mean, partial=TRUE)

scale_hongkong <- c(
  roll_total/tail(roll_total, 1),
  rep(1, nrow(data_hongkong_scale)-length(roll_total))
)

data_hongkong_fecal_scaled <- data_hongkong_fecal %>%
  group_by(key) %>%
  arrange(key, year, month) %>%
  mutate(
    cases=round(cases/scale_hongkong)
  )
