data_hongkong_resp <- read.csv("../data_processed/data_hongkong_resp.csv") %>%
  mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  group_by(year, week) %>%
  summarize(
    adeno=round(mean(adeno)),
    rsv=round(mean(rsv)),
    rvev=round(mean(rvev)),
    rv=round(mean(rv)),
    ev=round(mean(ev)),
    hmpv=round(mean(hmpv)),
    test=round(mean(test)),
    age=unique(age)
  ) %>%
  gather(key, value, -year, -week, -age, -test) %>%
  group_by(year, week, key) %>%
  summarize(
    value=sum(value), test=sum(test)
  ) %>%
  rename(
    cases=value
  )

data_hongkong_scale <- data_hongkong_resp %>%
  filter(key=="adeno")

roll_total <- zoo::rollapply(filter(data_hongkong_scale, year < 2020)$test,
                             width=2*52, mean, partial=TRUE)

scale_hongkong <- c(
  roll_total/tail(roll_total, 1),
  rep(1, nrow(data_hongkong_scale)-length(roll_total))
)

data_hongkong_resp_scaled <- data_hongkong_resp %>%
  group_by(key) %>%
  arrange(key, year, week) %>%
  mutate(
    cases=round(cases/scale_hongkong)
  )
  