data_hongkong_piv <- read.csv("../data_processed/data_hongkong_piv.csv") %>%
  mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  group_by(year, week) %>%
  summarize(
    piv1=round(mean(piv1)),
    piv2=round(mean(piv2)),
    piv3=round(mean(piv3)),
    piv4=round(mean(piv4)),
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
  ) %>%
  mutate(
    key=factor(key, labels=c("PIV 1", "PIV 2", "PIV 3", "PIV 4"))
  )

data_hongkong_scale_piv <- data_hongkong_piv %>%
  filter(key=="PIV 1")

roll_total_piv <- zoo::rollapply(filter(data_hongkong_scale_piv, year < 2020)$test,
                             width=2*52, mean, partial=TRUE)

scale_hongkong_piv <- c(
  roll_total_piv/tail(roll_total_piv, 1),
  rep(1, nrow(data_hongkong_scale_piv)-length(roll_total_piv))
)

data_hongkong_piv_scaled <- data_hongkong_piv %>%
  group_by(key) %>%
  arrange(key, year, week) %>%
  mutate(
    cases=round(cases/scale_hongkong_piv)
  )
