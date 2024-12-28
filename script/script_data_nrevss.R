data_nrevss_resp <- read.csv("../data_processed/data_nrevss_resp.csv") %>%
  dplyr::mutate(
    date=as.Date(date)
  ) %>%
  dplyr::mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  dplyr::group_by(region, type, year, week) %>%
  dplyr::summarize(
    tests=round(mean(tests)),
    positive=round(mean(positive))
  )

data_nrevss_resp$tests[15791] <- 3644
data_nrevss_resp$tests[22850] <- 5846

data_nrevss_resp_comb <- data_nrevss_resp %>%
  group_by(year, week, type) %>%
  summarize(
    tests=sum(tests),
    positive=sum(positive)
  ) %>%
  arrange(type, year, week)

data_ili_raw <- vroom::vroom("../data_nrevss/ILINet.csv", delim =",", skip=1)

data_ili_new <- data_ili_raw %>%
  dplyr::filter(`REGION TYPE`=="National", YEAR+WEEK/52 >= 2013) %>%
  dplyr::select(YEAR, WEEK, `% WEIGHTED ILI`) %>%
  rename(year=YEAR, week=WEEK) %>%
  mutate(
    week=ifelse(week==53, 52, week)
  ) %>%
  group_by(year, week) %>%
  summarize(
    ili=mean(as.numeric(`% WEIGHTED ILI`))
  )

data_nrevss_resp_comb_proxy <- data_nrevss_resp_comb %>%
  merge(data_ili_new) %>%
  mutate(
    proxy=positive/tests*ili
  )
