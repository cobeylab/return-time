library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
source("../script/script_data_korea.R")
source("../script/script_data_korea_int.R")
source("../script/script_data_canada_resp.R")
source("../script/script_data_hongkong_resp.R")
source("../script/script_data_hongkong_piv.R")
source("../script/script_data_hongkong_fecal.R")

data_hongkong_piv_all_scaled <- data_hongkong_piv_scaled %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases),
    key="PIV"
  )

data_hongkong_all_scaled <- bind_rows(
  data_hongkong_resp_scaled,
  data_hongkong_piv_all_scaled,
  data_hongkong_fecal_scaled %>% filter(key %in% c("noro"))
)

hongkong_pred <- lapply(split(data_hongkong_all_scaled, data_hongkong_all_scaled$key), function(x) {
  if (x$key[1]=="hmpv") {
    out <- x %>%
      filter(!is.na(cases)) %>%
      group_by(key, week) %>%
      mutate(
        mean=mean(cases[year < 2020]),
        time=year+week/52
      )
  } else if (x$key[1]=="noro") {
    out <- x %>%
      group_by(month) %>%
      mutate(
        mean=mean(cases[year < 2020]),
        lwr=t.test(cases[year < 2020])[[4]][1],
        upr=t.test(cases[year < 2020])[[4]][2],
        time=year+match(month, month.name)/12
      )
  } else {
    out <- x %>%
      filter(!is.na(cases)) %>%
      group_by(key, week) %>%
      mutate(
        mean=mean(cases[year < 2020]),
        lwr=t.test(cases[year < 2020])[[4]][1],
        upr=t.test(cases[year < 2020])[[4]][2],
        time=year+week/52
      )
    
    out 
  }
}) %>%
  bind_rows %>%
  filter(!key %in% c("rv", "ev"))  %>%
  bind_rows(
    data.frame(key=c("cov", "boca", "noro"))
  ) %>%
  mutate(
    key=factor(key,
               levels=c("adeno", "hmpv", "PIV", "rvev", "rsv", "cov", "boca", "noro"),
               labels=c("Adenovirus", "Human metapneumovirus",
                        "Parainfluenza virus", "Rhinovirus/Enterovirus",
                        "RSV", "Human coronavirus", "Bocavirus", "Norovirus"))
  )

data_korea_all_scaled <- bind_rows(
  data_korea_ari_scaled,
  filter(data_korea_int_scaled, key %in% c("Norovirus"))
)

korea_pred <- lapply(split(data_korea_all_scaled, data_korea_all_scaled$key), function(x) {
  y <- x %>% 
    filter(year < 2020) 
  
  x %>%
    group_by(key, week) %>%
    mutate(
      mean=mean(cases[year < 2020]),
      lwr=t.test(cases[year < 2020])[[4]][1],
      upr=t.test(cases[year < 2020])[[4]][2]
    )
}) %>%
  bind_rows %>%
  mutate(
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus",
                         "Parainfluenza virus", "Rhinovirus",
                         "RSV", "Human coronavirus", "Bocavirus", "Norovirus"))
  )

canada_pred <- lapply(split(data_canada_resp_scaled, data_canada_resp_scaled$key), function(x) {
  y <- x %>% 
    filter(year < 2020) 
  
  x %>%
    group_by(key, week) %>%
    mutate(
      mean=mean(cases[year < 2020]),
      lwr=t.test(cases[year < 2020])[[4]][1],
      upr=t.test(cases[year < 2020])[[4]][2]
    )
}) %>%
  bind_rows() %>%
  bind_rows(
    data.frame(key=c("boca", "noro"))
  ) %>%
  mutate(
    key=factor(key,
                levels=c("AdV", "HMPV", "PIV", "RV/EV", "RSV",
                         "Flu A", "Flu B", "CoV", "boca", "noro"),
                labels=c("Adenovirus", "Human metapneumovirus",
                         "Parainfluenza virus", "Rhinovirus/Enterovirus",
                         "RSV", "Flu A", "Flu B", "Human coronavirus",
                         "Bocavirus", "Norovirus"))
  ) %>% 
  filter(key %in% c("Adenovirus", "Human coronavirus", "Human metapneumovirus",
                          "Parainfluenza virus", "Rhinovirus/Enterovirus",
                          "RSV", "Bocavirus", "Norovirus"))

g1 <- ggplot(canada_pred) +
  geom_vline(xintercept=2014:2024, lty=3, col="gray70", lwd=0.5) +
  geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), fill="gray", col="gray90", alpha=0.3, lwd=0.5) +
  geom_line(aes(year+week/52, mean), col="gray", lwd=0.7) +
  geom_line(data=filter(canada_pred, year<2020), aes(year+week/52, cases), col="#EF6351", lwd=0.5) +
  geom_point(data=filter(canada_pred, year<2020), aes(year+week/52, cases), col="#EF6351", size=0.7) +
  geom_line(data=filter(canada_pred, year>=2020), aes(year+week/52, cases), col="#224B95", lwd=0.5) +
  geom_point(data=filter(canada_pred, year>=2020), aes(year+week/52, cases), col="#224B95", size=0.7) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2013:2024) +
  scale_y_continuous("Scaled cases", expand=c(0, 0)) +
  facet_wrap(~key, scale="free", ncol=1) +
  ggtitle("Canada") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

g2 <- ggplot(hongkong_pred) +
  geom_vline(xintercept=2014:2024, lty=3, col="gray70", lwd=0.5) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), fill="gray", col="gray90", alpha=0.3, lwd=0.5) +
  geom_line(aes(time, mean), col="gray", lwd=0.7) +
  geom_line(data=filter(hongkong_pred, year<2020), aes(time, cases), col="#EF6351", lwd=0.5) +
  geom_point(data=filter(hongkong_pred, year<2020), aes(time, cases), col="#EF6351", size=0.7) +
  geom_line(data=filter(hongkong_pred, year>=2020), aes(time, cases), col="#224B95", lwd=0.5) +
  geom_point(data=filter(hongkong_pred, year>=2020), aes(time, cases), col="#224B95", size=0.7) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2013:2024) +
  scale_y_continuous("Scaled cases", expand=c(0, 0)) +
  facet_wrap(~key, scale="free", nrow=8) +
  ggtitle("Hong Kong") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

g3 <- ggplot(korea_pred) +
  geom_vline(xintercept=2016:2024, lty=3, col="gray70", lwd=0.5) +
  geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), fill="gray", col="gray90", alpha=0.3, lwd=0.5) +
  geom_line(aes(year+week/52, mean), col="gray", lwd=0.7) +
  geom_line(data=filter(korea_pred, year<2020), aes(year+week/52, cases), col="#EF6351", lwd=0.5) +
  geom_point(data=filter(korea_pred, year<2020), aes(year+week/52, cases), col="#EF6351", size=0.7) +
  geom_line(data=filter(korea_pred, year>=2020), aes(year+week/52, cases), col="#224B95", lwd=0.5) +
  geom_point(data=filter(korea_pred, year>=2020), aes(year+week/52, cases), col="#224B95", size=0.7) +
  scale_x_continuous("Year", expand=c(0, 0),
                     breaks=2015:2024) +
  scale_y_continuous("Scaled cases", expand=c(0, 0)) +
  facet_wrap(~key, scale="free", ncol=1) +
  ggtitle("Korea") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(linewidth=1),
    axis.title.x = element_blank()
  )

gcomb <- ggarrange(g1, g2, g3, ncol=3, labels=c("A", "B", "C"))

ggsave("figure1.pdf", gcomb, width=12, height=8)
