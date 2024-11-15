library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_takens_acf_2d/analysis_korea_takens_acf_2d.rda")
load("../analysis_takens_acf_2d/analysis_hongkong_takens_acf_2d.rda")
load("../analysis_takens_acf_2d/analysis_hongkong_piv_takens_acf_2d.rda")
load("../analysis_takens_acf_2d/analysis_hongkong_noro_takens_acf_2d.rda")
load("../analysis_takens_acf_2d/analysis_canada_takens_acf_2d.rda")

summ_korea_takens_acf_2d <- analysis_korea_takens_acf_2d %>%
  group_by(key) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Korea"
  )

summ_hongkong_takens_acf_2d <- analysis_hongkong_takens_acf_2d %>%
  group_by(key) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Hong Kong",
    key=factor(key,
               levels=c("adeno", "hmpv", "rsv", "rvev"),
               labels=c("Adenovirus", "Human metapneumovirus",
                        "RSV", "Rhinovirus/Enterovirus"))
  )

summ_hongkong_piv_takens_acf_2d <- analysis_hongkong_piv_takens_acf_2d %>%
  group_by(key) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Hong Kong",
    key="Parainfluenza virus"
  )

summ_hongkong_noro_takens_acf_2d <- analysis_hongkong_noro_takens_acf_2d %>%
  group_by(key) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Hong Kong"
  )

summ_canada_takens_acf_2d <- analysis_canada_takens_acf_2d %>%
  group_by(key) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Canada",
    key=factor(key,
               levels=c("AdV", "CoV", "Flu A", "Flu B", "HMPV", "PIV", "RSV", "RV/EV"),
               labels=c("Adenovirus", "Human coronavirus", "Influenza A", "Influenza B",
                        "Human metapneumovirus",
                        "Parainfluenza virus",
                        "RSV", "Rhinovirus/Enterovirus"))
  )

summ_resilience_acf_2d <-
  bind_rows(
    summ_korea_takens_acf_2d,
    summ_hongkong_takens_acf_2d,
    summ_hongkong_piv_takens_acf_2d,
    summ_hongkong_noro_takens_acf_2d,
    summ_canada_takens_acf_2d
  )

g1 <- ggplot(summ_resilience_acf_2d) +
  geom_point(aes(key, -resilience, col=country),
             position=position_dodge(width=0.5)) +
  geom_errorbar(aes(key, ymin=-resilience_lwr, ymax=-resilience_upr, col=country), width=0,
                position=position_dodge(width=0.5)) +
  scale_color_viridis_d("Country", end=0.8) +
  scale_y_continuous("Resilience (1/years)") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

ggsave("figure_takens_acf_2d.pdf", g1, width=8, height=6)

analysis_all <- bind_rows(
  analysis_korea_takens_acf_2d %>%
    mutate(country="Korea"),
  analysis_hongkong_takens_acf_2d  %>%
    mutate(
      country="Hong Kong",
      key=factor(key,
                 levels=c("adeno", "hmpv", "rsv", "rvev"),
                 labels=c("Adenovirus", "Human metapneumovirus",
                          "RSV", "Rhinovirus/Enterovirus"))
    ),
  analysis_hongkong_piv_takens_acf_2d %>%
    mutate(
      country="Hong Kong",
      key="Parainfluenza virus"
    ),
  analysis_hongkong_noro_takens_acf_2d %>%
    mutate(
      country="Hong Kong"
    ),
  analysis_canada_takens_acf_2d %>%
    mutate(
      country="Canada",
      key=factor(key,
                 levels=c("AdV", "CoV", "Flu A", "Flu B", "HMPV", "PIV", "RSV", "RV/EV"),
                 labels=c("Adenovirus", "Human coronavirus", "Influenza A", "Influenza B",
                          "Human metapneumovirus",
                          "Parainfluenza virus",
                          "RSV", "Rhinovirus/Enterovirus"))
    )
)

g2 <- ggplot(analysis_all) +
  geom_line(aes(year+week/52, dist_takens)) +
  scale_x_continuous("Year") +
  scale_y_log10("Distance from attractor") +
  facet_grid(country~key)

ggsave("figure_takens_acf_2d_dist.pdf", width=12, height=6)