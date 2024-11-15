library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_takens_acf_nd/analysis_korea_takens_acf_nd.rda")
load("../analysis_takens_acf_nd/analysis_hongkong_takens_acf_nd.rda")
load("../analysis_takens_acf_nd/analysis_hongkong_piv_takens_acf_nd.rda")
load("../analysis_takens_acf_nd/analysis_hongkong_noro_takens_acf_nd.rda")
load("../analysis_takens_acf_nd/analysis_canada_takens_acf_nd.rda")

summ_korea_takens_acf_nd <- analysis_korea_takens_acf_nd %>%
  group_by(key, d) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Korea"
  )

summ_hongkong_takens_acf_nd <- analysis_hongkong_takens_acf_nd %>%
  group_by(key, d) %>%
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

summ_hongkong_piv_takens_acf_nd <- analysis_hongkong_piv_takens_acf_nd %>%
  group_by(key, d) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Hong Kong",
    key="Parainfluenza virus"
  )

summ_hongkong_noro_takens_acf_nd <- analysis_hongkong_noro_takens_acf_nd %>%
  group_by(key, d) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Hong Kong"
  )

summ_canada_takens_acf_nd <- analysis_canada_takens_acf_nd %>%
  group_by(key, d) %>%
  summarise(
    resilience=resilience[1],
    resilience_lwr=resilience_lwr[1],
    resilience_upr=resilience_upr[1]
  ) %>%
  mutate(
    country="Canada",
    key=factor(key,
               levels=c("AdV", "CoV", "Flu A", "Flu B", "HMPV", "PIV", "RSV", "RV/EV"),
               labels=c("Adenovirus", "Influenza A", "Influenza B", "Human coronavirus",
                        "Human metapneumovirus",
                        "Parainfluenza virus",
                        "RSV", "Rhinovirus/Enterovirus"))
  )


summ_resilience_acf_nd <-
  bind_rows(
    summ_korea_takens_acf_nd,
    summ_hongkong_takens_acf_nd,
    summ_hongkong_piv_takens_acf_nd,
    summ_hongkong_noro_takens_acf_nd,
    summ_canada_takens_acf_nd
  )

g1 <- ggplot(summ_resilience_acf_nd) +
  geom_point(aes(d, -resilience, col=country),
             position=position_dodge(width=0.5)) +
  geom_errorbar(aes(d, ymin=-resilience_lwr, ymax=-resilience_upr, col=country), width=0,
                position=position_dodge(width=0.5)) +
  facet_grid(country~key) +
  scale_x_continuous("Dimensions") +
  scale_y_continuous("Resilience (1/years)") +
  scale_color_viridis_d("Country", end=0.8) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

ggsave("figure_takens_acf_nd.pdf", g1, width=12, height=6)
