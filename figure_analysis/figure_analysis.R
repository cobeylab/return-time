library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_return_time/analysis_korea_growth.rda")
load("../analysis_return_time/analysis_hongkong_resp_growth.rda")
load("../analysis_return_time/analysis_hongkong_piv_growth.rda")

summ_korea_resilience <- analysis_korea_growth %>%
  group_by(key, id) %>%
  summarize(
    return_rC=unique(return_rC),
    return_takens=unique(return_takens)
  ) %>%
  gather(
    method, value, -id, -key
  ) %>%
  mutate(
    method=factor(method,
                  levels=c("return_rC", "return_takens"),
                  labels=c("growth rate", "Takens' theorem")),
    country="Korea"
  )

summ_hongkong_resp_resilience <- analysis_hongkong_resp_growth %>%
  group_by(key, id) %>%
  summarize(
    return_rC=unique(return_rC),
    return_takens=unique(return_takens)
  ) %>%
  gather(
    method, value, -id, -key
  ) %>%
  mutate(
    method=factor(method,
                  levels=c("return_rC", "return_takens"),
                  labels=c("growth rate", "Takens' theorem")),
    country="Hong Kong",
    key=factor(key,
               levels=c("adeno", "hmpv", "rsv", "rvev"),
               labels=c("Adenovirus", "Human metapneumovirus",
                        "RSV", "Rhinovirus/Enterovirus"))
  )

summ_hongkong_piv_resilience <- analysis_hongkong_piv_growth %>%
  group_by(key, id) %>%
  summarize(
    return_rC=unique(return_rC),
    return_takens=unique(return_takens)
  ) %>%
  gather(
    method, value, -id, -key
  ) %>%
  mutate(
    method=factor(method,
                  levels=c("return_rC", "return_takens"),
                  labels=c("growth rate", "Takens' theorem")),
    country="Hong Kong",
    key="Parainfluenza virus"
  )

summ_resilience <-
  bind_rows(
    summ_korea_resilience,
    summ_hongkong_resp_resilience,
    summ_hongkong_piv_resilience
  )

g1 <- ggplot(summ_resilience) +
  geom_violin(aes(key, value, fill=country)) +
  scale_y_log10("Return time (years)") +
  facet_wrap(~method) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_blank(),
    legend.position = "top"
  )

predall <- bind_rows(
  analysis_hongkong_growth_pred %>%
    mutate(
      country="Hong Kong",
      key=factor(key,
                 levels=c("adeno", "hmpv", "rsv", "rvev"),
                 labels=c("Adenovirus", "Human metapneumovirus",
                          "RSV", "Rhinovirus/Enterovirus"))
    ),
  analysis_hongkong_piv_growth_pred %>% mutate(country="Hong Kong", labels="Parainfluenza virus"),
  analysis_korea_growth_pred %>% mutate(country="Korea")
)

g2 <- ggplot(predall %>% filter(id <= 50, method=="Takens' theorem")) +
  geom_line(aes(time, pred, group=id)) +
  scale_y_log10("Distance from attractor") +
  scale_color_viridis_d() +
  coord_cartesian(ylim=c(0.1, 50)) +
  facet_grid(country~key) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )

gcomb <- ggarrange(g1, g2, nrow=2)

ggsave("figure_analysis.pdf", gcomb, width=12, height=8)
