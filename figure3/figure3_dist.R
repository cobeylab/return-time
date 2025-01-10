library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
load("../analysis_takens_acf_fnn/analysis_korea_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_piv_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_noro_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_canada_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_nrevss_takens_acf_fnn.rda")

analysis_all <- bind_rows(
  analysis_korea_takens_acf_fnn %>%
    mutate(country="Korea",
           time=year+week/52),
  analysis_hongkong_takens_acf_fnn  %>%
    mutate(
      country="Hong Kong",
      key=factor(key,
                 levels=c("adeno", "hmpv", "rsv", "rvev"),
                 labels=c("Adenovirus", "Human metapneumovirus",
                          "RSV", "Rhinovirus/Enterovirus")),
      time=year+week/52
    ),
  analysis_hongkong_piv_takens_acf_fnn %>%
    mutate(
      country="Hong Kong",
      key="Parainfluenza virus",
      time=year+week/52
    ),
  analysis_hongkong_noro_takens_acf_fnn %>%
    mutate(
      country="Hong Kong",
      time=year+match(month,month.name)/12
    ),
  analysis_canada_takens_acf_fnn %>%
    filter(
      !key %in% c("Flu A", "Flu B")
    ) %>%
    mutate(
      country="Canada",
      key=factor(key,
                 levels=c("AdV", "CoV", "Flu A", "Flu B", "HMPV", "PIV", "RSV", "RV/EV"),
                 labels=c("Adenovirus", "Human coronavirus", "Influenza A", "Influenza B",
                          "Human metapneumovirus",
                          "Parainfluenza virus",
                          "RSV", "Rhinovirus/Enterovirus")),
      time=year+week/52
    ),
  analysis_nrevss_takens_acf_fnn %>%
    mutate(
      country="US",
      key=factor(type,
                 levels=c("Adenovirus", "CoV", "Human metapneumovirus", "Rhinovirus", "PIV", "RSV"),
                 labels=c("Adenovirus", "Human coronavirus", "Human metapneumovirus",
                          "Rhinovirus",
                          "Parainfluenza virus",
                          "RSV")),
      time=year+week/52
    ) %>%
    select(-type)
) %>%
  mutate(
    key=ifelse(key=="Rhinovirus", "Rhinovirus/Enterovirus", key),
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus",
                        "Parainfluenza virus", "Rhinovirus/Enterovirus",
                        "RSV", "Human coronavirus", "Bocavirus", "Norovirus"))
  )

analysis_all_summ <- analysis_all %>%
  group_by(key, country) %>%
  filter(!is.na(dist_takens)) %>%
  arrange(key, country, year, week) %>%
  summarize(
    mean_pre=mean(dist_takens[year<2020], na.rm=TRUE),
    mean_recent=ifelse(key[1]=="Norovirus" & country[1]=="Hong Kong", mean(tail(dist_takens, 12)), mean(tail(dist_takens, 52)))
  )

g1 <- ggplot(analysis_all) +
  geom_vline(xintercept=2013:2027, lty=3, col="gray") +
  geom_hline(data=analysis_all_summ, aes(yintercept=mean_pre), lty=2) +
  geom_line(aes(time, dist_takens)) +
  # geom_ribbon(data=analysis_all_lm, aes(time, ymin=pred_lwr, ymax=pred_upr), fill="red", alpha=0.2) +
  # geom_line(data=analysis_all_lm, aes(time, pred), col="red") +
  scale_x_continuous("Year",
                     breaks=seq(2014, 2030, by=2),
                     expand=c(0, 0)) +
  scale_y_log10("Nearest neighbor distance from attractor", expand=c(0, 0)) +
  coord_cartesian(xlim=c(2013, 2025), ylim=c(5e-2, 20)) +
  facet_grid(country~key) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  )

g2 <- ggplot(analysis_all_summ) +
  geom_smooth(aes(key, mean_recent/mean_pre, group=country, col=country, fill=country), 
              method="lm", formula=y~1, lty=2, lwd=0.7,
              fullrange=TRUE) +
  geom_point(aes(key, mean_recent/mean_pre, col=country), size=3) +
  scale_y_log10("Relative distance from the attractor") +
  facet_wrap(~country, ncol=1) +
  coord_flip() +
  theme(
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("figure3_dist.pdf", g1, width=12, height=6)
ggsave("figure3_dist_rel.pdf", g2, width=4, height=6)
