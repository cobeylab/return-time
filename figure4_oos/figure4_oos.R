library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_seir.R")
load("../figure4/figure4_summ.rda")

load("../analysis_takens_acf_fnn/analysis_korea_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_piv_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_noro_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_canada_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_nrevss_takens_acf_fnn.rda")

load("../figure1/figure1.rda")

analysis_all <- bind_rows(
  analysis_korea_takens_acf_fnn %>%
    mutate(country="Korea",
           time=year+week/52),
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
                          "Rhinovirus/Enterovirus",
                          "Parainfluenza virus",
                          "RSV")),
      time=year+week/52
    ) %>%
    select(-type)
) %>%
  ungroup %>%
  mutate(
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus",
                        "Parainfluenza virus", "Rhinovirus/Enterovirus",
                        "RSV", "Human coronavirus", "Bocavirus", "Norovirus"))
  ) %>%
  filter(
    key != "Bocavirus", key != "Norovirus"
  )

analysis_all_lm <- lapply(split(analysis_all, analysis_all[,c("key", "country")]), function(x) {
  if (nrow(x) > 0) {
    x_filter <- x %>%
      filter(!is.na(dist_takens),
             dist_takens > 0) %>%
      arrange(year, week) %>%
      filter(year+week/52 < 2023+26/52)
    
    x_pre <- filter(x, year <2020)
    
    pre_mean <- mean(log(x_pre$dist_takens), na.rm=TRUE)
    
    dist <- x_filter$dist_takens
    time <- x_filter$time
    
    loesspred <- exp(predict(loess(log(dist)~time, span=0.2)))
    
    maxdist <- max(loesspred)
    maxdist_time <- time[which.max(loesspred)]
    meandist <- mean(loesspred[time<2020])
    
    target_first <- meandist * (maxdist/meandist)^(17/19)
    target_second <- meandist * (maxdist/meandist)^(2/19)
    
    time_first <- time[which(time > maxdist_time & loesspred < target_first)[1]]
    time_last <- time[tail(which(time > maxdist_time & loesspred > target_second), 1)]
    
    if (x$key[1]=="Norovirus" & x$country[1]=="Hong Kong") {
      time_first <- 2022.5
      time_last <- 2024
    } else if (x$key[1]=="Norovirus" & x$country[1]=="Korea") {
      time_first <- 2022.5
      time_last <- 2024
    } else if (x$key[1]=="Rhinovirus/Enterovirus" & x$country[1]=="US") {
      time_first <- 2021+21/52
      time_last <- 2022
    }
    
    # plot(time, dist, log="y")
    # abline(v=time_first)
    # abline(v=time_last)
    # lines(filter(x_filter, time >= time_first, time < time_last)$time, exp(predict(lfit)), col=2)
    
    lfit <- try(lm(log(dist_takens)~time, data=x_filter %>% filter(time >= time_first, time < time_last)))
    
    if (!inherits(lfit, "try-error")) {
      
      est <- -coef(lfit)[[2]]
      lwr <- -confint(lfit)[2,2]
      upr <- -confint(lfit)[2,1]
      
      time_pred <- seq(time_first, 2100, by=0.01)
      
      pp <- predict(lfit, newdata=data.frame(time=time_pred),
                    levels=0.95,
                    interval="prediction")
      
      when <- time_pred[which(pp[,1] < pre_mean)[1]]
      
      when_lwr <- time_pred[which(pp[,2] < pre_mean)[1]]
      
      when_upr <- time_pred[which(pp[,3] < pre_mean)[1]]
      
      data.frame(
        pred=exp(pp[,1]),
        pred_lwr=exp(pp[,2]),
        pred_upr=exp(pp[,3]),
        time=time_pred,
        key=x$key[1],
        country=x$country[1],
        tstart=time_first,
        tend=time_last,
        pre_mean=pre_mean,
        reldist=mean(tail(dist, 52))/mean(dist[time < 2020]),
        when=when,
        when_lwr=when_lwr,
        when_upr=when_upr,
        resilience=-coef(lfit)[[2]],
        resilience_lwr=-confint(lfit)[2,2],
        resilience_upr=-confint(lfit)[2,1]
      )
    } else {
      NULL
    }
    
  }
}) %>%
  bind_rows

analysis_all_summ <- analysis_all_lm %>%
  group_by(key, country) %>%
  filter(1:n()==1) %>%
  ungroup %>%
  mutate(
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus",
                        "Parainfluenza virus", "Rhinovirus/Enterovirus",
                        "RSV", "Human coronavirus", "Bocavirus", "Norovirus"))
  )

analysis_all_summ_filter_oos <- analysis_all_summ %>%
  filter(
    reldist < exp(2)
  )

analysis_all_lm_filter <- analysis_all_lm %>%
  merge(select(analysis_all_summ_filter_oos, key, country))

g1 <- ggplot(analysis_all) +
  geom_vline(xintercept=2013:2027, lty=3, col="gray") +
  geom_hline(data=analysis_all_summ, aes(yintercept=exp(pre_mean)), lty=2) +
  geom_line(aes(time, dist_takens)) +
  geom_ribbon(data=analysis_all_lm_filter, aes(time, ymin=pred_lwr, ymax=pred_upr), fill="red", alpha=0.2) +
  geom_line(data=analysis_all_lm_filter, aes(time, pred), col="red") +
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

ggsave("figure4_dist_oos.pdf", g1, width=12, height=6)

measles_resilience <- -eigen_seir()

analysis_all_summ_filter_oos %>%
  summarize(
    min=min(resilience),
    max=max(resilience)
  )

analysis_all_summ_filter_oos %>%
  filter(key != "Norovirus") %>%
  summarize(
    mean=mean(resilience),
    lwr=t.test(resilience)[[4]][1],
    upr=t.test(resilience)[[4]][2]
  ) %>%
  mutate(
    rel=mean/measles_resilience
  )

analysis_all_summ_filter_oos %>%
  filter(key == "Norovirus") %>%
  select(country, resilience, resilience_lwr, resilience_upr)

lfit <- lm(resilience~country+key, data=analysis_all_summ_filter_oos)

summary(lfit)

afit <- aov(resilience~country+key, data=analysis_all_summ_filter_oos)

summary(afit)

g2 <- ggplot(analysis_all_summ_filter_oos) +
  geom_hline(yintercept=measles_resilience, lty=1, col="gray", lwd=2) +
  annotate("text", x=-Inf, y=0.18, label="Prevaccination measles",
           hjust=-0.05, family="Times") +
  geom_errorbar(aes(key, ymin=resilience_lwr, ymax=resilience_upr, col=country), width=0, 
                position = position_dodge(width=0.5))+
  geom_point(aes(key, resilience, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  scale_y_log10("Resilience (1/year)", breaks=c(0.125, 0.25, 0.5, 1, 2)) +
  scale_color_viridis_d("Country", end=0.8) +
  scale_shape_manual("Country", values=c(15:18)) +
  coord_cartesian(ylim=c(0.09, 2.5)) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

analysis_all_summ_filter_oos_rename <- analysis_all_summ_filter_oos %>%
  mutate(
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus", "Parainfluenza virus",
                        "Rhinovirus/Enterovirus", "RSV", "Human coronavirus", "Norovirus"),
               labels=c("Adenovirus", "Human\nmetapneumovirus", "Parainfluenza\nvirus",
                        "Rhinovirus/\nEnterovirus", "RSV", "Human\ncoronavirus", "Norovirus"))
  )

g3 <- ggplot(analysis_all_summ_filter_oos_rename) +
  geom_hline(yintercept=2022:2028, lty=3, alpha=0.4) +
  geom_errorbar(aes(key, ymin=when_lwr, ymax=when_upr, col=country), width=0, 
                position = position_dodge(width=0.5)) +
  geom_point(aes(key, when, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  # geom_hline(yintercept = 2020:2039, lty=3, col="gray") +
  geom_hline(yintercept = 2025, lty=2) +
  scale_y_continuous("Expected\nreturn time",
                     breaks=seq(2020, 2030, by=2)) +
  scale_color_viridis_d("Country", end=0.8) +
  scale_shape_manual("Country", values=c(15:18)) +
  coord_cartesian(ylim=c(2020, 2030)) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

all_summ_comb <- analysis_all_summ_filter %>%
  select(key, country, when, when_lwr, when_upr, resilience, resilience_lwr, resilience_upr) %>%
  rename(
    when_true=when,
    when_lwr_true=when_lwr,
    when_upr_true=when_upr,
    resilience_true=resilience,
    resilience_lwr_true=resilience_lwr,
    resilience_upr_true=resilience_upr
  ) %>%
  merge(
    analysis_all_summ_filter_oos %>%
      select(key, country, when, when_lwr, when_upr, resilience, resilience_lwr, resilience_upr) 
  )

cor.test(all_summ_comb$when, all_summ_comb$when_true)

g4 <- ggplot(all_summ_comb) +
  geom_errorbar(aes(when, ymin=when_lwr_true, ymax=when_upr_true), width=0, alpha=0.4, lwd=1) +
  geom_errorbarh(aes(xmin=when_lwr, xmax=when_upr, y=when_true), height=0, alpha=0.4, lwd=1) +
  geom_point(aes(when, when_true), size=3, shape=21, fill="white", stroke=1) +
  geom_abline(intercept=0, slope=1, lty=2) +
  scale_x_continuous("Expected return time\nbased on partial fits",
                     breaks=c(2022, 2024, 2026, 2028, 2030, 2032, 2034)) +
  scale_y_continuous("Expected return time\nbased on full fits") +
  theme(
    panel.grid = element_blank()
  )

cor.test(all_summ_comb$resilience, all_summ_comb$resilience_true)

g5 <- ggplot(all_summ_comb) +
  geom_errorbar(aes(resilience, ymin=resilience_lwr_true, ymax=resilience_upr_true), width=0, alpha=0.4, lwd=1) +
  geom_errorbarh(aes(xmin=resilience_lwr, xmax=resilience_upr, y=resilience_true), height=0, alpha=0.4, lwd=1) +
  geom_point(aes(resilience, resilience_true), size=3, shape=21, fill="white", stroke=1) +
  geom_abline(intercept=0, slope=1, lty=2) +
  scale_x_continuous("Resilience estimates\nbased on partial fits (1/year)") +
  scale_y_continuous("Resilience estimates\nbased on partial fits (1/year)") +
  theme(
    panel.grid = element_blank()
  )

gcomb <- ggarrange(g2, g3,
                   labels=c("A", "B"),
                   heights=c(1, 1))

gcomb2 <- ggarrange(g4, g5, nrow=1, labels=c("C", "D"))

gfinal <- arrangeGrob(gcomb, gcomb2, nrow=2, heights=c(2, 1))

ggsave("figure4_oos.pdf", gfinal, width=8, height=8)
save("analysis_all_summ_filter_oos", file="figure4_summ_oos.rda")
