library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_seir.R")

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

analysis_all_lm <- lapply(split(analysis_all, analysis_all[,c("key", "country")]), function(x) {
  if (nrow(x) > 0) {
    # print(x$key[1])
    x_filter <- x %>%
      filter(!is.na(dist_takens),
             dist_takens > 0) %>%
      arrange(year, week)
    
    x_pre <- filter(x, year <2020)
    
    pre_mean <- mean(log(x_pre$dist_takens), na.rm=TRUE)
    
    dist <- x_filter$dist_takens
    time <- x_filter$time
    
    loesspred <- exp(predict(loess(log(dist)~time, span=0.2)))
    
    maxdist <- max(loesspred)
    maxdist_time <- time[which.max(loesspred)]
    meandist <- mean(loesspred[time<2020])
    
    target_first <- meandist * (maxdist/meandist)^(9/10)
    target_second <- meandist * (maxdist/meandist)^(1/10)
    
    time_first <- time[which(time > maxdist_time & loesspred < target_first)[1]]
    time_last <- time[tail(which(time > maxdist_time & loesspred > target_second), 1)]
    
    # plot(time, dist, log="y")
    # abline(v=time_first)
    # abline(v=time_last)
    # lines(filter(x_filter, time >= time_first, time < time_last)$time, exp(predict(lfit)), col=2)
    
    lfit <- lm(log(dist_takens)~time, data=x_filter %>% filter(time >= time_first, time < time_last))
    
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

analysis_all_summ_filter <- analysis_all_summ %>%
  filter(
    reldist < exp(2)
  )

analysis_all_lm_filter <- analysis_all_lm %>%
  merge(select(analysis_all_summ_filter, key, country))

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

ggsave("figure4_dist.pdf", g1, width=12, height=6)

measles_resilience <- -eigen_seir()

analysis_all_summ_filter %>%
  filter(key != "Norovirus", key != "Bocavirus",
         !(key == "Rhinovirus/Enterovirus" & country=="US")) %>%
  summarize(
    min=min(resilience),
    max=max(resilience),
    mean=mean(resilience),
    lwr=t.test(resilience)[[4]][1],
    upr=t.test(resilience)[[4]][2]
  )

analysis_all_summ_filter %>%
  filter(key == "Norovirus")

lfit <- lm(resilience~country+key, data=analysis_all_summ_filter)

summary(lfit)

afit <- aov(resilience~country+key, data=analysis_all_summ_filter)

summary(afit)

g2 <- ggplot(analysis_all_summ_filter) +
  geom_hline(yintercept=measles_resilience, lty=1, col="gray", lwd=3) +
  annotate("text", x=-Inf, y=0.3, label="Prevaccination measles",
           hjust=-0.05, family="Times") +
  geom_errorbar(aes(key, ymin=resilience_lwr, ymax=resilience_upr, col=country), width=0, 
                position = position_dodge(width=0.5))+
  geom_point(aes(key, resilience, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  scale_y_continuous("Resilience (1/year)") +
  scale_color_viridis_d("Country") +
  scale_shape_discrete("Country") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

g3 <- ggplot(analysis_all_summ_filter) +
  geom_hline(yintercept=2020:2030, lty=3, alpha=0.4) +
  geom_errorbar(aes(key, ymin=when_lwr, ymax=when_upr, col=country), width=0, 
                position = position_dodge(width=0.5)) +
  geom_point(aes(key, when, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  # geom_hline(yintercept = 2020:2039, lty=3, col="gray") +
  geom_hline(yintercept = 2025, lty=2) +
  scale_y_continuous("Expected return time",
                     breaks=seq(2020, 2030, by=2)) +
  scale_color_viridis_d("Country") +
  scale_shape_discrete("Country") +
  coord_cartesian(ylim=c(2020, 2030)) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

gcomb <- ggarrange(g2, g3,
                   labels=c("A", "B"))

ggsave("figure4.pdf", gcomb, width=6, height=6)
save("analysis_all_summ_filter", file="figure4_summ.rda")

g3b <- ggplot(analysis_all_summ_filter) +
  geom_errorbar(aes(key, ymin=when_lwr, ymax=when_upr, col=country), width=0, 
                position = position_dodge(width=0.5)) +
  geom_point(aes(key, when, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  # geom_hline(yintercept = 2020:2039, lty=3, col="gray") +
  geom_hline(yintercept = 2025, lty=2) +
  scale_y_continuous("Expected return time") +
  scale_color_viridis_d("Country") +
  scale_shape_discrete("Country") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

ggsave("figure4_unzoom.pdf", g3b, width=6, height=3)
