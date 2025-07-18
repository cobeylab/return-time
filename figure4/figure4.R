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

load("../figure1/figure1.rda")

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
  ) %>%
  filter(
    key != "Bocavirus"
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
  summarize(
    min=min(resilience),
    max=max(resilience)
  )

analysis_all_summ_filter %>%
  filter(key != "Norovirus") %>%
  summarize(
    mean=mean(resilience),
    lwr=t.test(resilience)[[4]][1],
    upr=t.test(resilience)[[4]][2]
  ) %>%
  mutate(
    rel=mean/measles_resilience
  )

analysis_all_summ_filter %>%
  filter(key == "Norovirus") %>%
  select(country, resilience, resilience_lwr, resilience_upr)

lfit <- lm(resilience~country+key, data=analysis_all_summ_filter)

summary(lfit)

afit <- aov(resilience~country+key, data=analysis_all_summ_filter)

summary(afit)

g2 <- ggplot(analysis_all_summ_filter) +
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
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

analysis_all_summ_filter_rename <- analysis_all_summ_filter %>%
  mutate(
    key=factor(key,
               levels=c("Adenovirus", "Human metapneumovirus", "Parainfluenza virus",
                        "Rhinovirus/Enterovirus", "RSV", "Human coronavirus", "Norovirus"),
               labels=c("Adenovirus", "Human\nmetapneumovirus", "Parainfluenza\nvirus",
                        "Rhinovirus/\nEnterovirus", "RSV", "Human\ncoronavirus", "Norovirus"))
  )

g3 <- ggplot(analysis_all_summ_filter_rename) +
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
  # coord_cartesian(ylim=c(2020, 2030)) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

predall <- bind_rows(
  canada_pred %>% 
    select(key, year, week, cases, mean, lwr, upr) %>% 
    mutate(country="Canada",
           time=year+week/52),
  hongkong_pred %>% 
    select(key, time, month, cases, mean, lwr, upr) %>% 
    mutate(country="Hong Kong"),
  korea_pred %>% 
    select(key, year, week, cases, mean, lwr, upr) %>% 
    mutate(country="Korea",
           time=year+week/52),
  nrevss_pred %>% 
    select(type, year, week, proxy, mean, lwr, upr) %>% 
    rename(cases=proxy, key=type) %>%
    mutate(country="US",
           time=year+week/52)
)

analysis_all_summ_filter %>%
  summarize(
    before=sum(when<2025),
    after=sum(when>=2025)
  )

analysis_all_summ_filter2 <- analysis_all_summ_filter %>%
  filter(when < 2025) %>%
  select(key, country, when, when_lwr, when_upr)

predall_filter <- predall %>%
  merge(analysis_all_summ_filter2) %>%
  mutate(
    group=paste0(key, ", ", country)
  )

g4 <- ggplot(predall_filter) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), fill="gray", col="gray90", alpha=0.3, lwd=0.5) +
  geom_line(aes(time, mean), col="gray", lwd=0.7) +
  geom_line(aes(time, cases), col="#EF6351", lwd=0.5) +
  geom_point(aes(time, cases), col="#EF6351", size=0.7) +
  geom_vline(data=analysis_all_summ_filter2 %>% mutate(group=paste0(key, ", ", country)),
             aes(xintercept=when), lty=2, col="orange", lwd=1) +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Cases/incidence proxy", expand=c(0, 0)) +
  coord_cartesian(xlim=c(2021, 2025), ylim=c(0, NA)) +
  facet_wrap(~group, scale="free", ncol=3) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

predall_filter_no <- predall %>%
  merge(analysis_all_summ_filter %>%
          filter(when >= 2025) %>%
          select(key, country, when, when_lwr, when_upr)) %>%
  mutate(
    group=paste0(key, ", ", country)
  )

g5 <- ggplot(predall_filter_no) +
  geom_ribbon(aes(time, ymin=lwr, ymax=upr), fill="gray", col="gray90", alpha=0.3, lwd=0.5) +
  geom_line(aes(time, mean), col="gray", lwd=0.7) +
  geom_line(aes(time, cases), col="#EF6351", lwd=0.5) +
  geom_point(aes(time, cases), col="#EF6351", size=0.7) +
  scale_x_continuous("Year", expand=c(0, 0)) +
  scale_y_continuous("Cases/incidence proxy", expand=c(0, 0)) +
  coord_cartesian(xlim=c(2021, 2025), ylim=c(0, NA)) +
  facet_wrap(~group, scale="free", ncol=3) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

gcomb <- ggarrange(g2, g3, g4,
                   labels=c("A", "B", "C"),
                   heights=c(1, 1, 3))

ggsave("figure4.pdf", gcomb, width=8, height=8)
ggsave("figure4_noreturn.pdf", g5, width=8, height=6)
save("analysis_all_summ_filter", file="figure4_summ.rda")
