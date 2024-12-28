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
    key=ifelse(key=="Rhinovirus", "Rhinovirus/Enterovirus", key)
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
    
    loess_fit <- loess(log(dist_takens)~time, data=x_filter, span=0.2)
    pred <- c(exp(predict(loess_fit)))
    
    time <- x_filter$time
    ## dt
    dpred <- c(NA, diff(pred))
    
    deriv_max <- c(NA, (head(dpred, -1) > 0 & tail(dpred, -1) < 0)) & time > 2020 &
      pred > 0.5 * max(pred)
    
    tmax <- tail(which(deriv_max), 1)
    max_pred <- pred[tmax]
    
    deriv_zero <- (c(NA, (head(dpred, -1) < 0 & tail(dpred, -1) > 0))) & 
      time > time[tmax]
    
    if (any(deriv_zero)) {
      final_pred <- pred[which(deriv_zero)[1]]
      
      ratio <- max_pred/final_pred
      
      tstart <- time[tail(which(pred > ratio^(8/10) * final_pred & time < time[which(deriv_zero)[1]]), 1)]
      
      tend <- time[tail(which(pred > ratio^(0.5/10) * final_pred  & time <= time[which(deriv_zero)[1]]), 1)]
    } else {
      final_pred <- tail(pred, 1)
      
      ratio <- max_pred/final_pred
      
      tstart <- time[tail(which(pred > ratio^(8/10) * final_pred), 1)]
      
      tend <- time[tail(which(pred > ratio^(0.5/10) * final_pred), 1)]
    }
    
    x_filter2 <- x_filter %>%
      filter(
        time >= tstart, time <= tend
      )
    
    lfit <- lm(log(dist_takens)~time, data=x_filter2)
    
    time_pred <- seq(tstart, 2100, by=0.01)
    
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
      tstart=tstart,
      tend=tend,
      pre_mean=pre_mean,
      when=when,
      when_lwr=when_lwr,
      when_upr=when_upr,
      resilience=-coef(lfit)[[2]],
      resilience_lwr=-confint(lfit)[2,2],
      resilience_upr=-confint(lfit)[2,1]
    )
  }
}) %>%
  bind_rows %>%
  mutate(
    key=ifelse(key=="Rhinovirus", "Rhinovirus/Enterovirus", key)
  )

analysis_all_summ <- analysis_all_lm %>%
  group_by(key, country) %>%
  filter(1:n()==1)

g1 <- ggplot(analysis_all) +
  geom_vline(xintercept=2013:2027, lty=3, col="gray") +
  geom_hline(data=analysis_all_summ, aes(yintercept=exp(pre_mean)), lty=2) +
  geom_line(aes(time, dist_takens)) +
  geom_ribbon(data=analysis_all_lm, aes(time, ymin=pred_lwr, ymax=pred_upr), fill="red", alpha=0.2) +
  geom_line(data=analysis_all_lm, aes(time, pred), col="red") +
  scale_x_continuous("Year",
                     breaks=seq(2014, 2030, by=2),
                     expand=c(0, 0)) +
  scale_y_log10("Nearest neighbor distance from attractor", expand=c(0, 0)) +
  coord_cartesian(xlim=c(2013, 2027), ylim=c(5e-2, 20)) +
  facet_grid(country~key) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  )

g2 <- ggplot(analysis_all_summ) +
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

g3 <- ggplot(analysis_all_summ) +
  geom_errorbar(aes(key, ymin=when_lwr, ymax=when_upr, col=country), width=0, 
                position = position_dodge(width=0.5)) +
  geom_point(aes(key, when, col=country, shape=country), 
             position = position_dodge(width=0.5), size=3) +
  geom_hline(yintercept = 2020:2039, lty=3, col="gray") +
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

gcomb <- ggarrange(g2, g3,
                   labels=c("A", "B"))

ggsave("figure3.pdf", gcomb, width=6, height=6)
