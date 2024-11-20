library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
load("../analysis_takens_acf_fnn/analysis_korea_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_piv_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_hongkong_noro_takens_acf_fnn.rda")
load("../analysis_takens_acf_fnn/analysis_canada_takens_acf_fnn.rda")

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
    )
)

analysis_all_lm <- lapply(split(analysis_all, analysis_all[,c("key", "country")]), function(x) {
  if (nrow(x) > 0 & x$key[1] != "Influenza B") {
    x_filter <- x %>%
      filter(!is.na(dist_takens),
             dist_takens > 0)
    
    loess_fit <- loess(log(dist_takens)~time, data=x_filter)
    pred <- exp(predict(loess_fit))
    
    time <- x_filter$time
    wmax <- which.max(pred)
    
    t_max_pred <- x_filter$time[wmax]
    t_min_pred <- x_filter$time[which.min(pred[(1:length(pred))>wmax])+wmax]
    
    tstart <- (t_min_pred-t_max_pred)*1/4+t_max_pred
    tend <- (t_min_pred-t_max_pred)*3/4+t_max_pred
    
    x_filter2 <- x_filter %>%
      filter(
        time >= tstart, time <= tend
      )
    
    lfit <- lm(log(dist_takens)~time, data=x_filter2)
    
    time_pred <- seq(tstart, 2030, by=0.01)
    
    pp <- predict(lfit, newdata=data.frame(time=time_pred),
                  levels=0.95,
                  interval="prediction")
    
    data.frame(
      pred=exp(pp[,1]),
      pred_lwr=exp(pp[,2]),
      pred_upr=exp(pp[,3]),
      time=time_pred,
      key=x$key[1],
      country=x$country[1],
      tstart=tstart,
      tend=tend
    )
  }
}) %>%
  bind_rows

analysis_all_filter <- analysis_all %>%
  group_by(key, country) %>%
  filter(
    time >= maxt_takens
  )

ggplot(analysis_all) +
  geom_line(aes(time, dist_takens)) +
  # geom_smooth(aes(time, dist_takens), method="loess") +
  geom_ribbon(data=analysis_all_lm, aes(time, ymin=pred_lwr, ymax=pred_upr), fill="red", alpha=0.2) +
  geom_line(data=analysis_all_lm, aes(time, pred), col="red") +
  scale_x_continuous("Year",
                     breaks=seq(2020, 2030, by=2)) +
  scale_y_log10("Distance from attractor") +
  coord_cartesian(xlim=c(2020, 2030), ylim=c(1e-1, 20)) +
  facet_grid(country~key) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
