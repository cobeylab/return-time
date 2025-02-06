library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times"))
library(egg)
library(gridExtra)
source("../R/simulate_sirs.R")
source("../R/simulate_sirs_stoch.R")
source("../R/npi_random.R")
source("../R/takens.R")

npifun_random <- npifun_random_generate(duration=4,
                                        npimin=0.5,
                                        seed=11)

npidata <- data.frame(
  time=seq(2015, 2030, by=1/52/7),
  npi=sapply(seq(2015, 2030, by=1/52/7), npifun_random)
)

ss <- simulate_SIRS_stoch(
  R0=2,
  theta=0.1,
  delta=1/52/7/1,
  npifun=npifun_random)

ee <- eigen_sirs(b_1=2*(365/7+1/50),
                 mu=1/50,
                 gamma=365/7,
                 delta=1/1)

ss_weekly <- ss %>%
  group_by(year, week) %>%
  summarize(
    cases=sum(cases)
  ) %>%
  mutate(
    time=year+week/52
  ) %>%
  filter(
    year > 2014, year < 2030
  )

logcases_pre <- log(filter(ss_weekly, year < 2020)$cases+1)

acfout <- acf(logcases_pre, lag.max=(52*2), plot=FALSE)

tau <- which(head(acfout$acf, -1) > 0 & tail(acfout$acf, -1) < 0)[1]

n.fnn <- fnn(logcases_pre, dmax=25, tau=tau, R_tol=10)

d <- which(n.fnn==0)[1]+1

tauvec <- 10:15
dvec <- 2:6

paramdata <- expand.grid(tauvec, dvec)
reslist <- vector('list', nrow(paramdata))

for (i in 1:nrow(paramdata)) {
  pp <- paramdata[i,]
  taunew <- pp[[1]]
  dnew <- pp[[2]]
  
  takens_perturb <- takens(log(ss_weekly$cases+1), d=dnew, tau=taunew)
  takens_unperturb <- takens(logcases_pre, d=dnew, tau=taunew)
  
  takens_data <- as.data.frame(takens_perturb)
  takens_data$time <- tail(ss_weekly$year+ss_weekly$week/52, -taunew*(dnew-1))
  
  dist <- sapply(1:nrow(takens_perturb), function(i) {
    dd <- sqrt(colSums((takens_perturb[i,] - t(takens_unperturb))^2))
    
    min(dd[dd>0])
  })
  
  reslist[[i]] <- data.frame(
    time=tail(ss_weekly$year+ss_weekly$week/52, -taunew*(dnew-1)),
    dist=dist,
    d=dnew,
    tau=taunew
  )
}

resdata <- reslist %>%
  bind_rows %>%
  mutate(
    tau=paste0(tau, " days"),
    d=paste0(d, " dimensions")
  )

g1 <- ggplot(resdata) +
  annotate("rect", xmin=2020, xmax=2024, ymin=0.08, ymax=Inf, fill="gray", alpha=0.4) +
  geom_point(aes(time, dist), shape=1, size=0.5) +
  scale_x_continuous("Year", limits=c(2015, 2030.5), expand=c(0, 0)) +
  scale_y_log10("Distance from attractor") +
  facet_grid(d~tau, scale="free") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    axis.text.x = element_text(angle=45, hjust=1)
  )

ggsave("figure3_sens.pdf", g1, width=8, height=6)