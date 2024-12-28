library(lubridate)
library(readxl)
library(dplyr)
library(tidyr)

nrevss_raw <- read_xlsx("../data_nrevss/NREVSS data Sarah.xlsx", sheet=3)

resp_us_rsv <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$RSVpos,
  positive=nrevss_raw$RSVtest,
  region=nrevss_raw$REGNAME,
  type="RSV"
)

resp_us_piv <- data.frame(
  year=rep(epiyear(nrevss_raw$RepWeekDate), 4),
  week=rep(epiweek(nrevss_raw$RepWeekDate), 4),
  date=rep(nrevss_raw$RepWeekDate, 4),
  tests=rep(nrevss_raw$PIVtest, 4),
  positive=c(
    nrevss_raw$PIV1pos,
    nrevss_raw$PIV2pos,
    nrevss_raw$PIV3pos,
    nrevss_raw$PIV4pos
  ),
  region=rep(nrevss_raw$REGNAME, 4),
  type=rep(c("PIV 1", "PIV 2", "PIV 3", "PIV 4"),
           each=nrow(nrevss_raw))
)

resp_us_piv_comb <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$PIVtest,
  positive=c(
    nrevss_raw$PIV1pos+
    nrevss_raw$PIV2pos+
    nrevss_raw$PIV3pos+
    nrevss_raw$PIV4pos
  ),
  region=nrevss_raw$REGNAME,
  type="PIV"
)

resp_us_adeno <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$RAdenotest,
  positive=nrevss_raw$RAdenopos,
  region=nrevss_raw$REGNAME,
  type="Adenovirus"
)

resp_us_hmpv <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$HMetapneumotest,
  positive=nrevss_raw$HMetapneumopos,
  region=nrevss_raw$REGNAME,
  type="Human metapneumovirus"
)

resp_us_rhino <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$Rhinotest,
  positive=nrevss_raw$Rhinopos,
  region=nrevss_raw$REGNAME,
  type="Rhinovirus"
)

resp_us_cov <- data.frame(
  year=rep(epiyear(nrevss_raw$RepWeekDate), 4),
  week=rep(epiweek(nrevss_raw$RepWeekDate), 4),
  date=rep(nrevss_raw$RepWeekDate, 4),
  tests=rep(nrevss_raw$CoVTest, 4),
  positive=c(
    nrevss_raw$CoVHKU1Pos,
    nrevss_raw$CoV229EPos,
    nrevss_raw$CovOC43Pos,
    nrevss_raw$CovNL63Pos
  ),
  region=rep(nrevss_raw$REGNAME, 4),
  type=rep(c("HKU1", "229E", "OC43", "NL63"),
           each=nrow(nrevss_raw))
)

resp_us_cov_comb <- data.frame(
  year=epiyear(nrevss_raw$RepWeekDate),
  week=epiweek(nrevss_raw$RepWeekDate),
  date=nrevss_raw$RepWeekDate,
  tests=nrevss_raw$CoVTest,
  positive=c(
    nrevss_raw$CoVHKU1Pos+
    nrevss_raw$CoV229EPos+
    nrevss_raw$CovOC43Pos+
    nrevss_raw$CovNL63Pos
  ),
  region=nrevss_raw$REGNAME,
  type="CoV"
)

data_nrevss_resp <- bind_rows(
  resp_us_rsv,
  resp_us_piv,
  resp_us_piv_comb,
  resp_us_adeno,
  resp_us_hmpv,
  resp_us_rhino,
  resp_us_cov,
  resp_us_cov_comb
)

write.csv(data_nrevss_resp, file="data_nrevss_resp.csv")
