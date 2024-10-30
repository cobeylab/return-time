library(lubridate)
library(readxl)
library(dplyr)
library(tidyr)

correct_date <- function(x) {
  out <- sapply(1:length(x), function(i) {
    if (grepl("-", x[i])) {
      as.Date(x[i], format="%d-%m-%Y")
    } else {
      if (i <= 103) {
        tmp <- as.Date(as.numeric(x[i]), origin="1899-12-30")
        
        as.Date(as.character(tmp), format="%Y-%d-%m")
      } else {
        as.Date(as.numeric(x[i]), origin="1899-12-30")
      }
    } 
  }, simplify = FALSE)
  
  as.Date(unlist(out))
}

pathogen <- c("flu", "rsv", "piv", 
              "adeno", "hmp")

data_flu <- read_xlsx("../data_canada/resp_canada.xlsx")

resp_canada_flu <- data.frame(
  week=data_flu$Week,
  date=correct_date(data_flu$`Week End`),
  tests=rep(unlist(data_flu[,grepl("Tests", names(data_flu))]), 2),
  positive=c(unlist(data_flu[,grepl("A%", names(data_flu))]), 
             unlist(data_flu[,grepl("B%", names(data_flu))])),
  province=gsub(" .*", "", names(unlist(data_flu[,grepl("Tests", names(data_flu))]))),
  type=rep(c("Flu A", "Flu B"), each=nrow(data_flu))
)

data_rsv <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=2)

resp_canada_rsv <- data.frame(
  week=data_rsv$Week,
  date=correct_date(data_rsv$`Week End`),
  tests=unlist(data_rsv[,grepl("Tests", names(data_rsv))]),
  positive=unlist(data_rsv[,grepl("RSV%", names(data_rsv))]),
  province=gsub(" .*", "", names(unlist(data_rsv[,grepl("Tests", names(data_rsv))]))),
  type="RSV"
)

data_piv <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=3)

resp_canada_piv <- data.frame(
  week=data_piv$Week,
  date=correct_date(data_piv$`Week End`),
  tests=unlist(data_piv[,grepl("Tests", names(data_piv))]),
  positive=unlist(data_piv[,grepl("Para%", names(data_piv))]),
  province=gsub(" .*", "", names(unlist(data_piv[,grepl("Tests", names(data_piv))]))),
  type="PIV"
)

data_adeno <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=4)

resp_canada_adeno <- data.frame(
  week=data_adeno$Week,
  date=correct_date(data_adeno$`Week End`),
  tests=unlist(data_adeno[,grepl("Tests", names(data_adeno))]),
  positive=as.numeric(gsub(",", ".", unlist(data_adeno[,grepl("Adeno%", names(data_adeno))]))),
  province=gsub(" .*", "", names(unlist(data_adeno[,grepl("Tests", names(data_adeno))]))),
  type="AdV"
)

data_hmpv <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=5)

resp_canada_hmpv <- data.frame(
  week=data_hmpv$Week,
  date=correct_date(data_hmpv$`Week End`),
  tests=as.numeric(gsub(",", ".", unlist(data_hmpv[,grepl("Tests", names(data_hmpv))]))),
  positive=as.numeric(gsub(",", ".", unlist(data_hmpv[,grepl("hMPV%", names(data_hmpv))]))),
  province=gsub(" .*", "", names(unlist(data_hmpv[,grepl("Tests", names(data_hmpv))]))),
  type="HMPV"
)

data_rhino <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=6)

resp_canada_rhino <- data.frame(
  week=data_rhino$Week,
  date=correct_date(data_rhino$`Week End`),
  tests=unlist(data_rhino[,grepl("Tests", names(data_rhino))]),
  positive=unlist(data_rhino[,grepl("RhV%", names(data_rhino))]),
  province=gsub(" .*", "", names(unlist(data_rhino[,grepl("Tests", names(data_rhino))]))),
  type="RV"
)

data_entero <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=7)

resp_canada_entero <- data.frame(
  week=data_entero$Week,
  date=correct_date(data_entero$`Week End`),
  tests=as.numeric(unlist(data_entero[,grepl("Tests", names(data_entero))])),
  positive=as.numeric(unlist(data_entero[,grepl("Entero/rhino%", names(data_entero))])),
  province=gsub(" .*", "", names(unlist(data_entero[,grepl("Tests", names(data_entero))]))),
  type="RV/EV"
)

data_cov <- read_xlsx("../data_canada/resp_canada.xlsx", sheet=8)

resp_canada_cov <- data.frame(
  week=data_cov$Week,
  date=correct_date(data_cov$`Week End`),
  tests=unlist(data_cov[,grepl("Tests", names(data_cov))]),
  positive=unlist(data_cov[,grepl("Coro%", names(data_cov))]),
  province=gsub(" .*", "", names(unlist(data_cov[,grepl("Tests", names(data_cov))]))),
  type="CoV"
)

data_canada_resp <- bind_rows(
  resp_canada_flu,
  resp_canada_rsv,
  resp_canada_piv,
  resp_canada_adeno,
  resp_canada_hmpv,
  resp_canada_rhino,
  resp_canada_entero,
  resp_canada_cov
)

rownames(data_canada_resp) <- NULL

write.csv(data_canada_resp, file="data_canada_resp.csv")
