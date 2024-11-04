library(dplyr)
library(readxl)

year <- 2013:2024

datalist <- vector('list', length(year))

for (i in 1:length(year)) {
  tmp <- read_xlsx("../data_hongkong/hongkong_fecal.xlsx", sheet=i, skip=1)
  
  date <- gsub(".* - ", "",  c(tmp[,2])[[1]])
  date <- gsub(".* -", "", date)
  
  if (year[i] <= 2016) {
    datalist[[i]] <- data.frame(
      year=year[i],
      month=c(tmp[,1])[[1]],
      noro=c(tmp[,3])[[1]],
      rota=c(tmp[,5])[[1]],
      test=c(tmp[,2])[[1]]
    ) 
  } else {
    datalist[[i]] <- data.frame(
      year=year[i],
      month=c(tmp[,1])[[1]],
      noro=c(tmp[,3])[[1]],
      rota=c(tmp[,5])[[1]],
      astro=as.numeric(c(tmp[,7])[[1]]),
      sapo=as.numeric(c(tmp[,9])[[1]]),
      adeno=as.numeric(c(tmp[,11])[[1]]),
      test=c(tmp[,2])[[1]]
    ) 
  }
}

data_hongkong_fecal <- datalist %>%
  bind_rows

write.csv(data_hongkong_fecal, file="data_hongkong_fecal.csv")
