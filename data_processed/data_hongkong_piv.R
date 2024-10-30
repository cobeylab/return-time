library(dplyr)
library(readxl)

year <- 2014:2024

datalist <- vector('list', length(year))

for (i in 1:length(year)) {
  tmp <- read_xlsx("../data_hongkong/hongkong_parainfluenza.xlsx", sheet=i, skip=1)
  
  date <- gsub(".* - ", "",  c(tmp[,2])[[1]])
  date <- gsub(".* -", "", date)
  
  if (year[i] %in% 2020:2022) {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      piv1=c(tmp[,5])[[1]],
      piv2=c(tmp[,7])[[1]],
      piv3=c(tmp[,9])[[1]],
      piv4=c(tmp[,11])[[1]],
      test=c(tmp[,4])[[1]],
      age=c(tmp[,3])[[1]]
    ) 
  } else {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      piv1=c(tmp[,4])[[1]],
      piv2=c(tmp[,6])[[1]],
      piv3=c(tmp[,8])[[1]],
      piv4=c(tmp[,10])[[1]],
      test=c(tmp[,3])[[1]],
      age=NA
    ) 
  }
}

data_hongkong_piv <- datalist %>%
  bind_rows %>%
  mutate(
    date=as.Date(date, format="%d/%m/%Y")
  )

write.csv(data_hongkong_piv, file="data_hongkong_piv.csv")
