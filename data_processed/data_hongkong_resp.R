library(dplyr)
library(readxl)

year <- 2014:2024

datalist <- vector('list', length(year))

for (i in 1:length(year)) {
  tmp <- read_xlsx("../data_hongkong/hongkong_other_resp.xlsx", sheet=i, skip=1)
  
  date <- gsub(".* - ", "",  c(tmp[,2])[[1]])
  date <- gsub(".* -", "", date)
  
  if (year[i] <= 2017) {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      adeno=c(tmp[,4])[[1]],
      rsv=c(tmp[,6])[[1]],
      rvev=c(tmp[,8])[[1]],
      rv=c(tmp[,11])[[1]],
      ev=c(tmp[,12])[[1]],
      test=c(tmp[,3])[[1]],
      age=NA
    ) 
  } else if (year[i] < 2020) {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      adeno=c(tmp[,4])[[1]],
      rsv=c(tmp[,6])[[1]],
      hmpv=c(tmp[,8])[[1]],
      rvev=c(tmp[,10])[[1]],
      rv=c(tmp[,13])[[1]],
      ev=c(tmp[,14])[[1]],
      test=c(tmp[,3])[[1]],
      age=NA
    )
  } else if (year[i] %in% 2020:2022) {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      adeno=c(tmp[,5])[[1]],
      rsv=c(tmp[,7])[[1]],
      hmpv=c(tmp[,9])[[1]],
      rvev=c(tmp[,11])[[1]],
      rv=c(tmp[,14])[[1]],
      ev=c(tmp[,15])[[1]],
      test=c(tmp[,4])[[1]],
      age=c(tmp[,3])[[1]]
    )
  } else {
    datalist[[i]] <- data.frame(
      year=year[i],
      week=c(tmp[,1])[[1]],
      date=paste0(date, "/", year[i]),
      adeno=c(tmp[,4])[[1]],
      rsv=c(tmp[,6])[[1]],
      hmpv=c(tmp[,8])[[1]],
      rvev=c(tmp[,10])[[1]],
      rv=c(tmp[,13])[[1]],
      ev=c(tmp[,14])[[1]],
      test=c(tmp[,3])[[1]],
      age=NA
    ) 
  }
}

data_hongkong_resp <- datalist %>%
  bind_rows %>%
  mutate(
    date=as.Date(date, format="%d/%m/%Y")
  )

write.csv(data_hongkong_resp, file="data_hongkong_resp.csv")
