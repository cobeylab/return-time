library(dplyr)
library(readxl)

year <- 2014:2024

datalist <- vector('list', length(year))

for (i in 1:length(year)) {
  tmp <- read_xlsx("../data_hongkong/hongkong_flu.xlsx", sheet=i)
  
  date <- gsub(".* - ", "",  c(tmp[,2])[[1]])
  date <- gsub(".* -", "", date)
  
  datalist[[i]] <- data.frame(
    year=year[i],
    week=c(tmp[,1])[[1]],
    date=paste0(date, "/", year[i]),
    flua=c(tmp[,3][[1]]),
    flua_h1=c(tmp[,4][[1]]),
    flua_h3=c(tmp[,5][[1]]),
    flub=c(tmp[,9][[1]]),
    flub_vic=c(tmp[,10][[1]]),
    flub_yam=c(tmp[,11][[1]])
  ) 
}

data_hongkong_flu <- datalist %>%
  bind_rows %>%
  mutate(
    date=as.Date(date, format="%d/%m/%Y")
  )

write.csv(data_hongkong_flu, file="data_hongkong_flu.csv")
