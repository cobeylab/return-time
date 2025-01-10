
lm_iterative <- function(dist,
                         time,
                         iter=10) {
  whichmax <- time[which.max(dist)]
  
  tmp <- data.frame(
    time=time,
    dist=dist
  ) %>%
    filter(time >= whichmax)
  
  reslist <- vector('list', iter)
  
  for (i in 1:iter) {
    if (i==1) {
      lfit <- lm(log(dist)~time, data=tmp)
      
      resilience <- -coef(lfit)[[2]]
      chartime <- 1/resilience
      window <- chartime
      
    } else {
      tmp2 <- tmp %>%
        filter(
          time < whichmax + window
        )
      
      lfit <- lm(log(dist)~time, data=tmp2)
      
      resilience <- -coef(lfit)[[2]]
      
      while (resilience < 0) {
        counter <- 2
        
        tmp2 <- tmp %>%
          filter(
            time < whichmax + window*counter
          )
        
        lfit <- lm(log(dist)~time, data=tmp2)
        
        resilience <- -coef(lfit)[[2]]
        
        counter <- counter + 1
      }
      
      chartime <- c(chartime, 1/resilience)
      window <- mean(chartime)
    }
    
    reslist[[i]] <- data.frame(
      iter=i,
      resilience=resilience,
      resilience_lwr=-confint(lfit)[2,2][[1]],
      resilience_upr=-confint(lfit)[2,1][[1]]
    )
  }
  
  out <- reslist %>%
    bind_rows
  
  out
}
