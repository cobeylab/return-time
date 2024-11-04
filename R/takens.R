takens <- function(x,
                   d=3,
                   tau=1) {
  out <-  matrix(0, nrow=length(x)-tau*(d-1), ncol=d)
  
  for (i in 1:d) {
    if (i==1) {
      coord <- head(x, -tau*(d-1))
    } else if (i==d) {
      coord <- tail(x, -tau*(d-1))
    } else {
      coord <- tail(head(x, -tau*(d-i)), -tau*(i-1))
    }
    
    out[,i] <- coord
  }
  out
}
