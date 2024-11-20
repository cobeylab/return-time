takens <- function(x,
                   d,
                   tau) {
  out <-  matrix(0, nrow=length(x)-tau*(d-1), ncol=d)
  
  for (i in 1:d) {
    if (i==1) {
      coord <- tail(x, -tau*(d-1))
    } else if (i==d) {
      coord <- head(x, -tau*(d-1))
    } else {
      coord <- head(tail(x, -tau*(d-i)), -tau*(i-1))
    }
    
    out[,i] <- coord
  }
  out
}

fnn_internal <- function(x,
                         d,
                         tau,
                         R_tol=10) {
  tmp_takens <- takens(x, d=d, tau=tau)
  
  nn_index <- c(sapply(1:nrow(tmp_takens), function(i) {
    dist <- sqrt(colSums((tmp_takens[i,] - t(tmp_takens))^2))
    
    which(rank(dist, ties.method="first")==2)
  }))
  
  dist_d <- sqrt(rowSums((tmp_takens - tmp_takens[nn_index,])^2))
  
  tmp_takens2 <- takens(x, d=d+1, tau=tau)
  
  nn_index2 <- tail(nn_index-tau, -tau)
  nn_index2[nn_index2 < 1] <- NA
  
  dist_d2 <- sqrt(rowSums((tmp_takens2 - tmp_takens2[nn_index2,])^2))
  
  R <- dist_d2/tail(dist_d, -tau)
  
  sum(R > R_tol, na.rm=TRUE)
}

fnn <- function(x, 
                dmax=10,
                tau,
                R_tol=10) {
  if (dmax*tau > length(x)) {
    dmax <- floor(length(x)/tau)
  }
  
  n.fnn <- sapply(2:dmax, function(d) {
    fnn_internal(x, d, tau, R_tol)
  })
  
  n.fnn
}
