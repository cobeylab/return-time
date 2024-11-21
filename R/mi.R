mi <- function(x,
               lag.max=208,
               bin=5) {
  ii <- sapply(1:lag.max, function(z) {
    mi_internal(x, lag=z, bin=bin)
  })
  
  ii
}

mi_internal <- function(x, lag,
                        bin=10) {
  xmin <- min(x)
  xmax <- max(x)
  
  xseq <- seq(xmin, xmax, length.out=bin+1)
  
  xcut <- head(x, -lag)
  xlag <- tail(x, -lag)
  
  xcut_group <- as.numeric(cut(xcut, xseq, include.lowest = TRUE))
  xlag_group <-as.numeric(cut(xlag, xseq, include.lowest = TRUE))
  
  P_hk <- as.matrix(table(xcut_group, xlag_group))
  P_hk <- P_hk/sum(P_hk)
  
  P_h <- c(table(xcut_group))
  P_h <- P_h/sum(P_h)
  
  P_k <- c(table(xlag_group))
  P_k <- P_k/sum(P_k)
  
  P_all <- P_hk * log(t(t(P_hk/P_k)/P_h))
  
  -sum(P_all[is.finite(P_all)])
}
