distfun <- function(mat_perturb,
                    mat_unperturb,
                    time, tbreak,
                    out="comb") {
  mat_sd <- apply(mat_unperturb, 2, sd)

  mat_perturb_scale <- t(t(mat_perturb)/mat_sd)
  mat_unperturb_scale <- t(t(mat_unperturb)/mat_sd)
    
  mat_dist <- abs(mat_perturb_scale-mat_unperturb_scale)
  
  if (out=="comb") {
    sqrt(rowSums(mat_dist^2))
  } else if (out=="all") {
    out <- as.data.frame(mat_dist)
    
    out$dist <- sqrt(rowSums(mat_dist^2))
    out$time <- time
    
    out
  }
}

distfun_nn <- function(mat_perturb,
                       mat_unperturb) {
  sapply(1:nrow(mat_perturb), function(i) {
    dd <- sqrt(colSums((mat_perturb[i,] - t(mat_unperturb))^2))
    
    min(dd)
  })  
}
