distfun <- function(mat_perturb,
                    mat_unperturb,
                    time, tbreak,
                    out="comb") {
  mat_sd <- apply(mat_perturb[time<tbreak,], 2, sd)

  mat_perturb_scale <- t(t(mat_perturb)/mat_sd)
  mat_unperturb_scale <- t(t(mat_unperturb)/mat_sd)
    
  mat_dist <- mat_perturb_scale-mat_unperturb_scale
  
  if (out=="comb") {
    sqrt(rowSums(mat_dist^2))
  } else if (out=="all") {
    out <- as.data.frame(mat_dist)
    
    out$dist <- sqrt(rowSums(mat_dist^2))
    out$time <- time
    
    out
  }
}
