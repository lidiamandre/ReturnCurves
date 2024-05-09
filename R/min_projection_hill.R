minproj_lambda <- function(data, w = seq(0, 1, by = 0.01), q = 0.95){
  t <- pmin(data[, 1]/w, data[, 2]/(1 - w))
  u <- quantile(t, q)
  lambda <- 1/(mean(t[t > u] - u))
  return(list("minproj" = t, "thresh" = u, "lambdahill" = lambda))
}