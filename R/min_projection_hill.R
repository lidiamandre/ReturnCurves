minproj_lambda <- function(data, w = seq(0, 1, by = 0.01), q = 0.95){
  Q <- pmin(data[, 1]/w, data[, 2]/(1 - w))
  u <- quantile(Q, q)
  lambda <- 1/(mean(Q[Q > u] - u))
  return(list("minproj" = Q, "thresh" = u, "lambdahill" = lambda))
}