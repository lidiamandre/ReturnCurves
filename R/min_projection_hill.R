minproj_lambda <- function(data, w, q_minproj){
  t <- pmin(data[, 1]/w, data[, 2]/(1 - w))
  u <- quantile(t, q_minproj)
  lambda <- 1/(mean(t[t > u] - u))
  return(list("minproj" = t, "thresh" = u, "lambdahill" = lambda))
}