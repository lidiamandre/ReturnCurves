HeffTawnNegLL <- function(X, Y, par){ 
  alpha <- par[1]
  beta <- par[2]
  sig <- par[3]
  mu <- par[4]
  if(alpha < 0 || alpha > 1 || beta > 1 || beta < 0 || sig <= 0){
    return(1e10)
  }
  negloglik <- -sum(dnorm(Y, alpha*X + mu*((X)^beta), sig*((X)^beta), log = T))
  if(is.finite(negloglik)){
    return(negloglik)
  }
  else{
    return(1e10)
  }
}

heff_tawn_alphas <- function(data, q = 0.95){
  u <- apply(data, 2, quantile, probs = q)
  excdata <- sapply(1:dim(data)[2], function(i) data[data[, i] > u[i], ], simplify = F)
  par <- rep(1/2, 4)
  Yopt <- optim(fn = HeffTawnNegLL, X = excdata[[2]][, 2], Y = excdata[[2]][, 1], par = par, control = list(maxit=100000)) 
  Ypar <- Yopt$par
  Xopt <- optim(fn = HeffTawnNegLL, X = excdata[[1]][, 1], Y = excdata[[1]][, 2], par = par, control = list(maxit=100000))
  Xpar <- Xopt$par
  return(c(Ypar[1], Xpar[1]))
}