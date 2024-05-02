bbp <- function(w, k = 7, a = 0, b = 1){
  if(a >= b) stop("The lower bound has to be smaller than the upper bound!")
  q <- sum(w >= a & w <= b)
  v <- seq(a, b, length.out = q)
  vnew <- (v - a)/(b - a)
  basis <- array(0, dim = c(q, k+1))
  for(i in 0:k){
    basis[, i + 1] <- (choose(k, i)) * vnew^i * (1 - vnew)^(k-i)
  }
  return(list("basis" = basis, "angles" = v)) 
}

est_beta <- function(par, basis, t, len_vec, w, lam_end = c(1, 1)){ 
  if(any(lam_end > 1)) stop("Value of the ADF is not supported!") # future me: change text!
  beta <- c(lam_end[1], exp(par), lam_end[2])
  lam <- basis %*% beta
  lam2 <- c() 
  for(i in 1:length(len_vec)){
    lam2[(length(lam2) + 1):(length(lam2) + len_vec[i])] <- rep(lam[i], len_vec[i])
  }
  loglike <- sum(log(lam2)) - sum(lam2 * t)   
  return(-loglike)
}

minfunction_mle <- function(w, data, a = 0, b = 1, lam_end = c(1, 1), k = 7, q = 0.95, tol = 0.0001){
  polynomials <- bbp(w = w, k = k, a = a, b = b)
  basis <- polynomials$basis
  angles <- polynomials$angles
  par_init  <- rep(0, k-1)
  min_proj <- sapply(angles, function(i) minproj_lambda(data = data, w = i, q = q))
  t <- c()
  len_vec <- c()
  for(i in 1:length(angles)){
    aux <- min_proj[1, ][[i]][min_proj[1, ][[i]] > min_proj[2, ][[i]]]
    len_vec[i] <- length(aux)
    t[(length(t) + 1):(length(t) + len_vec[i])] <- aux - min_proj[2, ][[i]]
  }
  results <- tryCatch(optim(par = par_init, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, method = "BFGS", control = list(maxit = 100000)),
                      error = function(e){1})
  if(is.list(results)){optim_output <- results}
  else{
    optim_output <- optim(par = par_init, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
  }
  results <- tryCatch(optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000)),
                      error = function(e){1})
  if(is.list(results)){optim_output2 <- results}
  else{
    optim_output2 <- optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
  }
  
  while(abs(optim_output2$val - optim_output$val) >= tol){
    optim_output <- optim_output2
    results <- tryCatch(optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000)),
                        error = function(e){1})
    if(is.list(results)){optim_output2 <- results}
    else{
      optim_output2 <- optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
    }
  }
  return(c(lam_end[1], exp(optim_output2$par), lam_end[2]))
}