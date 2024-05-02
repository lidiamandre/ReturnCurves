block_bootstrap_function <- function(data, k = 1, n = dim(data)[1]){ 
  data <- as.matrix(data)
  no_blocks <- ceiling(n/k)
  n_new <- no_blocks*k
  new_data <- matrix(NA, nrow = n_new, ncol = dim(data)[2])
  indices <- 1:(n-k+1)
  start_points <- sample(x = indices, size = no_blocks, replace = T)
  for(i in 1:no_blocks){
    new_data[((i-1)*k+1):(i*k), ] = data[(start_points[i]:(start_points[i] + k - 1)), ]
  }
  return(new_data[1:n, ])
} 