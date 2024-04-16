wads_tawn_curve_exp = function(data_exp,prob,q,w=seq(0.001,0.999,length.out = 200)){#function for estimating return curve for data on standard exponential margins using WT model
  #prob is curve survival probability (small)
  #q is quantile cdf probability for estimating hill estimator of the angular dependence function
  #w is a sequence of w values to fit the model at (spread out in interval [0,1])
  
  #property 3.1
  lim = qexp((1-prob))
  y = lim
  x = 0
  
  for(i in 1:length(w)){
    Min = pmin(data_exp[,1]/w[i],data_exp[,2]/(1-w[i])) #taking the min-projection at each w
    thresh=quantile(Min,q)
    lambda = 1/mean(Min[Min>thresh] - thresh) #hill estimator of the angular dependence function
    if(lambda<max(w[i],1-w[i])){#model constraint
      lambda=max(w[i],1-w[i])
    }
    r=thresh-log(prob/(1-q))/lambda #value corresponding to probability
    x[i+1]=r*w[i] #corresponding x and y values
    y[i+1]=r*(1-w[i])
    
    #imposing property 3.4
    if(x[i+1]<x[i]){
      x[i+1] = x[i]
    }
    if(y[i+1]>y[i]){
      y[i+1] = y[i]
    }
  }
  
  x[i+2] = lim #property 3.1
  y[i+2] = 0
  
  return(cbind(x,y)) #return matrix containing return curve estimate
}

block_bootstrap_function = function(data,k,n=length(as.matrix(data)[,1])){ #function for performing block bootstrapping
  #data is bivariate dataset
  #k is block length
  data = as.matrix(data)
  no_blocks = ceiling(n/k)
  n_new = no_blocks*k
  new_data = matrix(NA,nrow=n_new,ncol=dim(data)[2])
  indices = 1:(n-k+1)
  start_points = sample(x=indices,size=no_blocks,replace=TRUE)
  for(i in 1:no_blocks){
    new_data[((i-1)*k+1):(i*k),] = data[(start_points[i]:(start_points[i]+k-1)),]
  }
  return(new_data)
}

empirical_cdf = function(data,q){ #function for estimating empirical cdf, returns data on uniform margins
  #data is univariate dataset
  #q is marginal quantile level for fitting GPD
  u = c()
  thresh = quantile(data,q)
  par = gpd.fit(data,threshold = thresh,show=FALSE)$mle
  u[data<=thresh] = (rank(data)[data<=thresh])/(length(data)+1)
  u[data>thresh] = 1-(1-q)*pgpd(data[data>thresh],loc=thresh,scale=par[1],shape=par[2],lower.tail = FALSE)
  return(u)
}

curve_inverse_transform = function(vec,odata,q){ #function for moving marginal components of curve estimates on uniform margins back to original
  #vec is vector from return curve on uniform margins
  #odata is original marginal dataset
  #q is marginal quantile level for fitting GPD
  thresh = quantile(odata,q)
  par = gpd.fit(odata,threshold = thresh,show=FALSE)$mle
  nvec = c()
  nvec[vec>q] = qgpd((vec[vec>q]-q)/(1-q),loc=thresh, scale = par[1], shape = par[2])
  nvec[vec<=q] = quantile(odata,vec[vec<=q])
  return(nvec)
}