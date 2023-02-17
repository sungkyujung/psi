data.prep.JIVE <- function(model, sep = FALSE, r.list = NULL) {
  
  lambda = 0.5
  pvec = model$p
  K = model$K
  n = model$n
  
  
  d <- length(pvec)
  pcum <- cumsum(pvec)
  p.index = as.vector(cumsum(pvec))
  p.index.list = list()
  for(i in 1:d) 
  {
    p.index.list[[i]] = ifelse(i > 1, p.index[i-1] + 1, 1) : p.index[i]
  }
  
  
  
  
  if (is.null(r.list)) {
    X = model$X
  } else if (sep) {
    ### initial ranks for X_k separately
    X.temp = list()
    for(i in 1:model$K) {
      svd.X = svd(model$X[[i]])
      svd.X$d[(r.list[i] + 1):length(svd.X$d)] = 0
      X.temp[[i]] = svd.X$u %*% diag(svd.X$d) %*% t(svd.X$v)
    }
    X = X.temp
    svd.X = NULL
  } else {
    ### initial rank for row-concatenated X
    X = do.call(rbind, model$X)
    svd.X = svd(X)
    svd.X$d[(r.list + 1):length(svd.X$d)] = 0
    X = svd.X$u %*% diag(svd.X$d) %*% t(svd.X$v)
    
    
    X.out = list()
    for(k in 1:K) {
      X.out[[k]] = as.matrix(X[ p.index.list[[k]] ,])
    }
      
    X = X.out
    svd.X = NULL
  }
  
  
  
  
  
  return(X)
}
