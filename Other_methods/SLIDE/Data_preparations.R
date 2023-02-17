data.prep.SLIDE <- function(model) {

  lambda = 0.5
  pvec = model$p
  K = model$K
  n = model$n


  X = t(do.call(rbind, model$X))
  
    
  d <- length(pvec)
  pcum <- cumsum(pvec)
  p.index = as.vector(cumsum(pvec))
  p.index.list = list()
  for(i in 1:d) 
  {
    p.index.list[[i]] = ifelse(i > 1, p.index[i-1] + 1, 1) : p.index[i]
  }
  
  
  
  X.true = matrix(0,n, sum(pvec))
  for (i in 1:K){
    X.true[,p.index.list[[i]]] = model$V[[i]] %*% t(model$U[[i]])
  }
  
  
  X.true.svd = svd(X.true)
  U = X.true.svd$u
  V = t(diag(X.true.svd$d) %*% t(X.true.svd$v))
  
  return(list(X = X, U = U,V = V, pvec=pvec, 
              trueU = model$U, trueV = model$V,
              true.struct = model$true.struct))
}
