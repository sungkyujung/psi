


source("./Other_methods/JIVE/JIVE.R")

JIVE.parse <- function(model, X,  method){
  


  # apply JIVE
  JIVE.result = jive(data = X, method=method,  showProgress=FALSE) 

  r.joint = JIVE.result$rankJ
  r.A = JIVE.result$rankA
  
  
  p.list = model$p.list
  K = model$K
  
  
  r.ind = c()
  for (i in 1:K) {
    r.ind = c(r.ind, r.A[i])
  }
  
  r.ind.list = list()
  r.cumsum = cumsum(r.ind)
  
  if (r.ind[[1]] > 0) {
    r.ind.list[[1]] = (1:r.cumsum[1]) + r.joint
  } else {
    r.ind.list[[1]] = 0
  }
  
  if (K > 1) {
    for (i in 2:K) {
      if (r.ind[i] > 0) {
        r.ind.list[[i]] = ( (r.cumsum[i-1] + 1) : (r.cumsum[i]) ) + r.joint 
      } else {
        r.ind.list[[i]] = 0
      }
    }
  }
  
  JIVE.loading = matrix(0, sum(model$p),r.joint + sum(r.ind) )
  JIVE.score = matrix(0, model$n, r.joint + sum(r.ind) )
  JIVE.structure = matrix(0, K, r.joint + sum(r.ind))
  
  
  sc = JIVE.result$scale$`Scale Values`
  for(k in 1:K) {
    JIVE.result$joint[[k]] = JIVE.result$joint[[k]] * sc[k]
    JIVE.result$individual[[k]] = JIVE.result$individual[[k]] * sc[k]
  }
  
  
  if (r.joint > 0) {
    joint.svd = svd(do.call(rbind, JIVE.result$joint))
    joint.d = joint.svd$d[1:r.joint]
    joint.D = diag(joint.d, r.joint)
    
    JIVE.loading[,1:r.joint] = joint.svd$u[,1:r.joint] %*% joint.D
    JIVE.score[,1:r.joint] = joint.svd$v[,1:r.joint]
    JIVE.structure[,1:r.joint] = 1
  }
  
  for(i in 1:K) {
    if (length(r.ind.list[[i]]) > 0) {
      ind.svd = svd(JIVE.result$individual[[i]])
      ind.d = ind.svd$d[1:r.A[i]]
      ind.D = diag(ind.d, r.A[i])
      
      JIVE.loading[p.list[[i]], r.ind.list[[i]]] = ind.svd$u[,1:(r.A[i])] %*% ind.D
      JIVE.score[ , r.ind.list[[i]]] = ind.svd$v[,1:(r.A[i])]
      JIVE.structure[i,r.ind.list[[i]]] = 1
    }
  }
  
  
  

  return(list(JIVE.loading = JIVE.loading, JIVE.score = JIVE.score, JIVE.struct = JIVE.structure))
}


