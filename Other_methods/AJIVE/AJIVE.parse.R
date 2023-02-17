source("./Other_methods/AJIVE/ajive.R")


AJIVE.parse <- function(X, n, r.list, p.list, p, K){
  
  
  AJIVE.result = ajive(X, r.list)
  
  
  
  r.joint = AJIVE.result$block_decomps[[1]]$joint$rank
  
  r.ind = c()
  for (i in 1:K) {
    r.ind = c(r.ind, AJIVE.result$block_decomps[[i]]$individual$rank)
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

  
  
  AJIVE.loading = matrix(0, sum(p), sum(r.ind) + r.joint)
  AJIVE.score = matrix(0, n, sum(r.ind) + r.joint)
  AJIVE.struct = matrix(0, K, sum(r.ind) + r.joint)
  
  if (r.joint !=0 ) {
    AJIVE.score[,1:r.joint] = AJIVE.result$joint_scores
    for(i in 1:K) {
      # Joint.diag = AJIVE.result$block_decomps[[i]]$joint$d[1:r.joint]
      # Joint.Diag = diag(Joint.diag, r.joint)
      # AJIVE.loading[p.list[[i]], 1:r.joint] = AJIVE.result$block_decomps[[i]]$joint$v %*% Joint.Diag
      
      AJIVE.struct[i, 1:r.joint] = 1
      AJIVE.loading[p.list[[i]], 1:r.joint] = t(AJIVE.result$block_decomps[[i]]$joint$full) %*% AJIVE.result$joint_scores
    }
    
  }
  
  for(i in 1:K) {
    if (r.ind[i] != 0) {
      Ind.diag = AJIVE.result$block_decomps[[i]]$individual$d[1:r.ind[i] ]
      Ind.Diag = diag(Ind.diag, r.ind[i])
      AJIVE.loading[p.list[[i]], r.ind.list[[i]]] = AJIVE.result$block_decomps[[i]]$individual$v %*% Ind.Diag
      AJIVE.struct[i, r.ind.list[[i]]] = 1
      AJIVE.score[,r.ind.list[[i]]] = AJIVE.result$block_decomps[[i]]$individual$u
    }
  }
  
  return(list(AJIVE.loading = AJIVE.loading, AJIVE.score = AJIVE.score,AJIVE.struct = AJIVE.struct))
  
}