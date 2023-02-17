


source("./Other_methods/BIDIFAC.plus/bidifac.plus.R")
source("./Other_methods/BIDIFAC.plus/bidifac.plus.data.prep.R")


bidifac.plus.parse <- function(model) {
  res <- bidifac.plus(X0=model$X,p.ind=model$p.ind,n.ind=model$n.ind)
  
  
  if(is.null(res$p.ind.list)) {return ( list(bidifac.plus.loading = NA, 
                                           bidifac.plus.score = NA,
                                           bidifac.plus.struct = NA)   )}
  
  K = tail(res$p.ind[[length(res$p.ind)]],1)
  n = tail(res$n.ind[[length(res$n.ind)]],1)
  
  n.module = length(res$p.ind.list)
  
  structure = matrix(0,K,0)
  U = matrix(0,p,0)
  V = matrix(0,n,0)
  
  if (!is.na(res$Sums)){
    
    for(i in 1:n.module){
      svd.x = svd(res$S[[i]])
      d.temp = svd.x$d
      r.hat = length(which(d.temp>1e-10))
      if (r.hat > 0) {
        # structure
        col.temp = matrix(0,K,1)
        for(j in 1:K) {
          if (p.ind[[j]][1] %in% res$p.ind.list[[i]]) {col.temp[j,1] = 1}
        }
        for(j in 1:r.hat) {
          structure = cbind(structure, col.temp)
        }
        U = cbind(U,svd.x$u[,1:r.hat])
        V = cbind(V,svd.x$v[,1:r.hat])
      }  
      
    }
    
  }  
    return(list(bidifac.plus.loading = U, 
                bidifac.plus.score = V,
                bidifac.plus.struct = structure))
   
  
}