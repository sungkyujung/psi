
# FUNCTION main.fun(V.list, n, K, r.list, lambda = 0.01)
#
# the main algorithm
#
# INPUT  
#   V.list : the list of score matrices from datasets
#   n : the number of observations
#   K : the number of datasets
#   r.list : the vector of ranks to each dataset
#   ind.perm : index-set permutation (length = 2^K - 1)
#
# OUTPUT  
#   cluster : estimated structure
#   subspace : estimated partially joint scores
#   ranks : ranks of fully or partially-joint/individual/whole index-sets

main.fun <- function(V.list, n, K, r.list, lambda = 0.01, ind.perm)
{
  
  
  ncol.list = rep(0,K)
  for(i in  1:K) {
    V.list[[i]] = as.matrix(V.list[[i]])
    ncol.list[i] = ncol(V.list[[i]])
  }
  


  
  # Generating the power sets.
  power.set.list = power.set(K)
  
  # Shuffle the ordering of index-sets
  power.set.list = power.set.list[,ind.perm]
  
  
  # variables for the result of the algorithm
  clus.result = list(); clus.result.num = 0
  clus.subspace = list(); clus.subspace.num = 0
  clus.rank = rep(0,ncol(power.set.list))
  
  
  
  
  
  for (l in 1:(ncol(power.set.list) - K)){
    
    # get the current index-set (current.cluster)
    current.cluster = power.set.list[,l]
    
    while(TRUE){
      
      # check whether all the V spaces in the index-set are not empty
      if (!identical(which( ncol.list[current.cluster] == 0) < 1, logical(0))) {break}
      
      # STEP (a) : calculate the mean direction (current.prin.score) for the current step
      
      minr = min(ncol.list[current.cluster])
      V = do.call(cbind,V.list[current.cluster]) 
      
      
      current.prin.score = SVDavr(V)
      
      # STEP (b) : check whether the angle between mean direction and each V space are under the threshold
      edge.check = 0
      for(k in 1:length(current.cluster)) {
        if (current.cluster[k]){
          edge.check = max(edge.check,  acos(  pmin( pmax( 
            svd( t(current.prin.score) %*% V.list[[k]])$d[1] ,-1),1))  )
        }
      }
      if (edge.check > lambda) {break}
      
      
      # record the index-set and the mean direction
      clus.result.num = clus.result.num + 1
      clus.result[[clus.result.num]] = current.cluster
      clus.subspace.num = clus.subspace.num + 1
      clus.subspace[[clus.subspace.num]] = current.prin.score
      clus.rank[l] = clus.rank[l] + 1
      
      
      # STEP (c) : truncate each \widehat{V}_k by a Householder transformation
      for(k in 1:(length(current.cluster))){
        if (current.cluster[k]){
          current.subspace.projected = proj.QR(V.list[[k]], current.prin.score)
          V.list.tmp = as.matrix(Householder(V.list[[k]], current.subspace.projected))
          
          if (ncol(V.list.tmp) > 1 ){
            V.list[[k]] = as.matrix(V.list.tmp[, 2:ncol(V.list.tmp)])
            ncol.list[k] = ncol.list[k] - 1
          } else if (ncol(V.list.tmp) <= 1) {
            V.list[[k]] =  matrix(0,nrow(V.list[[k]]),0)
            if (ncol(V.list.tmp) == 1) {ncol.list[k] = ncol.list[k] - 1}
          } 
        }
      }
      
    }
  }
  
  # the sum of ranks from (partial) joint structures
  joint.rank = clus.result.num
  total.rank = clus.result.num
  
  # record the remaining individual scores
  
  for (i in 1:K){
    if (ncol(V.list[[i]]) > 0) {
      clus.result.num = clus.result.num + 1
      clus.result[[clus.result.num]] = c(i)
      
      clus.subspace.num = clus.subspace.num + 1
      clus.subspace[[clus.subspace.num]] = V.list[[i]]
      
      clus.rank[ncol(power.set.list) - K + i] = clus.rank[ncol(power.set.list) - K + i] + 1
      
      total.rank = total.rank + ncol(V.list[[i]])
    }
  }
  
  
  
  result = list(cluster = clus.result, subspace = clus.subspace, ranks = c(joint.rank, clus.result.num, total.rank))
  
  return (result)


}
