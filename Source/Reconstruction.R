



# FUNCTION mode.structure(result.parallel)
#
# find the modal structure from a list of partially-joint structures
#
# INPUT  
#   result.parallel : the results from estimation.fun (see Estimation_fun.R)
#
# OUTPUT  
#   result.optimal : the list of structure outputs
#   model.mode : the index of the mode on result.optimal
#   index.min.r : the index that has the least penalized risk within the mode structure (the index for optimal estimate)



mode.structure = function(result.parallel) {
  table.structure = list()
  table.index = list()
  str.count = c()
  
  table.structure[[1]] = result.parallel[[1]]$structure
  table.index[[1]] = c(1)
  str.count = c(1)
  
  
  count = 1
  nrep = length(result.parallel)
  
  for (j in 2:nrep) {
    rep.check = 0
    for (k in 1:count) {
      if (  identical(result.parallel[[j]]$structure, table.structure[[k]]) ) {
        str.count[k] = str.count[k] + 1
        table.index[[k]]= c(table.index[[k]] ,j)
        rep.check = 1
        break
      }
    }
    if (rep.check == 0){
      count = count + 1
      table.structure[[count]] = result.parallel[[j]]$structure
      str.count[count] = 1
      table.index[[count]] = c(j)
    }
  }
  
  
  model.mode = which.max(str.count)
  result.mode = list(structure = table.structure, count = str.count, index = table.index, model.mode = model.mode)
  
  
  
  
  
  
  return(result.mode)
}





# FUNCTION reconstruction.partial.cluster.mode(model.list, result.parallel, index.min.r, scale.param)
#
# reconstruction of partially-joint components on the mode structure
#
# INPUT  
#   model.list : list of parameters (n, K, p, r.list, X)
#   result.parallel : the result from estimation.fun (see Estimation_fun.R)
#   index.min.r : the index that has the least penalized risk within the mode structure
#   scale.param : scale attribute of each dataset
#
# OUTPUT  
#   pc : the list of principal directions from the optimal partial clustering structure estimate
#   pc.index : the number each principal direction appear on the optimal structure
#   pc.matrix : the reconstructed signal for each parincipal direcitions (rank can be over 1)





reconstruction.partial.cluster.mode = function(model.list, mode.str, scale.param, lambda) {
  
  
  r.list = model.list$r.list
  n = model.list$n
  K = model.list$K
  p = model.list$p
  
  
  V.list.untruncated = list()
  V.list = list()
  

  
  for (i in 1:K ) {
    V.list.untruncated[[i]] = svd(model.list$X[[i]])$v
    rv = min(nrow(V.list.untruncated[[i]]),ncol(V.list.untruncated[[i]]),r.list[i])
    if (rv > 0) {
      V.list[[i]] = V.list.untruncated[[i]][,1:rv  ]
    } else {
      V.list[[i]] = V.list.untruncated[[i]][,0]
    }
  }
  
  result.whole = list()
  result.str.whole = list()
  result.U.whole = list()
  result.diff.whole = rep(0,length(lambda))
  for (i in 1:length(lambda)) {
    result.whole[[i]] <- main.fun(V.list, n, K, r.list, lambda = lambda[i])  
    result.str.whole[[i]] <- get.structure(result.whole[[i]], K)
    result.U.whole[[i]] <- U.estimate(K, p, result.str = result.str.whole[[i]], 
                                      subspace = result.whole[[i]]$subspace, model.list$X)
    result.diff.whole[i] = distance.struct(mode.str,result.str.whole[[i]]$structure)
  }
  
  
  lam.hat = which.min(result.diff.whole)
  
  estim.str = result.str.whole[[lam.hat]]$structure
  
  

  
  pc = list()
  pc.count = c()
  pc.index = list()
  
  pc[[1]] = estim.str[,1]
  pc.count = c(1)
  pc.index[[1]] = c(1)
  
  count = 1
  
  for (j in 2:ncol(estim.str)) {
    rep.check = 0
    for (k in 1:count) {
      if (  identical(estim.str[,j], pc[[k]]) ) {
        pc.count[k] = pc.count[k] + 1
        pc.index[[k]] = c(pc.index[[k]] ,j)
        rep.check = 1
        break
      }
    }
    if (rep.check == 0){
      count = count + 1
      pc[[count]] = estim.str[,j]
      pc.count[count] = 1
      pc.index[[count]] = c(j)
    }
  }
  
  
  
  K = model.list$K
  n = model.list$n
  p = model.list$p
  
  p.cumsum = cumsum(p)
  p.list = list()
  p.list[[1]] = 1:p.cumsum[1]
  if (K > 1) {
    for (i in 2:K) {
      p.list[[i]] = (p.cumsum[i-1] + 1) : (p.cumsum[i])  
    }
  }
  
  
  U = result.U.whole[[lam.hat]]$U
  V = result.U.whole[[lam.hat]]$V
  
  
  pc.matrix = list()
  pc.split = list()
  for(i in 1:count) {
    index = pc.index[[i]]
    mat = U[,pc.index[[i]]] %*% t(V[,pc.index[[i]]])
    for(j in 1:K) { 
      pc.split[[j]] = t(mat[p.list[[j]],] )
    }
    pc.matrix[[i]] = pc.split
  }
  
  
  return (list(pc = pc,pc.index = pc.index,pc.matrix = pc.matrix, opt.lambda = lam.hat) )
}

