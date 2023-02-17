



# FUNCTION estimation.fun(model, lambda, set.training, light.ver = FALSE)
#
# estimate the optimal structure from the training set that minimizes the penalized risk
#
# INPUT  
#   model : a list of model data sets
#   lambda : a grid of lambda values
#   set.training : splitting of n observations (indices for the training set)
#   r.list : the ranks of the estimated signal matrices 
#   index.perm : permutated index-set order. if null, generate one.
#   light.ver : return only optimal structure and risk values when TRUE
#
# OUTPUT  
# LIGHT VERSION
#   structure : the estimated structure from the optimal lambda
#   r.test.list : a vector of risk values
#   r.pen.test.list : a vector of penalized risk values
#   r.list : a vector of (estimated) signal ranks
# FULL VERSION
#   r.test : risk from the optimal lambda 
#   r.pen.test : penalized risk from the optimal lambda 
#   result.U : loading from the optimal lambda
#   structure.list : a list of the computed structures from the grid of lambdas
#   result.U.list : a list of loadings from the grid of lambdas



estimation.fun <- function (model,  lambda,  set.training = seq(1:(floor(0.5 * n))), 
                            r.list = NULL, ind.perm = NULL, light.ver = FALSE) {
  
  
  
  n = model$n
  K = model$K
  if (is.null(r.list)) {r.list = model$r.list}
  p = model$p
  p.list = get.list(p,K)
  
  
  if (is.null(ind.perm)) {ind.perm = ind.perm.fun(K)}
  
  
  result.seq = list()
  result.str.seq = list()
  result.lambda.seq = list()
  
  
  # Split the data matrix into training and test sets.
  
  set.test = seq(1:n)[-set.training]
  
  X.training = list()
  X.test = list()
  for(i in 1:K){
    X.training[[i]] = model$X[[i]][,set.training]
    X.test[[i]] = model$X[[i]][,set.test]
  }
  
  
  
  # list of scores from training sets
  
  V.list.untruncated = list()
  V.list = list()
  
  for (i in 1:K ) {
    V.list.untruncated[[i]] = svd(X.training[[i]])$v
    rv = min(nrow(V.list.untruncated[[i]]),ncol(V.list.untruncated[[i]]),r.list[i])
    if (rv > 0) {
      V.list[[i]] = V.list.untruncated[[i]][,1:rv  ]
    } else {
      V.list[[i]] = V.list.untruncated[[i]][,0]
    }
  }
  
  
  # Analysis on the training sets
  # fit loadings for the training sets
  
  
  result.seq = list()
  result.str.seq = list()
  result.U = list()
  

  
  r.test = rep(0,length(lambda))
  r.pen.test = rep(0,length(lambda))
  
  for (i in 1:length(lambda)) {
    
    result.seq[[i]] <- main.fun(V.list, n, K, r.list, lambda = lambda[i], ind.perm)  
    result.str.seq[[i]] <- get.structure(result.seq[[i]], K)
    result.U[[i]] <- U.estimate(K, p, result.str = result.str.seq[[i]], 
                                subspace = result.seq[[i]]$subspace, X.training)
    
    

    
    X.test.con <- do.call(rbind, X.test)
    
    svd.test <- svd(t(X.test.con) %*% result.U[[i]]$U)
    V.hat.test <- svd.test$u %*% t(svd.test$v)
    
    
    # V.svd = svd(result.U[[i]]$V)
    # V.C = V.svd$v; V.D = diag(V.svd$d^2, length(V.svd$d))
    # svd.test <- svd(t(X.test.con) %*% result.U[[i]]$U %*% t(V.C) %*% V.D  )
    # V.hat.test <- svd.test$u %*% t(svd.test$v) %*% V.D %*% V.C
    

    
    r.test[i] = 0
    ER = rep(0,K)

    for(j in 1:K) {
      ER[j] = norm( X.test[[j]] - result.U[[i]]$U[p.list[[j]], ] %*% t(V.hat.test),"F")^2 / norm(X.test[[j]], "F")^2
      r.test[i] = r.test[i] + ER[j]
    }
    
    
    
    V.hat = do.call(cbind, result.seq[[i]]$subspace)
    pen = 0
    r.pen.test[i] <- r.test[i] 
    
  }
  
  # optimal lambda and its corresponding structure 
  
  opt.lam = which.min(r.test)
  structure.hat = result.str.seq[[opt.lam]]$structure
  
  
  # scores from the whole datasets
  
  V.list.untruncated = list()
  V.list = list()
  
  for (i in 1:K ) {
    V.list.untruncated[[i]] = svd(model$X[[i]])$v
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
    result.whole[[i]] <- main.fun(V.list, n, K, r.list, lambda = lambda[i], ind.perm)  
    result.str.whole[[i]] <- get.structure(result.whole[[i]], K)
    result.U.whole[[i]] <- U.estimate(K, p, result.str = result.str.whole[[i]], 
                                      subspace = result.whole[[i]]$subspace, model$X)
    result.diff.whole[i] = distance.struct(structure.hat,result.str.whole[[i]]$structure)
  }
  
  # find a structure from the whole datasets that is closest to the structure.hat
  
  lam.hat = which.min(result.diff.whole)
  
  
  if (light.ver == TRUE) {
    return( list( structure = result.str.whole[[lam.hat]]$structure,
                  r.test.list = r.test,
                  r.pen.test.list = r.pen.test,
                  result.U = result.U.whole[[lam.hat]],
                  r.list = r.list,
                  ind.perm = ind.perm)
    )  
  } else {
    return( list( structure = result.str.whole[[lam.hat]]$structure,
                  r.test = r.test[opt.lam], r.pen.test = r.pen.test[opt.lam],
                  result.U = result.U.whole[[lam.hat]],
                  structure.list = result.str.whole,
                  diff.list = result.diff.whole,
                  r.test.list = r.test,
                  r.pen.test.list = r.pen.test,
                  result.U.list = result.U.whole,
                  r.list = r.list,
                  ind.perm = ind.perm)
    )
  }
  
}



