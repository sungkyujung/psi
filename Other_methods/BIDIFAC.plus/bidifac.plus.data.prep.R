


bidifac.plus.data.prep <- function(model, r.list = NULL) {
  
  n.ind = list(1:model$n)
  
  if (is.null(r.list)) {X = do.call(rbind, model$X)}
  else {
    X.temp = list()
    for(i in 1:model$K) {
      svd.X = svd(model$X[[i]])
      svd.X$d[(r.list[i] + 1):length(svd.X$d)] = 0
      X.temp[[i]] = svd.X$u %*% diag(svd.X$d) %*% t(svd.X$v)
    }
    X = do.call(rbind, X.temp)
    X.temp = NULL
  }
  
  return(list(X = X, n.ind = n.ind, p.ind = model$p.list))
}