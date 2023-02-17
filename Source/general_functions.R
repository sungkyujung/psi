

# generating a power set

power.set <- function(K){
  
  a = seq(1:K)
  res = matrix(FALSE,length(a),0)
  for (i in seq(length(a),1,-1)) {
    A = combn(a,i)
    for (j in 1:ncol(A)) {
      b = matrix(FALSE,length(a),1)
      b[A[,j] ,] = TRUE
      res = cbind(res, b)
    }
  }
  
  res
}


# index-set permutation

ind.perm.fun <- function(K) {
  cnt = 0
  facK = factorial(K)
  result = c()
  for(i in 0:(K-1)) {
    perm.size = facK / factorial(i) / factorial(K-i) 
    result = c(result,sample(1:perm.size) + cnt)
    cnt = cnt + perm.size
  }
  
  return(result)
}


# get the index list

get.list <- function(p, K) {
  p.cumsum = cumsum(p)
  p.list = list()
  p.list[[1]] = 1:p.cumsum[1]
  if (K > 1) {
    for (i in 2:K) {
      p.list[[i]] = (p.cumsum[i-1] + 1) : (p.cumsum[i])  
    }
  }
  
  return (p.list)
}



# Householder transform

Householder <- function(A, B) {
  A = as.matrix(A); B = as.matrix(B); n.b = ncol(B)
  A.house = A
  for (i in 1:n.b) {
    d = (A.house[,i] - B[,i])
    d = d / norm(d,"2")
    A.house = (diag(nrow(A)) - 2 * d %*% t(d)) %*% A.house
  }
  return (A.house)
}


# principal angles btw two subspaces

p.ang <- function(A, B) {
  
  if (ncol(A) == 0 || ncol(B) == 0) {return(NA)}
  
  return(acos(pmin(pmax(svd(t(A) %*% B)$d,-1),1)) * 180 / pi)
}


# projection of B on the subspace made from A and orthonormalize by QR decomposition

proj.QR <- function(A, B) {
  C = A %*% solve(t(A) %*% A ) %*% t(A) %*% B
  return ( qr.Q(qr(C)) )
}


# subspace average via SVD (wrt projection F-norm)


# SVDavr <- function(V.list) {
#   V = do.call(cbind,V.list)
  
SVDavr <- function(V) {
  
  
  ### due to instability of irlba, we use RSpectra when ncol(V) >= 3
  ### suppressWarnings( { svd.V = irlba(V, nu = 1, nv = 1) })
  
  svd.V <- tryCatch(
    {  if (ncol(V) >= 3 ) {
        if(require(RSpectra)) svd.V = svds(V,1,nu = 1, nv = 1)
      } else {
        if(require(irlba)) suppressWarnings( { svd.V = irlba(V, nu = 1, nv = 1) })
      }
    }, 
    error = function(ee){
      svd.V = svd(V, nu = 1, nv = 1)
    },
    finally = {
      #pass
    }
  )
  
  # suppressWarnings( { svd.V = irlba(V, nu = 1, nv = 1) })
  
  return(svd.V$u)
}


# partial jointness structure of the result

get.structure <- function(r, K) {
  
  rank.joint = r$ranks[1]
  len.clus = r$ranks[2]
  rank.total = r$ranks[3]
  
  # Extract the (partial) joint structure in the score subspace.
  X = matrix(0, rank.total, K)
  
  if (rank.joint != 0) {
    for (i in 1:rank.joint) {
      X[i, r$cluster[[i]] ] = 1
    }
  }
  if (rank.joint != len.clus) {
    j = rank.joint 
    for (i in (rank.joint + 1):len.clus ) {
      X[ (j + 1) : (j + ncol(r$subspace[[i]]) ) , r$cluster[[i]] ] = 1
      j = j + ncol(r$subspace[[i]])
    }
  }
  
  # Sort the rows of X in lexicographic order
  X = as.data.frame(X)
  lex.order = do.call(order, c(X[1:K], decreasing = TRUE))
  X.sort = t(unname(as.matrix(X[lex.order, ])))
  
  

  
  return(list(structure = X.sort, lex.order = lex.order))
}


# sort the columns of structure matrix in a lexicographic order

sort.columns <- function(X, K, transpose = FALSE) {
  
  if (transpose == TRUE) {X = t(X)}
  
  # Sort the rows of X in lexicographic order
  X = as.data.frame(X)
  lex.order = do.call(order, c(X[1:K], decreasing = TRUE))
  X.sort = t(unname(as.matrix(X[lex.order, ])))
  
  return(list(structure = X.sort, lex.order = lex.order))
}


# distance between two structure matrix

distance.struct.f <- function(X, Y) {
  X = as.matrix(X); Y = as.matrix(Y)
  
  nrow.X = nrow(X); ncol.X = ncol(X)
  nrow.Y = nrow(Y); ncol.Y = ncol(Y)
  nrow.max = max(nrow.X, nrow.Y); ncol.max = max(ncol.X, ncol.Y)
  
  X.new = matrix(0, nrow.max, ncol.max)
  X.new[1:nrow.X, 1:ncol.X] = X
  
  Y.new = matrix(0, nrow.max, ncol.max)
  Y.new[1:nrow.Y, 1:ncol.Y] = Y
  
  Z = X.new - Y.new
  Z[which(Z != 0 )] = 1
  
  return(norm(Z,"F"))
}




# distance between two structure matrix

distance.struct <- function(X, Y) {
  X = as.matrix(X); Y = as.matrix(Y)
  
  nrow.X = nrow(X); ncol.X = ncol(X)
  nrow.Y = nrow(Y); ncol.Y = ncol(Y)
  nrow.max = max(nrow.X, nrow.Y); ncol.max = max(ncol.X, ncol.Y)
  
  out.X = rep(0,ncol.X) ; out.Y = rep(0,ncol.Y)
  for (i in 1:ncol.X) {
    for (j in 1:ncol.Y) {
      if (out.X[i] > 0 || out.Y[j] > 0) {next}
      if (identical(X[,i], Y[,j])) {out.X[i] = 1; out.Y[j] = 1; break}
    }
  }
  
  X.res = as.matrix(X[,which(out.X == 0)]); Y.res = as.matrix(Y[,which(out.Y == 0)])
  
  
  if (ncol(X.res) > 0 && ncol(Y.res) > 0) {
    dis.X = rep(nrow.X^2, ncol(X.res))
    dis.Y = rep(nrow.Y^2, ncol(Y.res))
    
    for (i in 1:ncol(X.res)) {
      for (j in 1:ncol(Y.res)) {
        dis.X[i] = min(dis.X[i], distance.struct.f(X.res[,i],Y.res[,j])^2 )
      }
    }
    
    for (j in 1:ncol(Y.res)) {
      for (i in 1:ncol(X.res)) {
        dis.Y[j] = min(dis.Y[j], distance.struct.f(X.res[,i],Y.res[,j])^2 )
      }
    }
    
    dis.X.sum = sum(dis.X)
    dis.Y.sum = sum(dis.Y)
  } else {
    dis.X.sum = sum(X.res); dis.Y.sum = sum(Y.res)
  } 
  
  
  dis.sum = dis.X.sum + dis.Y.sum
  
  return(dis.sum)
}






# Initial Values for U


U.estimate <- function(K, p, result.str, subspace, X.training) {
  
  str = result.str$structure
  lex.order = result.str$lex.order
  
  p.list = get.list(p, K)
  
  V0 = as.matrix(do.call(cbind,subspace)[ ,result.str$lex.order]) 


  U = matrix(0, nrow = sum(p), ncol = ncol(V0))
  
  for (i in 1:K) {
    ls = which(str[i,] == 1 )
    if (length(ls) > 0) {
      X_ = X.training[[i]]
      V_ = V0[,ls] 
      U[ p.list[[i]] ,ls ] = X_ %*% V_ %*% solve(t(V_) %*% V_)
    } 
  }
  
  return(list(U = U, V = V0))
}


str.trunc <- function(structure){
  K = nrow(structure)
  cs = colSums(structure)
  indv.indx = which(cs == 1)
  if(length(indv.indx) == ncol(structure)){
    joint.indx = numeric(0)
    return(list(structure = matrix(0,K,0), index = joint.indx))  
  } else if (length(indv.indx) == 0){
    return(list(structure = structure, index = seq(1:ncol(structure))))
  } else {
    joint.indx = seq(1:ncol(structure))[-indv.indx]
    return(list(structure = structure[,joint.indx], index = joint.indx))
  }
  
}





## Other functions

normalize.vector <- function(x) {x / norm(x, "2")}
normalize.matrix <- function(x) {x / norm(x, "2")}

