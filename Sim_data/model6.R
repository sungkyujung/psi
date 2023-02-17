

source("./Source/general_functions.R")

model.generate <- function(n = 200, K = 3, p = c(100,100,100), weight = c(1,1,1), 
                           U = NULL, V.var = c(1.8, 0.8, 1.7, 0.7, 1.6, 0.6, 1.5, 0.5, 1.4, 0.4, 1.3, 0.3, 1.2, 0.2),
                           SNR = 0,  rank = c(2,2,2,2,2,2,2))
{
  
  
  
  p.total <- sum(p)
  p.list = get.list(p,K)
  
  
  
  # List of index for r
  # For column indices, the order is  r123, r12, r13, r23, r1, r2, r3
  
  r.list = c(rank[1] + rank[2] + rank[3] + rank[5], 
             rank[1] + rank[2] + rank[4] + rank[6], 
             rank[1] + rank[3] + rank[4] + rank[7])
  r.index.list = get.list(rank, length(rank))
  
  r.structure = list()
  # r.structure[[1]] = c(1,2,3,4,5,6,9,10)
  # r.structure[[2]] = c(1,2,3,4,7,8,11,12)
  # r.structure[[3]] = c(1,2,5,6,7,8,13,14)
  
  r.structure[[1]] = c(r.index.list[[1]],r.index.list[[2]],r.index.list[[3]],r.index.list[[5]])
  r.structure[[2]] = c(r.index.list[[1]],r.index.list[[2]],r.index.list[[4]],r.index.list[[6]])
  r.structure[[3]] = c(r.index.list[[1]],r.index.list[[3]],r.index.list[[4]],r.index.list[[7]])
  
  # Frame for X
  X <- matrix(0, p.total, n)
  
  
  # Generate loadings
  
  r.U = list()
  r.U[[1]] = list(r.index.list[[1]],r.index.list[[2]],r.index.list[[3]],r.index.list[[5]])
  r.U[[2]] = list(r.index.list[[1]],r.index.list[[2]],r.index.list[[4]],r.index.list[[6]])
  r.U[[3]] = list(r.index.list[[1]],r.index.list[[3]],r.index.list[[4]],r.index.list[[7]])
  
  
  
  
  if(is.null(U)) {
    U <- list()
    for(k in 1:K){
      cnt = 0
      U[[k]] = matrix(0,p[k],length(r.structure[[k]]))
      for(i in 1:length(r.U[[k]])){
        temp <- matrix(rnorm(p[k] *length(r.U[[k]][[i]]),0,sqrt(V.var[r.U[[k]][[i]] ]) ), p[k], length(r.U[[k]][[i]]))
        U[[k]][,cnt + 1:length(r.U[[k]][[i]])] = scale( temp , scale = F)
        cnt = cnt + length(r.U[[k]][[i]])
      }
    }
    
  }
  
  
  # Generate score components
  # For column indices, the order is r12, r13, r23
  
  V.comp <- list()
  for (i in 1:length(rank) ){
    temp = matrix(0, n, rank[i])
    for (j in 1:length(r.index.list[[i]])) {
      temp[,j] <- matrix(rnorm(n,0,1), n, 1)
    }
    V.comp[[i]] = qr.Q(qr(scale( temp , scale = F)))
  }
  
  
  
  # Put the score components together into each scores matrix V.i
  # For column indices, the order is   r123, r12, r13, r23, r1, r2, r3
  V <- list()
  
  V[[1]] <- cbind(V.comp[[1]], V.comp[[2]], V.comp[[3]], V.comp[[5]])
  V[[2]] <- cbind(V.comp[[1]], V.comp[[2]], V.comp[[4]], V.comp[[6]])  
  V[[3]] <- cbind(V.comp[[1]], V.comp[[3]], V.comp[[4]], V.comp[[7]])
  
  
  # Summing up all the loadings and scores into the dataset X
  X <- list()
  for(i in 1:K) {
    X[[i]] = tcrossprod(U[[i]], V[[i]]) * weight[i]
  }
  
  
  # Add noise
  if (SNR != 0){
    SNR = SNR * (n / 200) 
    for(i in 1:K) {
      sigma <- sqrt(  1 / SNR  )
      X[[i]] = X[[i]] +   matrix(rnorm(n * p[i], sd = sigma),p[i],n)
    }
  }
  
  
  # Return the model
  
  true.struct = matrix(c(1,1,1, 
                         1,1,0, 
                         1,0,1, 
                         1,0,0,
                         0,1,1,
                         0,1,0,
                         0,0,1),3,7)
  true.struct = true.struct[,c(rep(1,rank[1]),
                               rep(2,rank[2]),
                               rep(3,rank[3]),
                               rep(4,rank[4]),
                               rep(5,rank[5]),
                               rep(6,rank[6]),
                               rep(7,rank[7]))]
  
  
  # Concatenated version of U, V
  
  U.con = matrix(0, sum(p), sum(rank))
  for (k in 1:K) {
    U.con[p.list[[k]], r.structure[[k]] ] = U[[k]]
  }
  V.con = do.call(cbind, V.comp)
  
  
  
  variable_names = c("n", "K","p","weght","SNR",
                     "X","V","U","V.comp","V.con", "U.con",
                     "r.list","r.structure","r","p.list","true.struct")
  model = list(n, K, p, weight, SNR, 
               X, V, U, V.comp, V.con, U.con, 
               r.list, r.structure, rank, p.list, true.struct)
  model = setNames(model, variable_names)
  
  return(model)
}

