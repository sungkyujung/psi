

source("./Source/general_functions.R")

model.generate <- function(n = 200, K = 3, p = c(100,100,100), weight = c(1,1,1), 
                           U = NULL, V.var = c(1.0, 0.9),
                           SNR = 0,  rank = c(2))
{
  
  
  
  p.total <- sum(p)
  p.list = get.list(p,K)
  
  
  
  # List of index for r
  # For column indices, the order is  r1, r1, r1
  
  r.list = c(rank[1] ,rank[1] ,rank[1] )
  r.index.list = get.list(rank, length(rank))
  
  r.structure = list()
  # r.structure[[1]] = c(1,2)
  # r.structure[[2]] = c(1,2)
  # r.structure[[3]] = c(1,2)
  
  r.structure[[1]] = c(r.index.list[[1]])
  r.structure[[2]] = c(r.index.list[[1]])
  r.structure[[3]] = c(r.index.list[[1]])
  
  # Frame for X
  X <- matrix(0, p.total, n)
  
  
  # Generate loadings
  
  r.U = list()
  r.U[[1]] = list(r.index.list[[1]])
  r.U[[2]] = list(r.index.list[[1]])
  r.U[[3]] = list(r.index.list[[1]])
  
  
  
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
  
  V.comp <- list()
  for (i in 1:length(rank) ){
    temp = matrix(0, n, rank[i])
    for (j in 1:length(r.index.list[[i]])) {
      temp[,j] <- matrix(rnorm(n,0,1), n, 1)
    }
    V.comp[[i]] = qr.Q(qr(scale( temp , scale = F)))
  }
  
  
  
  
  
  # Put the score components together into each scores matrix V.i
  # For column indices, the order is  r1, r1, r1
  V <- list()
  
  V[[1]] <- cbind(V.comp[[1]])
  V[[2]] <- cbind(V.comp[[1]])  
  V[[3]] <- cbind(V.comp[[1]])

  
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
  
  true.struct = matrix(c(1,1,1),3,1)
  true.struct = true.struct[,c(rep(1,rank[1]))]
  
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

