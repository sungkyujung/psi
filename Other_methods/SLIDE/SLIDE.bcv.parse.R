

source("./Other_methods//SLIDE/SLIDE.R")


SLIDE.bcv.parse <- function(model, nl = 30, n_fold = 3, 
                        p_fold = 3, k_max = 1000, eps = 1e-6, reduced = F, 
                        rank_total = NULL, center = T){
  
  
  
  
  out_s <- standardizeX(model$X, model$pvec, center = TRUE)
  
  trueU = model$trueV
  trueV = model$trueU
  
  
  X <- out_s$X
  svec <- out_s$svec
  pvec = model$pvec
  pcum <- cumsum(pvec)
  d <- length(pvec)
  
  
  
  # Solve for each lambda
  out <- solve_optim1_seq(X = X, pvec = pvec, n_lambda = 30, reduced = FALSE, 
                          rank_total = rank_total,
                          lambda_max = max(svec), lambda_min = 0.05)
  lambda_all <- out$lambda
  
  
  
  # Calculate ranks for each lambda
  if (length(pvec) == 2){
    ranks <- getAllranks_d2(out, pvec)
  }else{
    ranks <- getAllranks_d3(out, pvec)
  }
  
  # Get all the structures
  out_struct <- get_structure_v2(out, pvec)
  
  if (sum(out_struct$s_used)==0){
    stop("Something went terribly wrong")
  }
  
  errors = rep(0, length(out_struct$Slist))
  for (l in 1:length(out_struct$Slist)){
    distUV = 0
    outV <- est_givenranks_v4(X = X, pvec = pvec, pattern = out_struct$Slist[[l]], k_max = k_max, eps = eps)
    Vadj <- outV$V
    for (i in 1:d){
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      Vadj[index,] <- outV$V[index,]*sqrt(out_s$norms[i])
      # Calculate Frobenius norm squared of the signal part
      signal = sum(tcrossprod(trueU[[i]],trueV[[i]])^2)
      # If nonzero, calculate error relative to signal otherwise just directly
      if (signal > 0){
        distUV <- distUV + sum((tcrossprod(trueU[[i]],trueV[[i]])-tcrossprod(outV$U, Vadj[index,]))^2)/sum(tcrossprod(trueU[[i]],trueV[[i]])^2)
      }else{
        distUV <- distUV + sum((tcrossprod(outV$U, Vadj[index,]))^2)
      }
    }
    errors[l] = distUV
  }
  
  # Find the best tuning parameter using bcv from the structures above
  outbcv <- bcv_optim1_structure_centering(X, pvec = pvec, structure_list = out_struct$Slist, n_fold = n_fold, p_fold=p_fold, k_max = k_max, eps = eps, center = center)
  
  
  # out of the structures considered, which one was picked
  bcv_id_used <- which(outbcv$error_mean == min(outbcv$error_mean))
  ##### modification
  if (length(bcv_id_used) > 1) {bcv_id_used = bcv_id_used[1]}
  
  # out of all the lambdas, which one does it correspond to
  bcv_id <- which(out_struct$s_used == 1)[bcv_id_used]
  
  

  bcv.struct = t(out_struct$Slist[[bcv_id_used]])
  
  bcv.loading = est_givenranks_v4(X = X, pvec = pvec, pattern = out_struct$Slist[[bcv_id_used]], k_max = k_max, eps = eps)$V
  
  bcv.score = est_givenranks_v4(X = X, pvec = pvec, pattern = out_struct$Slist[[bcv_id_used]], k_max = k_max, eps = eps)$U
  
  
  return(list( SLIDE.bcv.loading = bcv.loading, 
               SLIDE.bcv.score = bcv.score, 
              SLIDE.bcv.struct = bcv.struct))
}