
library(pracma)
library(irlba)
library(RSpectra)


source("./Source/Algorithm.R")
source("./Source/general_functions.R")
source("./Source/Estimation_fun.R")
source("./Source/Reconstruction.R")
source("./Source/getnfac.R")


PSI.parse <- function(model, lambda, r.list = NULL, center = TRUE) {
  
  
  n = model$n
  
  
  if (center) {
    for (k in 1:model$K) {
      model$X[[k]] = t(scale(t(model$X[[k]]), scale = FALSE))
    }  
  }
  
  if (is.null(r.list)) {
    r.list = c()

    for (k in 1:model$K) {
      nc = ncol(t(model$X[[k]]))
      testPCA = getnfac(t(model$X[[k]]),kmax = floor(nc/2), criteria = 'IC3')
      
      r.list = c(r.list, testPCA$ic)

    }
  } 
  
  

  
  PSI.result = estimation.fun(model,  lambda, set.training =  sample(n, ceiling(n/2)), r.list = r.list)

  return(list(PSI.loading = PSI.result$result.U$U, 
              PSI.struct = PSI.result$structure,
              PSI.score = PSI.result$result.U$V)) 
}
