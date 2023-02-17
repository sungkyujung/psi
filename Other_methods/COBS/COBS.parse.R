source("./Other_methods/COBS/get.initial.R")
source("./Other_methods/COBS/update.b.R")
source("./Other_methods/COBS/update.v.R")
source("./Other_methods/COBS/summary.cobs.R")
source("./Other_methods/COBS/cobs.tune.R")
source("./Other_methods/COBS/bisect.dec.R")
source("./Other_methods/COBS/Sigmas_Full.R")
source("./Other_methods/COBS/Sigma.est.R")
source("./Other_methods/COBS/COBS function in r.R")
source("./Other_methods/COBS/COBS.R")

library(dplyr)


COBS.parse <- function(model.COBS, n.comp, alpha.v, method, M.lambda.v, orth, p.list, K){

  
  
  COBS.tune.result = cobs.tune(dataY = model.COBS$Y, Y.bk.idx = model.COBS$Y.bk.idx,
                          n.comp = n.comp,
                          alpha.v = alpha.v,
                          method = method,
                          M.lambda.v = M.lambda.v,
                          orth=orth)
  
  COBS.loading = COBS.tune.result$V
  COBS.score = COBS.tune.result$F.est
  
  COBS.struct = matrix(0,K,n.comp)
  for(i in 1:n.comp) {
    for(j in 1:K) {
      if (prod(COBS.loading[p.list[[j]], i] == 0) != 1) {
        COBS.struct[j,i] = 1
      }
    }
  }
  
  return(list(COBS.loading = COBS.loading, COBS.score = COBS.score, COBS.struct = COBS.struct))
}