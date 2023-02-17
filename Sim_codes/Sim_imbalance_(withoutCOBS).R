
###############################################################################
# Source code for Section 5.3
# Results on Imbalanced Signal Strength between Joint and Individual Components
###############################################################################




source("./Source/general_functions.R")
source("./Source/parser.R")

source("./Other_methods/SLIDE/SLIDE.parse.R")
source("./Other_methods/SLIDE/Data_preparations.R")

source("./Other_methods/COBS/COBS.parse.R")

source("./Other_methods/AJIVE/AJIVE.parse.R")

source("./Other_methods/JIVE/JIVE.parse.R")
source("./Other_methods/JIVE/JIVE.data.R")



model.no = 4
n = 200
SNR.list = c(0, 10)
lambda = seq(0,90,1) * pi / 180

model.name = paste0("./Sim_data_New/model",model.no,".R")
source(model.name)



### (Load the saved model)

# save.name1 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/ind.data.Rdata")
# save.name1 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/joint.data.Rdata")
# load(file = save.name1)



### model generation 
### generate 110 data sets and use 100 of them, because of error-handling due to COBS

Nrep = 110
nrep = 100



### Joint-emphasized imbalanced model (for Table 3)
# VV = c(seq(15, 5.5, by = -1),
#        seq(0.15, 0.069, by = -0.009),
#        seq(0.147, 0.066, by = -0.009),
#        seq(0.144, 0.063, by = -0.009))

### Individual-emphasized imbalanced model (for Table 4)
VV = c(seq(1.5, 0.60, by = -0.1),
       seq(15, 6.9, by = -0.9),
       seq(14.7, 6.6, by = -0.9),
       seq(14.4, 6.3, by = -0.9))



set.seed(1001)

model.loading = model.generate(n = n, U = NULL, SNR = 0,
                               rank = c(10,10,10,10),
                               V.var = VV)$U


set.seed(1001)

model.set = list()
for (j in 1:length(SNR.list)) {
  model.list.SNR = list()
  for (i in 1:Nrep) {
    model.list.SNR[[i]] = model.generate(SNR = SNR.list[j], U = model.loading, n = n,
                                         rank = c(10,10,10,10),
                                         V.var = VV)
  }
  model.set[[j]] = model.list.SNR
}





### Function for Calling Methods

fun <- function(model, trueU, lambda){
  
  
  model.SLIDE = list()
  model.COBS = list()
  
  result.fun = list()
  
  ########## Data preparation for each method
  
  ## SLIDE
  model.SLIDE =   data.prep.SLIDE( model)
  
  
  ## COBS
  Y.bk.idx = c()
  for (k in 1:model$K) {
    Y.bk.idx = c(Y.bk.idx, rep(k, model$p[k]))
  }
  model.COBS =   list( Y = t(do.call(rbind,model$X)),
                       Y.bk.idx = Y.bk.idx)
  
  ## AJIVE
  XX = list()
  for(kk in 1:model$K) {
    XX[[kk]] = t(model$X[[kk]])
  }
  
  ## JIVE
  
  model.JIVE = data.prep.JIVE( model,  r.list = NULL)
  ###
  
  
  
  ########## Call parser functions
  
  
  
  PSI.result = PSI.parse(model, lambda, r.list = model$r.list)
  
  SLIDE.result = SLIDE.parse(model = model.SLIDE)
  
  # COBS.result = COBS.parse(model.COBS, 
  #                          n.comp = sum(model$r), 
  #                          alpha.v = 0.5, 
  #                          method = "BIC-Low", 
  #                          M.lambda.v = 40, 
  #                          orth = F, 
  #                          p.list = model$p.list, 
  #                          K = model$K)
  COBS.result = list(COBS.loading = matrix(rnorm(sum(model$p),0,1),sum(model$p),1), 
                     COBS.score = matrix(rnorm(model$n,0,1),model$n,1), 
                     COBS.struct = matrix(0, model$K, 1))
  
  AJIVE.result = AJIVE.parse(XX, model$r.list,
                             n = model$n,
                             p.list = model$p.list, 
                             p = model$p,
                             K = model$K)
  
  JIVE.result = JIVE.parse(model, model.JIVE, "perm")
  
  result.fun = list(PSI.result = PSI.result,
                    SLIDE.result = SLIDE.result,
                    COBS.result = COBS.result,
                    AJIVE.result = AJIVE.result,
                    JIVE.result = JIVE.result)
  
  
  
  return(result.fun)
}


###



library(doParallel)
library(doRNG)
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)


set.seed(1001)


nrep = 100

result.total = list()
result.oneSNR = list()

for (j in 1:length(SNR.list)){
  cat(system.time(result.oneSNR <- foreach(i=1:Nrep, .errorhandling='remove',
                                           .packages=c("pracma","irlba","RSpectra","dplyr","gtools")) %dopar% {
                                             fun(model = model.set[[j]][[i]], trueU = model.loading, lambda = lambda)
                                           }),'\n')
  result.total[[j]] = result.oneSNR[1:nrep]
}





stopCluster(cl)



### Save the datasets and the results.

# save.name1 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/joint.data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/joint.Rdata")

# save.name1 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/ind.data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_imbalance_data_(Sec5.3)/ind.Rdata")




# save(model.set, file = save.name1)
# save(result.total, file = save.name2)



