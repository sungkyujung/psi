
###############################################################################
# Source code for Section 5.2
# Results on Comparative Study
###############################################################################




source("./Source/general_functions.R")
source("./Source/parser.R")

source("./Other_methods/SLIDE/SLIDE.parse.R")
source("./Other_methods/SLIDE/Data_preparations.R")

source("./Other_methods/COBS/COBS.parse.R")

source("./Other_methods/AJIVE/AJIVE.parse.R")

source("./Other_methods/JIVE/JIVE.parse.R")
source("./Other_methods/JIVE/JIVE.data.R")



### *** SWITCH THE MODEL SETTINGS HERE ***


model.no = 6
n = 200
p = c(100,100,100)
SNR.list = c(15,10,5)
lambda = seq(0,90,1) * pi / 180

model.name = paste0("./Sim_data_New/model",model.no,".R")
source(model.name)



### (Load the saved model)

# save.name1 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".data.Rdata")
# load(file = save.name1)



### Example model generation

set.seed(1001)

model.loading = model.generate(n = n, p = p, U = NULL, SNR = 0)$U



set.seed(1001)

nrep = 100

model.set = list()
for (j in 1:length(SNR.list)) {
  model.list.SNR = list()
  for (i in 1:nrep) {
    model.list.SNR[[i]] = model.generate(SNR = SNR.list[j], U = model.loading, n = n, p = p)
  }
  model.set[[j]] = model.list.SNR
}




### A parser function for comparative methods

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
  
  # model.JIVE = data.prep.JIVE( model, sep = TRUE, r.list = model$r.list)
  model.JIVE = data.prep.JIVE( model,  r.list = NULL)
  
  ########## Call parser functions
  
  
  PSI.result = PSI.parse(model, lambda, r.list = model$r.list)
  
  # SLIDE.result = SLIDE.parse(model = model.SLIDE, rank_total = sum(model$r))
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


### Execute the parser function using "doParallel" package


library(doParallel)
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)




set.seed(1001)
nrep = 100


result.total = list()
result.oneSNR = list()


for (j in 1:length(SNR.list)){
  cat(system.time(result.oneSNR <- foreach(i=1:nrep,
                                           .packages=c("pracma","irlba","RSpectra","dplyr","gtools")) %dopar% {
                                             fun(model = model.set[[j]][[i]], trueU = model.loading, lambda = lambda)
                                           }),'\n')
  result.total[[j]] = result.oneSNR
}




stopCluster(cl)




### Save the generic datasets and the results

# save.name1 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".result.Rdata")

# save(model.set, file = save.name1)
# save(result.total, file = save.name2)






######################################################################
# Draw Table 1 and Table 2
# => See "./Sim_codes/Sim_comparisons_table_(Sec5.2).R"
######################################################################




