
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

source("./Other_methods/BIDIFAC.plus/bidifac.plus.data.prep.R")
source("./Other_methods/BIDIFAC.plus/bidifac.plus.parse.R")


model.no = 1
n = 200
# p = c(200,200,200)
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




###

fun <- function(model, trueU, lambda){
  
  result.fun = list()
  
  
  
  
  model.BIDIFAC = bidifac.plus.data.prep(model)
  
  BIDIFAC.result = bidifac.plus.parse(model.BIDIFAC)
  
  
  result.fun = list(BIDIFAC.result = BIDIFAC.result)
  
  
  
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
  cat(system.time(result.oneSNR <- foreach(i=1:nrep,
                                           .packages=c("pracma","irlba","RSpectra","dplyr","gtools")) %dopar% {
                                             fun(model = model.set[[j]][[i]], trueU = model.loading, lambda = lambda)
                                           }),'\n')
  result.total[[j]] = result.oneSNR
}




stopCluster(cl)




###

# save.name1 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".Rdata")
# 
# save(model.set, file = save.name1)
# save(result.total, file = save.name2)



######################################################################
# Draw Table 1 and Table 2
# => See "./Sim_codes/Sim_comparisons_table_(Sec5.2).R"
######################################################################



######################################################################
# Draw Table 1 and Table 2
######################################################################



library(xtable)
source("./Source/general_functions.R")
source("./Source/parser.R")
###



ncol = 1       #number of methods
SNR.list = c(10, 5, 2)
Nrep = 100 

mm = 4     #number of measures



loading.grass.angle.mean.mat = matrix(0, length(SNR.list), ncol)
loading.grass.angle.sd.mat = matrix(0, length(SNR.list), ncol)


score.grass.angle.mean.mat = matrix(0, length(SNR.list), ncol)
score.grass.angle.sd.mat = matrix(0, length(SNR.list), ncol)


str.mean.mat = matrix(0, length(SNR.list), ncol)
str.sd.mat = matrix(0, length(SNR.list), ncol)

true.str.mean.mat = matrix(0, length(SNR.list), ncol)
true.str.sd.mat = matrix(0, length(SNR.list), ncol)


for(i in 1:length(SNR.list)) {  
  
  
  
  BIDIFAC.var = matrix(0,Nrep,mm)
  
  
  for(j in 1:Nrep) {
    
    K = model.set[[i]][[j]]$K
    trueU =  model.set[[i]][[j]]$U 
    trueV =  model.set[[i]][[j]]$V.con 
    p.list = model.set[[i]][[j]]$p.list
    
    trueStructure = model.set[[i]][[j]]$true.struct
    
    X = do.call(rbind, model.set[[i]][[j]]$X)
    
    
    Z <- list()
    for(k in 1:K) {
      Z[[k]] = tcrossprod(model.set[[i]][[j]]$U[[k]], model.set[[i]][[j]]$V[[k]])     
    }
    
    
    rank_total = sum(model.set[[i]][[j]]$r.list)
    
    true.r = model.set[[i]][[j]]$r
    
    true.r.cumsum = cumsum(true.r)
    true.r.list = list()
    true.r.list[[1]] = 1:true.r.cumsum[1]
    if (length(true.r) > 1) {
      for (k in 2:length(true.r)) {
        true.r.list[[k]] = (true.r.cumsum[k-1] + 1) : (true.r.cumsum[k])  
      }
    }
    
    true.r.list.part = list()
    for (k in 1:K) {
      b = list(); b.num = 0; b.rank = 0
      for(kk in 1:(length(true.r.list))){
        if (trueStructure[k, true.r.list[[kk]][1]  ] == 1) { 
          b.num = b.num + 1
          b[[b.num]] = seq(b.rank + 1, b.rank + length(true.r.list[[kk]]))  
          b.rank = b.rank + length(true.r.list[[kk]])
        }
      }  
      true.r.list.part[[k]] = b
    }
    
    
    ######################
    
    
    
    BIDIFAC.loading = result.total[[i]][[j]]$BIDIFAC.result$BIDIFAC.loading
    BIDIFAC.structure = sort.columns(result.total[[i]][[j]]$BIDIFAC.result$BIDIFAC.struct, K,
                                     transpose = TRUE)$structure
    BIDIFAC.score = result.total[[i]][[j]]$BIDIFAC.result$BIDIFAC.score
    
    
    
    
    
    
    
    
    
    if (ncol(BIDIFAC.loading) > 0 && ncol(BIDIFAC.score) > 0) {
      
      
      
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(length(true.r.list))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,true.r.list[[kkk]] ])$u, svd(BIDIFAC.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,true.r.list[[kkk]]])$u, svd(BIDIFAC.score)$u)))
      } 
      
      
      
      
      
      
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        
        if (length(which(BIDIFAC.structure[k,] != 0)) > 0) {
          
          for(kkk in 1:(length(true.r.list.part[[k]]))) {
            
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ] )$u, 
                                    svd(BIDIFAC.loading[p.list[[k]],  which(BIDIFAC.structure[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ]  )$u, 
                                            svd(BIDIFAC.loading[p.list[[k]], which(BIDIFAC.structure[k,] != 0)])$u)))
          }
        }
      }
      
      
      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]-BIDIFAC.loading[p.list.c,] %*% t(BIDIFAC.score))^2) / sum(Z[[k]]^2) 
      }
      RSE = RSE / K
      
      
      
      BIDIFAC.var[j, 4] = mean(angle.loading) 
      BIDIFAC.var[j, 3] = mean(angle.score)  
      BIDIFAC.var[j, 2] = RSE
      BIDIFAC.var[j, 1] = identical(trueStructure, BIDIFAC.structure )
    } 
    else {
      BIDIFAC.var[j,] = rep(NA, mm)
    }
    
    
    
  }
  
  
  
  
  
  
  
  
  
  loading.grass.angle.mean.mat[i,1] = mean(BIDIFAC.var[,4])
  loading.grass.angle.sd.mat[i,1] = sd(BIDIFAC.var[,4])
  
  score.grass.angle.mean.mat[i,1] = mean(BIDIFAC.var[,3])
  score.grass.angle.sd.mat[i,1] = sd(BIDIFAC.var[,3])
  
  str.mean.mat[i,1] = mean(BIDIFAC.var[,2])
  str.sd.mat[i,1] = sd(BIDIFAC.var[,2])
  
  true.str.mean.mat[i,1] = mean(BIDIFAC.var[,1])
  true.str.sd.mat[i,1] = sd(BIDIFAC.var[,1])
  
  
  
  
  
  
}


names = c("BIDIFAC")

total.mean.mat = matrix(0, mm * length(SNR.list), ncol)
total.sd.mat = matrix(0, mm * length(SNR.list), ncol)

for(i in 1:length(SNR.list)) {
  total.mean.mat[ mm * (i-1) + 3 , ] = loading.grass.angle.mean.mat[i,]
  total.mean.mat[ mm * (i-1) + 4 , ] = score.grass.angle.mean.mat[i,]
  total.mean.mat[ mm * (i-1) + 2 , ] = str.mean.mat[i,]
  total.mean.mat[ mm * (i-1) + 1 , ] = true.str.mean.mat[i,] * 100
  
  total.sd.mat[ mm * (i-1) + 3 , ] = loading.grass.angle.sd.mat[i,]
  total.sd.mat[ mm * (i-1) + 4 , ] = score.grass.angle.sd.mat[i,]
  total.sd.mat[ mm * (i-1) + 2 , ] = str.sd.mat[i,]
  total.sd.mat[ mm * (i-1) + 1 , ] = true.str.sd.mat[i,] * 100
}



msd <- paste(round(total.mean.mat,2)," (",round(total.sd.mat,2),")",sep="")

tab.msd <- matrix(msd, mm * length(SNR.list), ncol)
colnames(tab.msd) = names

# Measure = rep(c("FFFF", "EEEE", "BBBB", "DDDD"), length(SNR.list))

Measure = rep(c("Accuracy  ", "RSE", "Loading  ", "Score     "), length(SNR.list))


SNR = c()
for (i in 1:length(SNR.list)) {
  SNR = c(SNR, as.character(SNR.list[i]))
  SNR = c(SNR, rep(" ", mm - 1))
}
# SNR = rep(as.character(SNR.list), each = mm)


tab.msd = cbind(SNR, Measure, tab.msd)


print(xtable(tab.msd),  include.rownames=FALSE)

