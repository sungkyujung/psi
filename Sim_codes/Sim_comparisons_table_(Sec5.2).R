######################################################################
# Draw Table 1 and Table 2
######################################################################


######################################################################
# Needs both "model.set" and "result.total" variables
# from Sim_comparisons_(Sec5.2).R
# (or "m#.data.Rdata" and "m#.result.Rdata")
######################################################################




library(xtable)
source("./Source/general_functions.R")
source("./Source/parser.R")
source("./Other_methods/SLIDE/Data_preparations.R")



ncol = 5       #number of methods
mm = 4     #number of measures

SNR.list = c(15, 10, 5)
Nrep = 100 



### (Load the saved model and results)

# save.name1 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".data.Rdata")
# save.name2 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".result.Rdata")
# 
# load(file = save.name1)
# load(file = save.name2)




loading.grass.angle.mean.mat = matrix(0, length(SNR.list), ncol)
loading.grass.angle.sd.mat = matrix(0, length(SNR.list), ncol)


score.grass.angle.mean.mat = matrix(0, length(SNR.list), ncol)
score.grass.angle.sd.mat = matrix(0, length(SNR.list), ncol)


str.mean.mat = matrix(0, length(SNR.list), ncol)
str.sd.mat = matrix(0, length(SNR.list), ncol)

true.str.mean.mat = matrix(0, length(SNR.list), ncol)
true.str.sd.mat = matrix(0, length(SNR.list), ncol)




### Outer loop (SNR.list)

for(i in 1:length(SNR.list)) {  
  
  
  PSI.var = matrix(0,Nrep,mm)
  SLIDE.best.var = matrix(0,Nrep,mm)
  SLIDE.bcv.var = matrix(0,Nrep,mm)
  COBS.var = matrix(0,Nrep,mm)
  AJIVE.var = matrix(0,Nrep,mm)
  JIVE.var = matrix(0,Nrep,mm)

  
  ### Inner loop (Nrep)
    
  for(j in 1:Nrep) {
    
    #########################################
    ### TRUE MODEL Specifications
    #########################################
    
    K = model.set[[i]][[j]]$K
    trueU =  model.set[[i]][[j]]$U 
    trueV =  model.set[[i]][[j]]$V.con 
    p.list = model.set[[i]][[j]]$p.list
    
    trueStructure = model.set[[i]][[j]]$true.struct
    
    X = do.call(rbind, model.set[[i]][[j]]$X)
    
    
    ## True signal
    Z <- list()
    for(k in 1:K) {
      Z[[k]] = tcrossprod(model.set[[i]][[j]]$U[[k]], model.set[[i]][[j]]$V[[k]])     
    }
    
    

    
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
    
    
    #########################################
    ### PSI
    #########################################
    
    
    
    PSI.loading = result.total[[i]][[j]]$PSI.result$PSI.loading
    PSI.structure = sort.columns(result.total[[i]][[j]]$PSI.result$PSI.struct, K,
                                          transpose = TRUE)$structure
    PSI.score = result.total[[i]][[j]]$PSI.result$PSI.score
    
    
    
    if (ncol(PSI.loading) > 0 && ncol(PSI.score) > 0) {
      
      
      
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(length(true.r.list))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,true.r.list[[kkk]] ])$u, svd(PSI.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,true.r.list[[kkk]]])$u, svd(PSI.score)$u)))
      } 
      
      
      
      
      
      
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        
        if (length(which(PSI.structure[k,] != 0)) > 0) {
          
          for(kkk in 1:(length(true.r.list.part[[k]]))) {
            
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ] )$u, 
                                    svd(PSI.loading[p.list[[k]],  which(PSI.structure[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ]  )$u, 
                                            svd(PSI.loading[p.list[[k]], which(PSI.structure[k,] != 0)])$u)))
          }
        }
      }
      
      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]-PSI.loading[p.list.c,] %*% t(PSI.score))^2) / sum(Z[[k]]^2) 
      }
      RSE = RSE / K
      
      
      PSI.var[j, 4] = mean(angle.loading) 
      PSI.var[j, 3] = mean(angle.score)  
      PSI.var[j, 2] = RSE
      PSI.var[j, 1] = identical(trueStructure, PSI.structure )
    } 
    else {
      PSI.var[j,] = rep(NA, mm)
    }
    
    
    #########################################
    ### SLIDE
    #########################################
    
    
    SLIDE.bcv.loading = result.total[[i]][[j]]$SLIDE.result$SLIDE.bcv.loading
    SLIDE.bcv.structure = sort.columns(result.total[[i]][[j]]$SLIDE.result$SLIDE.bcv.struct, K,
                                       transpose = TRUE)$structure
    SLIDE.bcv.score = result.total[[i]][[j]]$SLIDE.result$SLIDE.bcv.score
    
    
    

    
    if (ncol(SLIDE.bcv.loading) > 0 && ncol(SLIDE.bcv.score) > 0) {
      
      
      
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(length(true.r.list))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,true.r.list[[kkk]] ])$u, svd(SLIDE.bcv.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,true.r.list[[kkk]]])$u, svd(SLIDE.bcv.score)$u)))
      } 
      
      
      
      
      
      
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        
        if (length(which(SLIDE.bcv.structure[k,] != 0)) > 0) {
          
          for(kkk in 1:(length(true.r.list.part[[k]]))) {
            
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ] )$u, 
                                    svd(SLIDE.bcv.loading[p.list[[k]],  which(SLIDE.bcv.structure[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ]  )$u, 
                                            svd(SLIDE.bcv.loading[p.list[[k]], which(SLIDE.bcv.structure[k,] != 0)])$u)))
          }
        }
      }
      
      

      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]  - SLIDE.bcv.loading[p.list.c,] %*% t(SLIDE.bcv.score))^2) /sum(Z[[k]]^2)
      }
      RSE = RSE / K
      
      
      SLIDE.bcv.var[j, 4] = mean(angle.loading) 
      SLIDE.bcv.var[j, 3] = mean(angle.score)
      SLIDE.bcv.var[j, 2] = RSE
      SLIDE.bcv.var[j, 1] = identical(trueStructure, SLIDE.bcv.structure )
    } 
    else {
      SLIDE.bcv.var[j,] = rep(NA, mm)
    }
    
    
    #########################################
    ### COBS
    #########################################
    
    COBS.loading = result.total[[i]][[j]]$COBS.result$COBS.loading
    COBS.structure = sort.columns(result.total[[i]][[j]]$COBS.result$COBS.struct, K,
                                  transpose = TRUE)$structure
    COBS.score = result.total[[i]][[j]]$COBS.result$COBS.score
    
    
    
    if (ncol(COBS.loading) > 0 && ncol(COBS.score) > 0) {
      
      
      
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(length(true.r.list))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,true.r.list[[kkk]] ])$u, svd(COBS.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,true.r.list[[kkk]]])$u, svd(COBS.score)$u)))
      } 
      
      
      
      
      
      
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        
        if (length(which(COBS.structure[k,] != 0)) > 0) {
          
          for(kkk in 1:(length(true.r.list.part[[k]]))) {
            
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ] )$u, 
                                    svd(COBS.loading[p.list[[k]],  which(COBS.structure[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ]  )$u, 
                                            svd(COBS.loading[p.list[[k]], which(COBS.structure[k,] != 0)])$u)))
          }
        }
      }
      
      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]-COBS.loading[p.list.c,] %*% t(COBS.score))^2) / sum(Z[[k]]^2)
      }
      RSE = RSE / K
      
      
      COBS.var[j, 4] = mean(angle.loading)  
      COBS.var[j, 3] = mean(angle.score)
      COBS.var[j, 2] = RSE
      COBS.var[j, 1] = identical(trueStructure, COBS.structure )
    } 
    else {
      COBS.var[j,] = rep(NA, mm)
    }
    
    
    #########################################
    ### AJIVE
    #########################################
    
    AJIVE.loading = result.total[[i]][[j]]$AJIVE.result$AJIVE.loading
    AJIVE.structure = sort.columns(result.total[[i]][[j]]$AJIVE.result$AJIVE.struct, K,
                                   transpose = TRUE)$structure
    AJIVE.score = result.total[[i]][[j]]$AJIVE.result$AJIVE.score
    
    
    
    if (ncol(AJIVE.loading) > 0 && ncol(AJIVE.score) > 0) {
      
      
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(length(true.r.list))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,true.r.list[[kkk]] ])$u, svd(AJIVE.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,true.r.list[[kkk]]])$u, svd(AJIVE.score)$u)))
      } 
      
      
      
      
      
      
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        
        if (length(which(AJIVE.structure[k,] != 0)) > 0) {
          
          for(kkk in 1:(length(true.r.list.part[[k]]))) {
            
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ] )$u, 
                                    svd(AJIVE.loading[p.list[[k]],  which(AJIVE.structure[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,true.r.list.part[[k]][[kkk]] ]  )$u, 
                                            svd(AJIVE.loading[p.list[[k]], which(AJIVE.structure[k,] != 0)])$u)))
          }
        }
      }
      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]-AJIVE.loading[p.list.c,] %*% t(AJIVE.score))^2) / sum(Z[[k]]^2)
      }
      RSE = RSE / K
      
      AJIVE.var[j, 4] = mean(angle.loading) 
      AJIVE.var[j, 3] = mean(angle.score)
      AJIVE.var[j, 2] = RSE
      AJIVE.var[j, 1] = identical(trueStructure, AJIVE.structure )
    } 
    else {
      AJIVE.var[j,] = rep(NA, mm)
    }
    
    
    #########################################
    ### JIVE
    #########################################
    
    
    JIVE.loading = result.total[[i]][[j]]$JIVE.result$JIVE.loading
    JIVE.structure = sort.columns(result.total[[i]][[j]]$JIVE.result$JIVE.struct, K,
                                  transpose = TRUE)$structure
    JIVE.score = result.total[[i]][[j]]$JIVE.result$JIVE.score
    
    
    
    if (ncol(JIVE.loading) > 0 && ncol(JIVE.score) > 0) {
      
      angle.loading = c()
      angle.loading.max = c()
      for(k in 1:K){
        if (length(which(JIVE.structure[k,] != 0)) > 0) {
          for(kkk in 1:(ncol(trueU[[k]]))) {
            angle.loading = c(angle.loading, 
                              p.ang(svd(trueU[[k]][,kkk]  )$u, 
                                    svd(JIVE.loading[p.list[[k]], which(result.total[[i]][[j]]$JIVE.result$JIVE.struct[k,] != 0) ])$u))
            angle.loading.max = c(angle.loading.max, 
                                  max(p.ang(svd(trueU[[k]][,kkk]  )$u, 
                                            svd(JIVE.loading[p.list[[k]], which(result.total[[i]][[j]]$JIVE.result$JIVE.struct[k,] != 0)])$u)))
          }
        }
      }
      angle.score = c()
      angle.score.max = c()
      for(kkk in 1:(ncol(trueStructure))) {
        angle.score = c(angle.score, p.ang(svd(trueV[,kkk])$u, svd(JIVE.score)$u))
        angle.score.max = c(angle.score.max, max(p.ang(svd(trueV[,kkk])$u, svd(JIVE.score)$u)))
      } 
      
      
      RSE = 0
      for(k in 1:K){
        p.list.c = model.set[[i]][[j]]$p.list[[k]]
        RSE = RSE + sum((Z[[k]]-JIVE.loading[p.list.c,] %*% t(JIVE.score))^2) / sum(Z[[k]]^2)
      }
      RSE = RSE / K
      
      JIVE.var[j, 4] = mean(angle.loading)  
      JIVE.var[j, 3] = mean(angle.score)
      JIVE.var[j, 2] = RSE
      JIVE.var[j, 1] = identical(trueStructure, JIVE.structure )
    } 
    else {
      JIVE.var[j,] = rep(NA, mm)
    }
    
    
  }
  
  
  
  
  
  ### Put together the mean and sd of each measure
  
  
  loading.grass.angle.mean.mat[i,1] = mean(PSI.var[,4])
  loading.grass.angle.sd.mat[i,1] = sd(PSI.var[,4])
  
  score.grass.angle.mean.mat[i,1] = mean(PSI.var[,3])
  score.grass.angle.sd.mat[i,1] = sd(PSI.var[,3])
  
  str.mean.mat[i,1] = mean(PSI.var[,2])
  str.sd.mat[i,1] = sd(PSI.var[,2])
  
  true.str.mean.mat[i,1] = mean(PSI.var[,1])
  true.str.sd.mat[i,1] = sd(PSI.var[,1])
  
  
  
  
  loading.grass.angle.mean.mat[i,2] = mean(SLIDE.bcv.var[,4])
  loading.grass.angle.sd.mat[i,2] = sd(SLIDE.bcv.var[,4])
  
  score.grass.angle.mean.mat[i,2] = mean(SLIDE.bcv.var[,3])
  score.grass.angle.sd.mat[i,2] = sd(SLIDE.bcv.var[,3])
  
  str.mean.mat[i,2] = mean(SLIDE.bcv.var[,2])
  str.sd.mat[i,2] = sd(SLIDE.bcv.var[,2])
  
  true.str.mean.mat[i,2] = mean(SLIDE.bcv.var[,1])
  true.str.sd.mat[i,2] = sd(SLIDE.bcv.var[,1])
  
  
  
  
  
  
  loading.grass.angle.mean.mat[i,3] = mean(COBS.var[,4])
  loading.grass.angle.sd.mat[i,3] = sd(COBS.var[,4])
  
  score.grass.angle.mean.mat[i,3] = mean(COBS.var[,3])
  score.grass.angle.sd.mat[i,3] = sd(COBS.var[,3])
  
  str.mean.mat[i,3] = mean(COBS.var[,2])
  str.sd.mat[i,3] = sd(COBS.var[,2])
  
  true.str.mean.mat[i,3] = mean(COBS.var[,1])
  true.str.sd.mat[i,3] = sd(COBS.var[,1])
  
  
  
  
  
  
  loading.grass.angle.mean.mat[i,4] = mean(AJIVE.var[,4])
  loading.grass.angle.sd.mat[i,4] = sd(AJIVE.var[,4])
  
  score.grass.angle.mean.mat[i,4] = mean(AJIVE.var[,3])
  score.grass.angle.sd.mat[i,4] = sd(AJIVE.var[,3])
  
  str.mean.mat[i,4] = mean(AJIVE.var[,2])
  str.sd.mat[i,4] = sd(AJIVE.var[,2])
  
  true.str.mean.mat[i,4] = mean(AJIVE.var[,1])
  true.str.sd.mat[i,4] = sd(AJIVE.var[,1])
  
  
  
  
  
  
  loading.grass.angle.mean.mat[i,5] = mean(JIVE.var[,4])
  loading.grass.angle.sd.mat[i,5] = sd(JIVE.var[,4])
  
  score.grass.angle.mean.mat[i,5] = mean(JIVE.var[,3])
  score.grass.angle.sd.mat[i,5] = sd(JIVE.var[,3])
  
  str.mean.mat[i,5] = mean(JIVE.var[,2])
  str.sd.mat[i,5] = sd(JIVE.var[,2])
  
  true.str.mean.mat[i,5] = mean(JIVE.var[,1])
  true.str.sd.mat[i,5] = sd(JIVE.var[,1])
  
}


names = c("PSI", "SLIDE bcv", "COBS", "AJIVE", "JIVE")

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
