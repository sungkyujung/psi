

######################################################################
# Draw Table 3 (Joint-emphasized) and Table 4 (Individual-emphasized)
######################################################################

######################################################################
# Needs both "model.set" and "result.total" variables
# from Sim_imbalance_(Sec5.3).R
# (or "joint.data.Rdata" and "joint.Rdata")
# (or "ind.data.Rdata" and "ind.Rdata")
######################################################################



library(xtable)
source("./Source/general_functions.R")
source("./Source/parser.R")
###

ncl = 5       #number of methods
mm = 3     #number of measures
nrep = 100


SNR.list = c(0,10)

res = array(0, c(length(SNR.list), 3, ncl, nrep ))

for(i in 1:length(SNR.list)){
  for(j in 1:nrep){
    str = list()
    str[[1]] = result.total[[i]][[j]]$PSI.result$PSI.struct
    str[[2]] = result.total[[i]][[j]]$SLIDE.result$SLIDE.bcv.struct
    str[[3]] = result.total[[i]][[j]]$COBS.result$COBS.struct
    str[[4]] = result.total[[i]][[j]]$AJIVE.result$AJIVE.struct
    str[[5]] = result.total[[i]][[j]]$JIVE.result$JIVE.struct
    
    trueStructure = model.set[[i]][[j]]$true.struct
    
    for(l in 1:ncl){
      
      if (identical(trueStructure, str[[l]] )) {res[i,1,l,j] = 1}
      
      r = ncol(str[[l]])
      for(cn in 1:r){
        if (identical(str[[l]][,cn], c(1,1,1) )) { res[i,2,l,j] = res[i,2,l,j] + 1 }
        if (identical(str[[l]][,cn], c(1,0,0) )) { res[i,3,l,j] = res[i,3,l,j] + 1 }
        if (identical(str[[l]][,cn], c(0,1,0) )) { res[i,3,l,j] = res[i,3,l,j] + 1 }
        if (identical(str[[l]][,cn], c(0,0,1) )) { res[i,3,l,j] = res[i,3,l,j] + 1 }
        
      }
    }
  }

  

}





names = c("PSI", "SLIDE bcv", "COBS", "AJIVE", "JIVE")


total.mean.mat = matrix(0, mm * length(SNR.list), ncl)
total.sd.mat = matrix(0, mm * length(SNR.list), ncl)


for(i in 1:length(SNR.list)) {
  for(l in 1:ncl){
  total.mean.mat[ mm * (i-1) + 3 , l] = mean(res[i,3,l,])
  total.mean.mat[ mm * (i-1) + 2 , l] = mean(res[i,2,l,])
  total.mean.mat[ mm * (i-1) + 1 , l] = mean(res[i,1,l,]) * 100
  
  total.sd.mat[ mm * (i-1) + 3 , l] = sd(res[i,3,l,])
  total.sd.mat[ mm * (i-1) + 2 , l] = sd(res[i,2,l,])
  total.sd.mat[ mm * (i-1) + 1 , l] = sd(res[i,1,l,]) * 100
  }
}

msd <- paste(round(total.mean.mat,2)," (",round(total.sd.mat,2),")",sep="")

tab.msd <- matrix(msd, mm * length(SNR.list), ncl)
colnames(tab.msd) = names


Measure = rep(c("Accuracy  ", "Joint  ", "Individual"), length(SNR.list))


SNR = c()
for (i in 1:length(SNR.list)) {
  SNR = c(SNR, as.character(SNR.list[i]))
  SNR = c(SNR, rep(" ", mm - 1))
}


tab.msd = cbind(SNR, Measure, tab.msd)

print(xtable(tab.msd),  include.rownames=FALSE)
