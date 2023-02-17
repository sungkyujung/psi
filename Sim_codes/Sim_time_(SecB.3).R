
library(microbenchmark)


source("./Source/general_functions.R")
source("./Source/parser.R")

source("./Other_methods/SLIDE/SLIDE.bcv.parse.R")
source("./Other_methods/SLIDE/Data_preparations.R")

source("./Other_methods/COBS/COBS.parse.R")

source("./Other_methods/AJIVE/AJIVE.parse.R")

source("./Other_methods/JIVE/JIVE.parse.R")
source("./Other_methods/JIVE/JIVE.data.R")



model.no = 4
save.name1 = paste0("./Sim_result/Sim_comparisons_(Sec5.2)/m",model.no,".data.Rdata")
load(file = save.name1)



SNR.list = c(15, 10, 5)
lambda = seq(0,90,1) * pi / 180



bench.result = list()





for(j in 1:(length(SNR.list))){
  
  model = model.set[[j]][[1]]
  
  
  K = model$K  
  
  ###  
  
  
  
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
  
  model.JIVE = data.prep.JIVE( model, sep = TRUE, r.list = model$r.list)
  
  
  
  
  ###
  
  trueU = matrix(0, sum(model$p), sum(model$r))
  for (k in 1:(model$K)) {
    trueU[model$p.list[[k]], model$r.structure[[k]] ] = model$U[[k]]
  }
  trueV = do.call(cbind, model$V.comp)
  true.structure = model$true.struct
  
  
  lambda = seq(0,90,1) * pi / 180
  
  
  
  
  
  bench.result[[j]] = microbenchmark(
    "JIVE" = JIVE.parse(model, model.JIVE, "perm"),
    "AJIVE" =  AJIVE.parse(XX, model$r.list,
                           n = model$n,
                           p.list = model$p.list, 
                           p = model$p,
                           K = model$K),
    "COBS" = COBS.parse(model.COBS, 
                        n.comp = sum(model$r), 
                        alpha.v = 0.5, 
                        method = "BIC-Low", 
                        M.lambda.v = 40, 
                        orth = F, 
                        p.list = model$p.list, 
                        K = model$K),
    "SLIDE bcv" = SLIDE.bcv.parse(model = model.SLIDE),
    "PSI" = PSI.parse(model, lambda, model$r.list),
    times = 10L
  )
  
}



# save(bench.result, file = paste0("./Sim_result/Sim_time/m",model.no,".time.Rdata"))















###########################




library(xtable)

ncol= 5

time.mean = matrix(0, length(bench.result), ncol)
time.sd = matrix(0,length(bench.result), ncol)


for(i in 1:length(bench.result)){
  NewAlgorithm.time = bench.result[[i]]$time[ bench.result[[i]]$expr == "PSI" ] / (10^9)
  SLIDE.bcv.time = bench.result[[i]]$time[ bench.result[[i]]$expr == "SLIDE bcv" ] / (10^9)
  COBS.time = bench.result[[i]]$time[ bench.result[[i]]$expr == "COBS" ] / (10^9)
  AJIVE.time = bench.result[[i]]$time[ bench.result[[i]]$expr == "AJIVE" ] / (10^9)
  JIVE.time = bench.result[[i]]$time[ bench.result[[i]]$expr == "JIVE" ] / (10^9)
  
  time.mean[i,1] = mean(NewAlgorithm.time)
  time.mean[i,2] = mean(SLIDE.bcv.time)
  time.mean[i,3] = mean(COBS.time)
  time.mean[i,4] = mean(AJIVE.time)
  time.mean[i,5] = mean(JIVE.time)
  
  time.sd[i,1] = sd(NewAlgorithm.time)
  time.sd[i,2] = sd(SLIDE.bcv.time)
  time.sd[i,3] = sd(COBS.time)
  time.sd[i,4] = sd(AJIVE.time)
  time.sd[i,5] = sd(JIVE.time)
}




names = c("PSI", "SLIDE bcr", "COBS", "AJIVE", "JIVE")
SNR = rep(as.character(SNR.list), each = 1)


msd <- paste(round(time.mean,2)," (",round(time.sd,2),")",sep="")

tab.msd <- matrix(msd, length(SNR.list), ncol)
colnames(tab.msd) = names




tab.msd = cbind(SNR, tab.msd)
xtable(tab.msd)












