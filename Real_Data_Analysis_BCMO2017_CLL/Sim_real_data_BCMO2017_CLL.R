
###############################################################################
# Source code for Section 6
# Real Data Analysis
###############################################################################





source("./Source/Algorithm.R")
source("./Source/general_functions.R")
source("./Source/Estimation_fun.R")
source("./Source/Reconstruction.R")
source("visualization.R")

load(file="./Real_Data_Analysis_BCMO2017_CLL/BCMO2017_CLL.RData")
data = list(BCMO2017$drugs, BCMO2017$methylation,  BCMO2017$expression)



### Data Preprocessing


X = list()
scale.param = list()

for(i in 1:length(data)) {
  scale.X = scale(data[[i]])
  X[[i]] = t(scale.X)
  
  scale.center = attr(scale.X,"scaled:center")
  scale.scale = attr(scale.X,"scaled:scale")
  scale.param[[i]] = list(scale.center, scale.scale)
}

Data.name = c("Drug Viability", "Methylation", "Expression")





### Cut off the Signal Ranks (60% cuts of variations)

cut.off = 60.0
cut.upper = c()
cut.lower = c()
for(i in 1:length(data)) {
  data.pca = prcomp(data[[i]])
  cum.var = cumsum(100 * data.pca$sdev^2 / sum(data.pca$sdev^2))

  if (cut.off < 100 ){
    cut.upper = c(cut.upper, min(which((cum.var > cut.off ))) )
  } else {
    cut.upper = c(cut.upper, min(nrow(X[[i]]), ncol(X[[i]]) )  )
  }
  cut.lower = c(cut.lower, min(which((cum.var > 00.00))))
}







### Scree plots

scree.plot(data[[1]], Data.name[[1]], cut.upper[1])
scree.plot(data[[2]], Data.name[[2]], cut.upper[2])
scree.plot(data[[3]], Data.name[[3]], cut.upper[3])










### Setting Parameters


n = nrow(data[[1]]); K = length(data)
p = c()
for(i in 1:length(data)) {
  p  = c(p, nrow(X[[i]]))
}
r.list = cut.upper - cut.lower + rep(1,K)
lambda = seq(0,90,1) * pi / 180
model.list = list(n = n,K = K,p = p, r.list = r.list, X = X)





### Main Algorithm

library(doParallel)
library(doRNG)
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl) 


set.seed(1001)

Nrep = 110; nrep = 100

Sys.time()

cat( "   time = ",system.time(
  result.parallel.val  <- foreach(j=1:Nrep,.packages = c("irlba", "RSpectra")) %dopar% {
    estimation.fun(model.list,  lambda, set.training =  sample(n, ceiling(n/2)), light.ver = TRUE, r.list)
  }),"\n")
result.parallel = result.parallel.val[1:nrep]

Sys.time()

stopCluster(cl)





### The Estimated Structure by finding Mode Structure

mode.result = mode.structure(result.parallel)
mode.str = mode.result$structure[[mode.result$model.mode]]
mode.str

### Reconstructed Signal Matrix for each Index-set

pc.result = reconstruction.partial.cluster.mode(model.list,mode.str, scale.param, lambda) 
pc.matrix = pc.result$pc.matrix
pc.list = pc.result$pc



