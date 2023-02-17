
source("./Source/Algorithm.R")
source("./Source/general_functions.R")
source("./Source/Estimation_fun.R")
source("./Source/Reconstruction.R")
source("visualization.R")

source("./Other_methods/COBS/COBS.parse.R")

library(dplyr)

load(file="./Real_Data_Analysis_BCMO2017_CLL/BCMO2017_CLL.RData")
data = list(BCMO2017$drugs, BCMO2017$methylation,  BCMO2017$expression)



dataY <- cbind(scale(data[[1]],center=T,scale=T),
               scale(data[[2]],center=T,scale=T),
               scale(data[[3]],center=T,scale=T))

# Retrieve the dimensions of the data matrices
n=dim(dataY)[1]; p=dim(dataY)[2]
gp=ceiling(1:p/200) #group index of Y, multi-block 
pvec <- vector()  
for (i in 1:length(unique(gp))) {
  pvec[i] <- sum( gp == unique(gp)[i] ) 
}

K = length(data)
pp = c(dim(data[[1]])[2],dim(data[[2]])[2],dim(data[[3]])[2])
p.list = get.list(pp,K)

n.comp = 5 + 42 +3


time1 = Sys.time()

COBS.tune.result = cobs.tune(dataY = dataY, Y.bk.idx = gp,
                             n.comp = 5 + 42 +3,
                             alpha.v = 0.5,
                             method = "BIC-Low",
                             M.lambda.v = 40,
                             orth=FALSE)

time2 = Sys.time()

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

# save(COBS.tune.result, file="./Real_Data_Analysis_BCMO2017_CLL/Results/result_COBS_CLL.Rdata")
