source("./Other_methods/SLIDE/SLIDE.R")
source("./Other_methods/SLIDE/SLIDE.parse.R")


load(file="BCMO2017_CLL.RData")
data = list(BCMO2017$drugs, BCMO2017$methylation,  BCMO2017$expression)




# Concatenate the data and form pvec
pvec <- c(310,5000,5000)
Xall <- cbind(data[[1]],data[[2]],data[[3]])


# Do group lasso with structure first

out_s <- standardizeX(Xall, pvec, center = T)
X <- out_s$X
svec <- out_s$svec

# Form the list of candidate structures
nl = 70

time1 = Sys.time()

fit_seq<- solve_optim1_seq_restarts(X, lambda_seq = NULL, pvec = pvec, k_max = 1000, eps = 1e-8, reduced = F, rank_total = NULL, lambda_max = max(svec), lambda_min = 0.05)

time2 = Sys.time()

out_struct <- get_structure_v2(fit_seq, pvec)

time3 = Sys.time()

# Select the best structure from the list
outbcv <- bcv_optim1_structure_centering(X, pvec = pvec, structure_list = out_struct$Slist, n_fold = 3, p_fold=3, k_max = 2000, eps = 1e-8, center = T)

time4 = Sys.time()


save(fit_seq, out_struct, outbcv, file="./Real_Data_Analysis_BCMO2017_CLL/Results/result_SLIDE_CLL.Rda")


