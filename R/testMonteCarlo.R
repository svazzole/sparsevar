# 
# # source("mcSimulations.R")
# 
# ##
# # ALL
# ##
# 
# g <- expand.grid(nobs = c(50,250), N = c(30,50,100,200), rho = c(0.1,0.5,0.9), sparsity = c(0.05,0.1,0.15), 
#                  penalty = c("ENET", "SCAD", "MCP"), lambda = c("lambda.min", "lambda.1se"), stringsAsFactors = FALSE)
# save(file = "results/ENET/grid.RData", g)
# 
# opt <- list(parallel = TRUE, ncores = 2)
# 
# for (i in 1:nrow(g)) {
#   
#   nobs <- g$nobs[i]
#   N <- g$N[i]
#   rho <- g$rho[i]
#   sparsity <- g$sparsity[i]
#   penalty <- g$penalty[i]
#   
#   opt$lambda <- g$lambda[i] 
#   
#   results <- mcSimulations(N = N, nobs = nobs, rho = rho, sparsity = sparsity, penalty = penalty, options = opt)
#   
#   fileName <- paste0("results/ENET/results_", as.character(i), ".RData")
#   
#   save(file = fileName, results)
#   
# }
# 
