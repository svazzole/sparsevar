# # library(h5)
# 
# file <- h5file(name = "data/1987.h5", mode = "r")
# 
# rets <- file["rets"][]
# 
# h5close(file)
# 
# 
# library(rhdf5)
# 
# rets <- t(h5read("data/2011.h5","rets/block0_values"))
# rets <- rets[3:nrow(rets), 1:30]
# prices <- t(h5read("data/2011.h5","prices/block0_values"))
# prices <- prices[3:nrow(prices), 1:30]
# 
# 
# res <- estimateVAR(rets, options = list(parallel = TRUE, ncores = 2))
# 
# # source("estimateVAR.R")
# 
# for (y in 1926:2012) {
# 
#   data <- paste0("data/",as.character(y),".h5")
# 
#   rets <- t(h5read(data,"rets/block0_values"))
# 
#   rets <- rets[3:nrow(rets), ]
# 
#   res <- estimateVAR(rets, options = list(parallel = TRUE, ncores = 2))
# 
#   fileN <- paste0("fig/",as.character(y),".png")
#   png(filename = fileN)
#   image(res$A[100:1, ])
#   dev.off()
# 
# }
# 
# par(mfrow = c(1,2))
# plotMatrix(res$Pi)
# plotMatrix(res$G[[1]])
# 
# #
# resVAR <- estimateVAR(rets)
# resVAR$mse
# 
# res1 <- estimateVECM(prices, options = list(lambda = "lambda.min"))
# res1$mse
# res2 <- estimateVECM(prices, options = list(lambda = "lambda.1se"))
# res2$mse
# 
# par(mfrow = c(2,2))
# plotMatrix(res1$Pi)
# plotMatrix(res1$G[[1]])
# plotMatrix(res2$Pi)
# plotMatrix(res2$G[[1]])
# 
