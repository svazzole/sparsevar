library(h5)

file <- h5file(name = "data/1987.h5", mode = "r")

rets <- file["rets"][]

h5close(file)


library(rhdf5)

rets <- t(h5read("data/1987.h5","rets/block0_values"))
rets <- rets[3:nrow(rets), ]

res <- estimateVAR(rets, options = list(parallel = TRUE, ncores = 2))

source("estimateVAR.R")

for (y in 1926:2012) {
  
  data <- paste0("data/",as.character(y),".h5")
  
  rets <- t(h5read(data,"rets/block0_values"))
  
  rets <- rets[3:nrow(rets), ]
  
  res <- estimateVAR(rets, options = list(parallel = TRUE, ncores = 2))
  
  fileN <- paste0("fig/",as.character(y),".png")
  png(filename = fileN)
  image(res$A[100:1, ])
  dev.off()
  
}
