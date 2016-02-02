par(mfrow = c(3,2))

plot(results[,1], type="l", ylab = "Accuracy A", xlab = "MonteCarlo")
plot(results[,2], type="l", ylab = "Sparsity", xlab = "MonteCarlo")
plot(results[,3], type="l", ylab = "l2", xlab = "MonteCarlo")
plot(results[,4], type="l", ylab = "Frobenius", xlab = "MonteCarlo")
plot(results[,5], type="l", ylab = "MSE", xlab = "MonteCarlo")
