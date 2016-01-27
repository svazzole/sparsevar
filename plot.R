par(mfrow = c(4,2))

plot(results[,1], type="l", ylab = "Accuracy A", xlab = "MonteCarlo")
plot(results[,2], type="l", ylab = "Accuracy A2", xlab = "MonteCarlo")
plot(results[,3], type="l", ylab = "Sparsity A", xlab = "MonteCarlo")
plot(results[,4], type="l", ylab = "Sparsity A2", xlab = "MonteCarlo")

plot(results[,5], type="l", ylab = "l2(A)", xlab = "MonteCarlo")
plot(results[,6], type="l", ylab = "l2(A2)", xlab = "MonteCarlo")
plot(results[,7], type="l", ylab = "frob(A)", xlab = "MonteCarlo")
plot(results[,8], type="l", ylab = "frob(A2)", xlab = "MonteCarlo")