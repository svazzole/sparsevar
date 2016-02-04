#' Plot the matrix
#' 
#' @param M the matrix to plot
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)

plotMatrix <- function(M){
  
  nr <- nrow(M)
  nc <- ncol(M)
  colfunc<-colorRampPalette(c("red","white","royalblue"))  
  image(1:nr, 1:nc, M[, nr:1], col=colfunc(39), xlab = "Row number", ylab = "Column number")
  
}

