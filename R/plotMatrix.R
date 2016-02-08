#' @title Matrix plot
#' 
#' @description Plot a sparse matrix
#' 
#' @param M the matrix to plot
#' @return An \code{image} plot with a particular color palette
#' @examples
#' plotMatrix(M)
#' 
#' @author Simone Vazzoler
#'
#' @export
#' 
plotMatrix <- function(M){
  
  nr <- nrow(M)
  nc <- ncol(M)
  colfunc<-colorRampPalette(c("red","white","royalblue"))  
  image(1:nr, 1:nc, M[, nr:1], col=colfunc(39), xlab = "Row number", ylab = "Column number")
  
}

