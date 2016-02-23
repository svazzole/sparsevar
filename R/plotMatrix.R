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
  M <- Matrix(M, sparse = TRUE)
  image(M)
  #image(1:nr, 1:nc, t(M)[,nc:1], col=colfunc(39), xlab = "Column Number", ylab = "Row number", yaxt = "n", xaxt = "n")
  #axis(1, at = 1:nc, labels = as.character(1:nc)) 
  #axis(2, at = 1:nc, labels = as.character(nc:1), las = 2) 
  #axis(2, at = 1:nc, labels = as.character(nc:1)) 
  
}

