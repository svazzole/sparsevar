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
plotMatrix <- function(M) {
  
  nr <- nrow(M)
  nc <- ncol(M)
  M <- t(M)[, nc:1]
  ggplot(melt(M), aes(Var1,Var2, fill=value)) + geom_raster() + scale_fill_gradient2(low='red', high='green', mid='black') + xlab("Row") + ylab("Col")
  
  # Old plot style
  # nr <- nrow(M)
  # nc <- ncol(M)
  # colfunc<-colorRampPalette(c("red","white","royalblue"))  
  # # M <- Matrix(M, sparse = TRUE)
  # # image(M)
  # image(1:nr, 1:nc, t(M)[,nc:1], col=colfunc(39), xlab = "Column Number", ylab = "Row number", yaxt = "n", xaxt = "n")
  # axis(1, at = 1:nc, labels = as.character(1:nc))
  # # axis(2, at = 1:nc, labels = as.character(nc:1), las = 2)
  # axis(2, at = 1:nc, labels = as.character(nc:1))
  
}

#' @export
plotVAR <- function(A) {
  
  p <- length(A)
  
  pl <- list()
  # par(mfrow = c(1,p))
  for (i in 1:p) {
    pl[[i]] <- plotMatrix(A[[i]])
  }
  
  multiplot(plotlist = pl, cols = p)
}

#' @export
plotComparisonVAR <- function(est, sim) {
  
  p <- length(sim$A)
  
  pl <- list()
  # par(mfrow = c(2,p))
  for (i in 1:p) {
    pl[[i]] <- plotMatrix(est$A[[i]])
  }
  #title("Estimate")
  for (i in 1:p) {
    pl[[i + p]] <- plotMatrix(sim$A[[i]])
  }
  #title("Simulation")
  multiplot(plotlist = pl, cols = p, layout = matrix(1:(2*p), nrow = 2, byrow = TRUE))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# From R Cookbook
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}