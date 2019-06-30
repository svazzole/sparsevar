#' @title Matrix plot
#' 
#' @description Plot a sparse matrix
#' 
#' @param M the matrix to plot
#' @param colors dark or light
#' @return An \code{image} plot with a particular color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotMatrix(M, colors)
#' 
#' @export
plotMatrix <- function(M, colors = "dark") {
  
  if (!is.matrix(M)) { 
    stop("Input must be a matrix") 
  }
  
  nr <- nrow(M)
  nc <- ncol(M)
  M <- t(M)[, nr:1]
  if (colors == "dark") {
    ggplot2::ggplot(reshape2::melt(M), ggplot2::aes_string(x='Var1', y='Var2', fill='value')) + ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(low='red', high='green', mid='black') + ggplot2::xlab("Row") + ggplot2::ylab("Col") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1,vjust=1))    
  } else if (colors == "light") {
    ggplot2::ggplot(reshape2::melt(M), ggplot2::aes_string(x='Var1', y='Var2', fill='value')) + ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(low='red', high='blue', mid='white') + ggplot2::xlab("Row") + ggplot2::ylab("Col") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1,vjust=1))    
  } else {
    stop("Colors must be\"light\" or \"dark\".")
  }
   
  
}

#' @title Plot VARs
#' 
#' @description Plot all the matrices of a VAR model
#' 
#' @param ... a sequence of VAR objects (one or more
#' than one, as from \code{simulateVAR} or \code{fitVAR})
#' @param colors the gradient used to plot the matrix. It can be "light" (low =
#' red -- mid = white -- high = blue) or "dark" (low = red -- mid = black -- 
#' high = green)
#' @return An \code{image} plot with a specific color palette 
#' @usage plotVAR(..., colors)
#' 
#' @export
plotVAR <- function(..., colors = "dark") {
  
  vars <- list(...)
  l <- length(vars)
  
  for (i in 1:l) {
    if (!checkIsVar(vars[[i]])) {
      stop("Inputs must be var objects")
    } 
  }
  
  pl <- list()
  varorder <- length(vars[[1]]$A)
  differentVarOrder <- FALSE
  for (i in 1:l) {
    if (varorder != length(vars[[i]]$A)) {
      differentVarOrder <- TRUE
      varorder <- min(varorder, length(vars[[i]]$A))
    }  
  }
  
  if (differentVarOrder == TRUE) {
    warning("Different VAR orders: plotting up to the min one")
  }
  
  for (i in 1:l) {
    for (j in 1:varorder) {
      pl[[((i-1)*varorder)+j]] <- plotMatrix(vars[[i]]$A[[j]], colors = colors)
    }
  }
  
  multiplot(plotlist = pl, cols = varorder, layout = matrix(1:(l*varorder), nrow = l, byrow = TRUE))
  
}

#' @title Plot VECMs
#' 
#' @description Plot all the matrices of a VECM model
#' 
#' @param v a VECM object (as from \code{fitVECM})
#' @return An \code{image} plot with a specific color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotVECM(v)
#' 
#' @export
plotVECM <- function(v) {
  
  if (attr(v, "class") != "vecm") {
    stop("v must be a VECM object")
  }

  l <- length(v$G)
  pl <- list()
  
  pl[[1]] <- plotMatrix(v$Pi)
  
  if (l > 0) {
    for (i in 1:l) {
      pl[[i+1]] <- plotMatrix(v$G[[i]])
    }
  }
  
  multiplot(plotlist = pl, cols = l+1, layout = matrix(1:(l+1), nrow = 1, byrow = TRUE))
  
}

#' @title Multiplots with ggplot
#' 
#' @description Multiple plot function. ggplot objects can be passed in ..., or 
#' to plotlist (as a list of ggplot objects)
#' @param ... a sequence of ggplots to be plotted in the grid.
#' @param plotlist a list containing ggplots as elements.  
#' @param cols number of columns in layout
#' @param layout a matrix specifying the layout. If present, 'cols' is ignored.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' Taken from R Cookbook
#' 
#' @return A ggplot containing the plots passed as arguments 
#' @export
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
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
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

