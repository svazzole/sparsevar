#' @title IRF plot
#' 
#' @description Plot a IRF object
#' 
#' @param irf the matrix to plot
#' @return An \code{image} plot with a particular color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotIRF(M)
#' 
#' @export
plotIRF <- function(irf, bb, i, j, type = "irf") {
  
  if (attr(irf, "class") != "irf" | attr(bb, "class") != "irfBands") { 
    stop("Inputs must be an irf object and an irfBands object") 
  }
  
  nz <- dim(irf$irf)[3]
  t <- 0:(nz-1)
  
  if (type == "irf") {
    
    irfString <- paste0("IRF ", j, " -> ", i)
    ub <- irf$irf[i,j,] + bb$irfUB[i,j,]
    lb <- irf$irf[i,j,] + bb$irfLB[i,j,]
    d <- as.data.frame(cbind(t, irf$irf[i,j,], lb, ub, 0))
    
  } else if (type == "oirf") {
    
    irfString <- paste0("OIRF ", j, " -> ", i)
    ub <- irf$oirf[i,j,] + bb$oirfUB[i,j,]
    lb <- irf$oirf[i,j,] + bb$oirfLB[i,j,]
    d <- as.data.frame(cbind(t, irf$oirf[i,j,], lb, ub, 0))
    
  } else {
    stop("Unknown type")
  }

  ggplot2::ggplot(d, ggplot2::aes(x = d[,1], y = d[,2])) + ggplot2::geom_line() + ggplot2::xlab("Time") + ggplot2::ylab(irfString) +
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,3]), linetype = "dashed", color = "blue") + 
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,4]), linetype = "dashed", color = "blue") + 
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,5]), color = "red")
  
}

plotIRFGrid <- function(irf, bb, indexes, type = "irf") {
  
  n <- length(indexes)
  g <- expand.grid(indexes, indexes)
  nrgrid <- nrow(g)
  
  pl <- list()
  
  for (i in 1:nrgrid) {
    pl[[i]] <- plotIRF(irf, bb, g[i,1], g[i,2], type = type)
  }
  
  #title("Simulation")
  multiplot(plotlist = pl, cols = n, layout = matrix(1:nrgrid, nrow = n, byrow = TRUE))
  
}