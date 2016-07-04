#' @title IRF plot
#' 
#' @description Plot a IRF object
#' 
#' @param irf the irf object to plot
#' @param eb the errorbands to plot
#' @param i the first index
#' @param j the second index
#' @param type \code{type = "irf"} or \code{type = "oirf"}
#' @param bands \code{"quantiles"} or \code{"sd"}
#' @return An \code{image} plot relative to the impulse response function.
#' @usage plotIRF(irf, eb, i, j, type, bands)
#' 
#' @export
plotIRF <- function(irf, eb, i, j, type = "irf", bands = "quantiles") {
  
  if (attr(irf, "class") != "irf" | attr(eb, "class") != "irfBands") { 
    stop("Inputs must be an irf object and an irfBands object") 
  }
  
  nz <- dim(irf$irf)[3]
  t <- 0:(nz-1)
  
  ebs <- list()
  
  if (bands == "quantiles") {
    
    ebs$irfUB <- eb$irfQUB[i,j,] - irf$irf[i,j,]
    ebs$irfLB <- eb$irfQLB[i,j,] - irf$irf[i,j,]
    ebs$oirfUB <- eb$oirfQUB[i,j,] - irf$oirf[i,j,]
    ebs$oirfLB <- eb$oirfQLB[i,j,] - irf$oirf[i,j,]
  
  } else if (bands == "sd") {
    
    ebs$irfUB <- eb$irfUB[i,j,]
    ebs$irfLB <- eb$irfLB[i,j,]
    ebs$oirfUB <- eb$oirfUB[i,j,]
    ebs$oirfLB <- eb$oirfLB[i,j,]
    
  } else {
    stop("Possible values for bands are sd or quantiles")
  }
  
  if (type == "irf") {
    
    irfString <- paste0("IRF ", j, " -> ", i)
    ub <- irf$irf[i,j,] + ebs$irfUB
    lb <- irf$irf[i,j,] + ebs$irfLB
    d <- as.data.frame(cbind(t, irf$irf[i,j,], lb, ub, 0))
    
  } else if (type == "oirf") {
    
    irfString <- paste0("OIRF ", j, " -> ", i)
    ub <- irf$oirf[i,j,] + ebs$oirfUB
    lb <- irf$oirf[i,j,] + ebs$oirfLB
    d <- as.data.frame(cbind(t, irf$oirf[i,j,], lb, ub, 0))
    
  } else {
    stop("Unknown type")
  }

  ggplot2::ggplot(d, ggplot2::aes(x = d[,1], y = d[,2])) + ggplot2::geom_line() + ggplot2::xlab("Time") + ggplot2::ylab(irfString) +
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,3]), linetype = "dashed", color = "blue") + 
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,4]), linetype = "dashed", color = "blue") + 
    ggplot2::geom_line(data = d, ggplot2::aes(x = t, y = d[,5]), color = "red")
  
}

#' @title IRF grid plot
#' 
#' @description Plot a IRF grid object
#' 
#' @param irf the irf object computed using impulseResponse
#' @param eb the error bands estimated using errorBands
#' @param indexes a vector containing the indeces that you want to plot
#' @param type plot the irf (\code{type = "irf"} by default) or the orthogonal irf 
#' (\code{type = "oirf"})
#' @return An \code{image} plot relative to the impulse response function.
#' @usage plotIRFGrid(irf, eb, indexes, type)
#' 
#' @export
plotIRFGrid <- function(irf, eb, indexes, type = "irf") {
  
  n <- length(indexes)
  g <- expand.grid(indexes, indexes)
  nrgrid <- nrow(g)
  
  pl <- list()
  
  for (i in 1:nrgrid) {
    pl[[i]] <- plotIRF(irf, eb, g[i,1], g[i,2], type = type)
  }
  
  multiplot(plotlist = pl, cols = n, layout = matrix(1:nrgrid, nrow = n, byrow = TRUE))
  
}