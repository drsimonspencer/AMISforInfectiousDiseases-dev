#' Plot weighted simulated prevalences for a given an object of class 'amis'.
#'
#' @param x The output from the function \code{\link{amis}()}
#' @param cex.points Graphical parameter
#' @param location Location to be plotted
#' @param lwd.points Graphical parameter
#' @param pch Graphical parameter
#' @param lwd Graphical parameter
#' @param main Title for the plot
#' @param ... Graphical parameters passed to plot().
#' @importFrom  weights wtd.hist
#' @return A plot
#' @export
plot.amis <- function(x, location=1, main=NULL, cex.points=NULL, 
                      lwd.points=NULL, pch=NULL, lwd=NULL, ...){
  if(!inherits(x, "amis")){
    stop("'x' must be of type amis")
  }
  amis_params <- x$amis_params
  obj <- x$sample
  seeds <- obj$seeds
  W <- obj$W
  k <- obj$k
  kprev1 <- obj$prev1
  kprev2 <- obj$prev2
  kprev3 <- obj$prev3
  kprev4 <- obj$prev4
  
  # L <- 
  weights <- obj[,(ncol(obj)-L+1):ncol(obj)]
  
  print(amis_params)
  if(amis_params$log){weights <- exp(weights)}
  weights::wtd.hist(kprev1,breaks=500,weight=weights[,location],
                    probability=T,
                    xlab = "weighted prevalence",
                    main = paste0("Location: ", location), ...)
  
}
