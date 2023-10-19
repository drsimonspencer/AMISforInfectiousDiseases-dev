#' Plot weighted simulated prevalences for a given an object of class \code{amis}.
#'
#' @param x The output from the function \code{\link{amis}()}.
#' @param location Integer identifying the location. Default to 1.
#' @param time Integer identifying the timepoint. Default to 1.
#' @param breaks Argument passed to \code{\link{wtd.hist}()}. Default to 500.
#' @param xlim Argument passed to \code{\link{wtd.hist}()}. Default to NULL.
#' @param main Title for the plot.
#' @param ... Other graphical parameters passed to \code{\link{wtd.hist}()}.
#' @importFrom  weights wtd.hist
#' @return A plot.
#' @export
plot.amis <- function(x, location=1, time=1, main=NULL, 
                      breaks=500, xlim=NULL, ...){
  if(!inherits(x, "amis")){
    stop("'x' must be of type amis")
  }
  amis_params <- x$amis_params
  weights <- x$weight_matrix
  sim_prev <- x$simulated_prevalences[,time]

  if(amis_params$log){weights <- exp(weights)}
  
  if(is.null(xlim)){
    if(all(is.finite(amis_params$boundaries))){
      xlim <- amis_params$boundaries
    }else{
      xlim <- range(sim_prev)
    }
  }
  if(is.null(main)){
    main <- paste0("Location ", location, " at time ", time)
  }
  weights::wtd.hist(x=sim_prev, breaks=breaks, 
                    weight=weights[,location],
                    probability=T, xlim=xlim,
                    xlab="Weighted prevalence",
                    main=main, ...)
  
}
