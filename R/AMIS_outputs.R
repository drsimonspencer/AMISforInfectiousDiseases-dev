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

#' Print method for object of class \code{amis}
#'
#' @param x The output from the function \code{\link{amis}()}.
#' @param ... Additional printing options.
#' @return Brief description of data and model specifications used to run \code{\link{amis}()}.
#' @export
print.amis <- function(x, ...) {
  
  if(!inherits(x, "amis")){
    stop("'x' must be of type amis")
  }
  amis_params <- x$amis_params
  cat(paste0("=============================================================", "\n"))
  cat(paste0("Data description: \n"))
  cat(paste0("- Number of locations:  ", nrow(x$prevalence_map[[1]]$data),"\n"))
  cat(paste0("- Number of map samples for each location:  ",  ncol(x$prevalence_map[[1]]$data),"\n"))
  cat(paste0("- Number of time points:  ",  length(x$prevalence_map),"\n"))
  locs_no_data <- x$locations_with_no_data
  if(!is.null(locs_no_data)){
    if(length(locs_no_data)==1L){
      message(paste0("- Location ", shQuote(locs_no_data), " has no valid map samples. Its corresponding posterior samples were therefore driven by the prior."))
    }else{
      message(paste0("- Locations (", paste(shQuote(locs_no_data), collapse=","), ") have no valid map samples. Their corresponding posterior samples were therefore driven by the prior."))
    }
  }

  cat(paste0("=============================================================", "\n"))
  cat(paste0("Model and algorithm specifications: \n"))
  cat(paste0("- For the nonparametric estimation of the density of the likelihood: \n"))
  if(!is.null(amis_params[["breaks"]])){
    cat("     Histogram method was used with breaks supplied by the user. \n")
  }else{
    if(!is.null(amis_params[["sigma"]])){
      cat(paste0("     Gaussian kernel was used with smoothing parameter sigma = ",as.character(round(amis_params$sigma, digits = getOption("digits"))), "\n"))
    }else{
      cat(paste0("     Uniform kernel was used with smoothing parameter delta = ",as.character(round(amis_params$delta, digits = getOption("digits"))), "\n"))
    }
  }
  cat(paste0("- Lower and upper boundaries for prevalences:  ", 
             paste0(amis_params$boundaries, collapse=", "),"\n"))
  cat(paste0("- Number of new samples drawn within each AMIS iteration:  ",  amis_params$nsamples,"\n"))
  cat(paste0("- Maximum number of iterations:  ",  amis_params$max_iters,"\n"))
  cat(paste0("- Target effective sample size:  ",  amis_params$target_ess,"\n"))
  
}
  