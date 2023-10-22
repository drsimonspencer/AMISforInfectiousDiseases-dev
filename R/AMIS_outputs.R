#' Plot weighted simulated prevalences for a given an object of class \code{amis}.
#'
#' @param x The output from the function \code{\link{amis}()}.
#' @param what What posterior distribution should be plotted. 
#' It can be 'prev' (default) for plotting prevalences, or one of the parameter names. 
#' @param type Type of plot. It can be 'hist' (default) for histogram, 
#' or 'CI' for credible intervals
#' @param locations Integer identifying the location. Default to 1.
#' @param time Integer identifying the timepoint. Default to 1.
#' @param measure_central Measure of central tendency for credible interval plots. 
#' It can be 'mean' (default) or 'median'.
#' @param alpha Numeric value between 0 and 1 indicating the endpoints of the 
#' credible intervals, which are evaluated at (alpha/2, 1-alpha/2)% quantiles. 
#' Default (0.05) will create 95% credible intervals.
#' @param breaks Argument passed to \code{\link{wtd.hist}()} for histogram plots. 
#' Default to 500.
#' @param xlim The x limits of the plots. 
#' Default to NULL.
#' @param main Title for the plot.
#' @param mfrow A vector of the form \code{c(nrows, ncols)} describing the layout 
#' of multiple histogram plots. Default to \code{c(1,1)}.
#' @param ... Other graphical parameters passed to \code{\link{wtd.hist}()}.
#' @importFrom  weights wtd.hist
#' @importFrom graphics segments
#' @importFrom graphics par
#' @return A plot.
#' @export
plot.amis <- function(x, what="prev", type="hist", locations=1, time=1, 
                      measure_central="mean", alpha=0.05, 
                      breaks=500, xlim=NULL, main=NULL, mfrow=c(1,1), ...){
  if(!inherits(x, "amis")){
    stop("'x' must be of type amis")
  }
  
  param_names <- colnames(x$param)
  
  if(!type%in%c("hist","CI")){stop("Argument 'type' must be either 'hist' or 'CI'.")}
  
  if(type=="hist"){
    if(length(what)!=1){
      stop("If type = 'hist', 'what' must have length one")
    }
    if(!what%in%c("prev", param_names)){
      stop("Argument 'what' must be either 'prev' or one of the model parameter names.")
    }
  }else{
    if(!all(what%in%c("prev", param_names))){
      stop("All entries of argument 'what' must be 'prev' or model parameter names.")
    }
  }
  
  n_locs <- length(locations)
  location_names <- colnames(x$weight_matrix)[locations]
  
  # # ---------------------------------------------------
  # Histograms
  if(type=="hist"){
    par(mfrow=mfrow)
    for(location in locations){
      amis_params <- x$amis_params
      weights <- x$weight_matrix[,location]
      if(amis_params$log){weights <- exp(weights)}
      
      
      if(what=="prev"){
        statistic <- x$simulated_prevalences[,time]
      }else{
        if(!(what%in%colnames(x$param))){
          stop("'parameter' must be either 'prev' or have a valid parameter name.")
        }
        statistic <- x$param[,what]
      }


      if(is.null(xlim)){
        if(what=="prev"){
          if(all(is.finite(amis_params$boundaries))){
            xlim <- amis_params$boundaries
          }else{
            xlim <- range(statistic)
          }
        }else{
          xlim <- range(statistic)
        }
      }
      
      if(is.null(main)){
        if(what=="prev"){
          main_ <- paste0("Location '", location, "' at time ", time)
        }else{
          main_ <- paste0("Location '", location, "'")
        }
      }else{
        main_ <- NULL
      }
      hist_title <- ifelse(what=="prev", "Weighted prevalence", 
                           paste0("Weighted ", what))
      weights::wtd.hist(x=statistic, breaks=breaks, 
                        weight=weights,
                        probability=T, xlim=xlim,
                        xlab=hist_title,
                        main=main_, ...)
      
    }
    par(mfrow=c(1,1))
  }
  
  # # ---------------------------------------------------
  # Credible intervals
  if(type=="CI"){
    for(what_ in what){
      summaries <- calculate_summaries(x=x, what=what_, time=1, locations=locations, alpha=alpha)
      if(measure_central=="mean"){
        mu <- summaries[["mean"]]
      }else if(measure_central=="median"){
        mu <- summaries[["median"]]
      }
      names(mu) <- location_names
      lo <- summaries[["quantiles"]][1, ]
      up <- summaries[["quantiles"]][2, ]
      CItitle <- ifelse(what_=="prev", "Prevalences", what_)
      if(is.null(xlim)){
        xlim <- c(min(lo), max(up))
      }
      plot(mu, 1:n_locs, pch = 20,
           xlim = xlim,
           ylim = c(0.5,n_locs+0.5),
           xlab = paste0("Mean and ", 100-alpha*100,"% credible interval"),
           ylab = "Location",
           main = CItitle
      )
      for(l in 1:n_locs){
        graphics::segments(lo[l], l, up[l], l, lwd = 2)
      }
    }
  }
  
}

#' Print method for object of class \code{amis}
#'
#' @param x The output from the function \code{\link{amis}()}.
#' @param ... Other arguments to match the generic \code{print}() function
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

#' Summary method for object of class \code{amis}
#'
#' @param object The output from the function \code{\link{amis}()}.
#' @param ... Other arguments to match the generic \code{summary}() function
#' @return Summary statistics of the fitted model.
#' @export
summary.amis <- function(object, ...) {
  
  if(!inherits(object, "amis")){
    stop("'object' must be of type amis")
  }
  x <- object
  amis_params <- x$amis_params
  
  cat(paste0("=============================================================", "\n"))
  cat(paste0("Data description: \n"))
  cat(paste0("- Number of locations:  ", nrow(x$prevalence_map[[1]]$data),"\n"))
  cat(paste0("- Number of map samples for each location:  ",  ncol(x$prevalence_map[[1]]$data),"\n"))
  cat(paste0("- Number of time points:  ",  length(x$prevalence_map),"\n"))
  
  cat(paste0("=============================================================", "\n"))
  cat(paste0("Fitted model: \n"))
  
  n_sims_total <- nrow(x$simulated_prevalences)
  n_iters <- n_sims_total/amis_params$nsamples
  cat(paste0("- Number of iterations:  ",  n_iters, "\n"))
  cat(paste0("- Total number of simulated samples:  ",  n_sims_total, "\n"))
  cat(paste0("- Target effective sample size:  ",  amis_params$target_ess,"\n"))
  
  # ESS_by_location <- data.frame(ESS = round(x$ess, digits = getOption("digits")))
  # rownames(ESS_by_location) <- colnames(x$weight_matrix)
  # cat(paste0("- Effective sample size by location: \n"))
  # print(ESS_by_location)
  
  which_didnot_exceed_ESS <- which(x$ess < amis_params$target_ess)
  if(length(which_didnot_exceed_ESS)>0){
    message(paste0("- ESS of the following location(s) was lower than the target ESS: "))
    message(paste0(paste(shQuote(colnames(x$weight_matrix)[which_didnot_exceed_ESS]), collapse=", ")))
    ESS_by_location <- data.frame(ESS = round(x$ess[which_didnot_exceed_ESS], digits = getOption("digits")))
    rownames(ESS_by_location) <- colnames(x$weight_matrix)[which_didnot_exceed_ESS]
    print(ESS_by_location)
  }

}


#' Calculate summaries of weighted statistics
#'
#' @param x The output from the function \code{\link{amis}()}.
#' @param what What statistic should be calculated the summaries from. 
#' It must be either 'prev' or the name of one of the model parameters. Default to 'prev'.
#' @param time Time point. Only used if 'what' is set to 'prev'.
#' @param locations Integer vector or location names identifying locations where 
#' summaries should be calculated for. If not specified, summary statistics of 
#' all locations will be provided.
#' @param alpha Numeric value between 0 and 1. Calculations are for the (alpha/2, 1-alpha/2)% quantiles.
#' @return A list with mean, median, and quantiles of the weighted distribution
#' @importFrom  Hmisc wtd.mean
#' @importFrom  Hmisc wtd.quantile
#' @export
calculate_summaries <- function(x, what="prev", time=1, locations=NULL, alpha=0.05) {

  out <- vector(mode='list', length=3)
  names(out) <- c("mean","median","quantiles")
  
  if(is.null(locations)){
    n_locs <- nrow(x$prevalence_map[[1]]$data)
    locations <- 1:n_locs
  }
  
  if(alpha < 0 || alpha > 1){stop("'alpha' must be within 0 and 1.")}
  if(what=="prev"){
    statistic <- x$simulated_prevalences[,time]  #    #     #   
  }else{
    if(!(what%in%colnames(x$param))){
      stop("'parameter' must be either 'prev' or have a valid parameter name.")
    }
    statistic <- x$param[,what]
  }

  wtd <- x$weight_matrix[,locations,drop=F]
  if(x$amis_params$log){wtd <- exp(wtd)}
  
  location_names <- colnames(x$weight_matrix)[locations]
  # weighted mean
  out[[1]] <- sapply(1:length(locations), function(l) Hmisc::wtd.mean(statistic, weights=wtd[,l], normwt = T))
  names(out[[1]]) <- location_names
  # weighted median
  out[[2]] <- sapply(1:length(locations), function(l) Hmisc::wtd.quantile(statistic, weights=wtd[,l], probs=0.5, normwt = T))
  names(out[[2]]) <- location_names
  # weighted quantiles
  out[[3]] <- sapply(1:length(locations), function(l) Hmisc::wtd.quantile(statistic, weights=wtd[,l], probs=c(alpha/2, 1-alpha/2), normwt = T))
  colnames(out[[3]]) <- location_names
  
  return(out)
  
}

