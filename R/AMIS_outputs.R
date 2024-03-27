#' Sample parameters from their weighted distributions given an AMIS output
#'
#' @param x The output from the function \code{\link{amis}}.
#' @param n_samples Number of samples to draw. Default to 200.
#' @param locations Integer identifying the locations. Default to 1.
#' @return Matrix with parameter values and corresponding prevalences for each location.
#' @export
sample_parameters <- function(x, n_samples=200, locations=1) {
  log <- x$amis_params[["log"]]
  n_tims <- ncol(x$simulated_prevalences)
  sampled_pars <- data.frame(matrix(NA, 0, 1 + ncol(x$param) + n_tims))
  weight_matrix <- x$weight_matrix
  param <- x$param
  simulated_prevalences <- x$simulated_prevalences
  for (l in locations) {
    idx <- systematic_sample(n_samples, weight_matrix[,l], log)
    sampled_pars_p <- cbind(l, param[idx, , drop=F], simulated_prevalences[idx, , drop=F])
    sampled_pars <- rbind(sampled_pars, sampled_pars_p)
  }
  colnames(sampled_pars) <- c("location",colnames(x$param),paste0("prev_t",n_tims))
  return(sampled_pars)
}

#' Plot histogram or credible interval of weighted distributions for a given an object of class \code{amis}.
#'
#' @param x The output from the function \code{\link{amis}}.
#' @param what What posterior distribution should be plotted. 
#' It can be 'prev' (default) for plotting prevalences, or one of the parameter names. 
#' @param type Type of plot. It can be 'hist' (default) for histogram, 
#' or 'CI' for credible intervals
#' @param locations Integer identifying the locations. Default to 1.
#' @param time Integer identifying the timepoint. Default to 1.
#' @param measure_central Measure of central tendency for credible interval plots. 
#' It can be 'mean' (default) or 'median'.
#' @param alpha Numeric value between 0 and 1 indicating the endpoints of the 
#' credible intervals, which are evaluated at (alpha/2, 1-alpha/2)% quantiles. 
#' Default (0.05) will create 95% credible intervals.
#' @param breaks Argument passed to \code{\link{wtd.hist}} for histogram plots. 
#' Default to 500.
#' @param cex Argument passed to plots of credible intervals.
#' Default to 1
#' @param lwd Argument passed to plots of credible intervals.
#' Default to 1.
#' @param xlim The x limits of the plots. For for credible intervals of multiple 
#' statistics (i.e. length(what)>1), it must be either NULL or a list with 
#' the x limits for each statistic. Default to NULL.
#' @param main Title for the plot.
#' @param ... Other graphical parameters passed to \code{\link{wtd.hist}}.
#' @importFrom weights wtd.hist
#' @importFrom graphics segments
#' @importFrom graphics par
#' @return A plot.
#' @export
plot.amis <- function(x, what="prev", type="hist", locations=1, time=1, 
                      measure_central="mean", alpha=0.05, 
                      breaks=500, cex=1, lwd=1, xlim=NULL, main=NULL, ...){
  if(!inherits(x, "amis")){
    stop("'x' must be of type 'amis'")
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
    if(length(what)==1){
      xlim <- list(xlim)
    }else{
      if(is.null(xlim)){
        xlim <- replicate(length(what), NULL)
      }
      if(!(is.list(xlim) && length(xlim)==length(what))){
        stop("If length(what)>1, 'xlim' must be either NULL or a list of length(what) elements")
      }
    }
  }

  n_locs <- length(locations)
  location_names <- colnames(x$weight_matrix)[locations]
  
  # Histograms
  if(type=="hist"){
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
            xlim <- range(statistic[weights>0])
          }
        }else{
          xlim <- range(statistic[weights>0])
          xlim[2] <- xlim[2]+0.1*abs(xlim[2])
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
  }
  
  # Credible intervals
  if(type=="CI"){
    if(!measure_central%in%c("mean","median")){stop("Argument 'measure_central' must be either 'mean' or 'median'.")}
    name <- paste0(toupper(substr(measure_central, 1, 1)), substr(measure_central, 2, nchar(measure_central)))
    i <- 1
    for(what_ in what){
      xlim_ <- xlim[[i]]
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
      if(is.null(xlim_)){
        xlim_ <- c(min(lo), max(up))
      }
      plot(mu, 1:n_locs, pch=20, cex=cex,
           xlim = xlim_,
           ylim = c(0.5,n_locs+0.5),
           xlab = paste0(name, " and ", 100-alpha*100,"% credible interval"),
           ylab = "Location",
           main = CItitle
      )
      for(l in 1:n_locs){
        graphics::segments(lo[l], l, up[l], l, lwd = lwd)
      }
      i <- i + 1
    }
  }
  
}

#' Print method for object of class \code{amis}
#'
#' @param x The output from the function \code{\link{amis}}.
#' @param ... Other arguments to match the generic \code{print} function
#' @return Brief description of data and model specifications used to run \code{\link{amis}}.
#' @export
print.amis <- function(x, ...) {
  
  if(!inherits(x, "amis")){
    stop("'x' must be of type 'amis'")
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
  cat(paste0("- Number of new samples drawn within each AMIS iteration:  ",  amis_params$n_samples,"\n"))
  cat(paste0("- Maximum number of iterations:  ",  amis_params$max_iters,"\n"))
  cat(paste0("- Target effective sample size:  ",  amis_params$target_ess,"\n"))
  
}

#' Summary method for object of class \code{amis}
#'
#' @param object The output from the function \code{\link{amis}}.
#' @param ... Other arguments to match the generic \code{summary} function
#' @return Summary statistics of the fitted model.
#' @export
summary.amis <- function(object, ...) {
  
  if(!inherits(object, "amis")){
    stop("'object' must be of type 'amis'")
  }
  x <- object
  amis_params <- x$amis_params
  
  cat(paste0("=============================================================", "\n"))
  cat(paste0("Fitted model: \n"))
  
  n_locs <- nrow(x$prevalence_map[[1]]$data)
  n_sims_total <- nrow(x$simulated_prevalences)
  n_iters <- n_sims_total/amis_params$n_samples
  cat(paste0("- Number of iterations:  ",  n_iters, "\n"))
  cat(paste0("- Total number of simulated samples:  ",  n_sims_total, "\n"))
  cat(paste0("- Target effective sample size:  ",  amis_params$target_ess,"\n"))
  
  ESS_by_location <- data.frame(ESS = round(x$ess, digits = 0))
  if(n_locs<=10){
    cat(paste0("- Effective sample size by location: \n"))
    print(ESS_by_location)
  }else{
    which_didnot_exceed_ESS <- which(x$ess < amis_params$target_ess)
    num_below_ESS <- length(which_didnot_exceed_ESS)
    num_above_ESS <- n_locs - num_below_ESS
    cat(paste0("Number of locations whose ESS exceeded the target ESS:  ",  num_above_ESS, "\n"))
    cat(paste0("Number of locations whose ESS was lower the target ESS:  ",  num_below_ESS, "\n"))
    if(num_below_ESS>0){
      if(num_below_ESS<=10){
        message(paste0("  ESS for the following location(s) was lower than the target ESS: "))
        below_ESS <- data.frame(ESS = round(x$ess[which_didnot_exceed_ESS], digits = 0))
        rownames(below_ESS) <- colnames(x$weight_matrix)[which_didnot_exceed_ESS]
        print(below_ESS)
      }else{
        message("   ESS of more than 10 locations was lower than the target ESS. To see all of them, run 'out$ess[(out$ess < out$amis_params$target_ess)]', where 'out' is an output returned by amis().")
      }
    }
  }
}


#' Calculate summaries of weighted statistics
#'
#' @param x The output from the function \code{\link{amis}}.
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



#' Wrapper function for \code{\link{plot.Mclust}}
#'
#' @param x The output from the function \code{\link{amis}}.
#' @param what A string specifying the type of plot requested:
#' \describe{
#' \item{\code{"uncertainty"}}{A plot of classification uncertainty (default)}
#' \item{\code{"density"}}{A plot of estimated density}
#' \item{\code{"BIC"}}{A plot showing BIC values used to choose the number of components}
#' }
#' @param iteration Integer indicating which iteration the plot should be about. 
#' If NULL (default), the plot will be for the final iteration.
#' See more details in \code{\link{plot.Mclust}}.
#' @param datapoints A string specifying what the datapoints should represent in the 
#' plot of classification uncertainty: 
#' \describe{
#'   \item{\code{"proposed"}}{datapoints will represent the samples simulated from the mixture.
#'   The colours indicate which mixture components the samples were simulated from.}
#'   \item{\code{"fitted"}}{datapoints will show the samples that the mixture model was fitted to, 
#'   i.e. weighted samples from the previous iteration. 
#'   The colour of a datapoint indicates the most likely mixture component the sample belongs to.}
#' }
#' @param main Title of the plot. If NULL, the default title will be displayed. Set to NA for omitting title.
#' @param ... Other arguments to match the \code{plot.Mclust} function
#' @importFrom graphics title
#' @return A plot for model-based clustering results.
#' @export
plot_mixture_components <- function(x, what="uncertainty", iteration=NULL, 
                                    datapoints="proposed", main=NULL, xlim=NULL, ylim=NULL, ...) {
  if(!inherits(x, "amis")){
    stop("'x' must be of type 'amis'")
  }
  if(!(what%in%c("uncertainty", "density", "BIC"))){
    stop("'what' must be 'uncertainty', 'density', or 'BIC'.")
  }
  if(what=="uncertainty"){
    if(!(datapoints%in%c("proposed", "fitted"))){
      stop("'datapoints' must be either 'proposed' or 'fitted'.")
    }
  }
  total_n_samples <- nrow(x$simulated_prevalences)
  last_iteration <- total_n_samples/x$amis_params$n_samples
  if(last_iteration==1){
    stop("Algorithm stoped after one iteration. Therefore, no mixture model has been fitted.")
  }
  if(is.null(iteration)){
    iteration <- last_iteration
    if(iteration==1){
      stop("Algorithm stoped after one iteration. Therefore, no mixture model has been fitted.")
    }
  }else{
    if(iteration==1){
      stop("At iteration 1, no mixture model has been fitted.")
    }
    if(!(iteration%in%c(2:last_iteration))){
      stop("'iteration' must be an integer that does not exceed the maximum number of iterations that have been run.")
    }
  }
  clustering <- x$components_per_iteration[[iteration]]
  
  param_names <- colnames(x$param)
  d <- length(param_names)
  
  if(what=="uncertainty"){
    if(datapoints=="proposed"){
      cat(paste0("Plotting components of the mixture model (and samples simulated from) at iteration ", iteration, "...\n"))
      clustering$data <- clustering$data_proposed
      n <- length(clustering$compon_proposal)
      clustering$n <- n
      comp_idx <- clustering$compon_proposal
      G <- clustering$G
      # forcing same point sizes
      z <- matrix(0, n, G)
      for(i in 1:n){
        z[i,comp_idx[i]] <- 1
      }
      # not displaying samples that will be given zero weight
      # ... because of zero prior density:
      wh_delete <- NULL
      n_samples <- x$amis_params$n_samples
      idx <- 1:n_samples + n_samples*(iteration-1)
      prior_densities <- x$prior_density[idx]
      if(x$amis_params$log){
        wh_delete <- which(prior_densities==-Inf)
      }else{
        wh_delete <- which(prior_densities==0)
      }
      # ... because they are out of parameter boundaries:
      if(!is.null(x$amis_params$boundaries_param)){
        nrows <- nrow(clustering$data)
        ncols <- ncol(clustering$data)
        for(j in 1:ncols){
          boundaries <- x$amis_params$boundaries_param[j,]
          for(i in 1:nrows){
            if(clustering$data[i,j]<boundaries[1] || clustering$data[i,j]>boundaries[2]){
              wh_delete <- c(wh_delete, i)
            }
          }
        }
      }
      wh_delete <- unique(wh_delete)
      if(length(wh_delete)>0){
        clustering$data <- clustering$data[-wh_delete,]
        z <- z[-wh_delete,]
        n <- nrow(z)
        clustering$n <- n  
      }
      clustering$z <- z
    }else if(datapoints=="fitted"){
      # forcing same point sizes
      z <- clustering$z
      ncols <- ncol(z)
      for(i in 1:nrow(z)){
        which_comp <- which.max(z[i,])
        rest <- (1:ncols)[-which_comp]
        z[i,which_comp] <- 1
        z[i,rest] <- 0
      }
      clustering$z <- z
      cat(paste0("Plotting components of the fitted mixture model at iteration ", iteration, " and weighted samples of the previous iteration...\n"))
    }
    colnames(clustering$data) <- colnames(x$param)
    if(d==1){
      mclust::plot.Mclust(clustering, what = what, xlim=xlim, ylim=ylim, xlab=colnames(x$param), ...)
    }else{
      mclust::plot.Mclust(clustering, what = what, xlim=xlim, ylim=ylim, ...)
    }
    if(datapoints=="proposed"){
      default_main <- "Proposed Parameter Values"
    }else if(datapoints=="fitted"){
      default_main <- "Parameter Values that the Mixture Model Was Fitted To"
    }
    main <- ifelse(is.null(main), default_main, main)
    title(main = main)
  }else if(what=="density"){
    # xlab = colnames(x$param)[1]
    # ylab = colnames(x$param)[2]
    cat(paste0("Plotting density of the mixture model at iteration ", iteration, "...\n"))
    # mclust::plot.Mclust(clustering, what=what, xlab=xlab, ylab=ylab, ...)
    colnames(clustering$data) <- colnames(x$param)
    if(d==1){
      mclust::plot.Mclust(clustering, what=what, xlab=param_names, xlim=xlim, ylim=ylim, ...)
    }else{
      mclust::plot.Mclust(clustering, what=what, xlim=xlim, ylim=ylim, ...)
    }
    main <- ifelse(is.null(main), "Proposal Density", main)
    title(main = main)
  }else if(what=="BIC"){
    cat(paste0("Plotting BIC against number of components at iteration ", iteration, "...\n"))
    mclust::plot.Mclust(clustering, what = what, 
                        legendArgs = list(plot=FALSE), xlab="Number of Components", ...)
    main <- ifelse(is.null(main), "BIC", main)
    title(main = main)
  }
  
}


