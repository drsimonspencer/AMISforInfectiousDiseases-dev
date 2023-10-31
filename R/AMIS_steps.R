#' @useDynLib AMISforInfectiousDiseases
#' @importFrom Rcpp sourceCpp
NULL

#' Environment where intermittent outputs are saved
#' @export
amis_env <- new.env()

#' Produce list containing the default AMIS parameters
#' 
#' For description of AMIS parameters, see argument \code{amis_params} in \code{\link{amis}}.
#' @return List containing the default AMIS parameters.
#' @export
default_amis_params <- function() {
  amis_params<-list(nsamples=500, boundaries=c(0,1), boundaries_param=NULL, 
                    bayesian=FALSE,
                    mixture_samples=1000, df=3,
                    target_ess=500, log=FALSE, max_iters=12,
                    intermittent_output=FALSE, 
                    delta=0.01, sigma=NULL, breaks=NULL)
  return(amis_params)
}

#' Check inputs of \code{amis} function
#' 
#' Check whether all the inputs of \code{\link{amis}} function are as expected.
#' @inheritParams amis
#' @export
check_inputs <- function(prevalence_map, transmission_model, prior, amis_params, seed) {

  if (!(is.matrix(prevalence_map) || is.data.frame(prevalence_map) || is.list(prevalence_map))) {
    stop(("prevalence_map must be either \n - a matrix or data frame of size #locations by #samples (for one timepoint); or \n - a list with n_tims timepoints, each one with a matrix named 'data'."))
  }
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  dims <- lapply(prevalence_map, dim)
  if(!all(sapply(dims, FUN = identical, dims[[1]]))){
    stop("'prevalence_map' must have the same dimension (number of spatial units and number of samples) at each time point. If data for some locations are missing at a timepoint, set to NA.")
  }
  stopifnot("'transmission_model' must be a function." = is.function(transmission_model))
  stopifnot("'prior' must be a list." = is.list(prior))
  stopifnot("'prior' must be a list of two elements." = (length(prior)==2))
  stopifnot("'prior'  must be a list with two elements called 'dprior' and 'rprior'." = all(sort(names(prior))==c("dprior","rprior")))
  stopifnot("'prior$dprior' must be a function." = is.function(prior$dprior))
  stopifnot("'prior$rprior' must be a function." = is.function(prior$rprior))
  stopifnot("'log' must be a single logical value" = (length(amis_params$log)==1 && is.logical(amis_params$log)))
  delta <- amis_params$delta
  sigma <- amis_params$sigma
  nsamples <- amis_params$nsamples
  mixture_samples <- amis_params$mixture_samples
  df <- amis_params$df
  target_ess <- amis_params$target_ess
  max_iters <- amis_params$max_iters
  bayesian <- amis_params$bayesian
  breaks <- amis_params$breaks
  boundaries <- amis_params$boundaries
  boundaries_param <- amis_params$boundaries_param
  boundaries <- as.numeric(boundaries)
  if(length(boundaries)!=2){stop("'boundaries' must be a vector of length 2.")}
  if(!(diff(boundaries)>0)){stop("The second element of 'boundaries' must be larger than the first one.")}
  if(!is.null(breaks)){
    if(any(breaks!=sort(breaks))){
      stop("'breaks' should be a vector with increasing values.")
    }
    if(length(breaks)!=length(unique(breaks))){
      stop("'breaks' must not have repeated values.")
    }
    if(is.finite(boundaries[1])){
      if(breaks[1]!=boundaries[1]){
        stop("The first entry of 'breaks' must be equal to the left boundary if the left boundary is finite.")
      }
    }
    if(is.finite(boundaries[2])){
      if(breaks[length(breaks)]!=boundaries[2]){
        stop("The last entry of 'breaks' must be equal to the right boundary if the right boundary is finite.")
      }
    }
  }
  if(!is.null(boundaries_param)){
    stopifnot("'boundaries_param' must be a (#parameters x 2) matrix" = (is.matrix(boundaries_param)) && 
                (ncol(boundaries_param)==2) && (nrow(boundaries_param)==ncol(rprior(1))))
  }
  stopifnot("'delta' must be either NULL or a single positive numeric value" = ((length(delta)==1 && is.numeric(delta) && delta>0) || is.null(delta)))
  stopifnot("'sigma' must be either NULL or a single positive numeric value" = ((length(sigma)==1 && is.numeric(sigma) && sigma>0) || is.null(sigma)))
  stopifnot("'nsamples' must be a single numeric value greater than 1" = (length(nsamples)==1 && is.numeric(nsamples) && nsamples>1))
  stopifnot("'mixture_samples' must be a single numeric value" = (length(mixture_samples)==1 && is.numeric(mixture_samples)))
  stopifnot("'df' must be a single numeric value greater than 0" = (length(df)==1 && is.numeric(df) && df>0))
  stopifnot("'target_ess' must be a single numeric value greater than 0" = (length(target_ess)==1 && is.numeric(target_ess) && target_ess>0))
  stopifnot("'max_iters' must be a single numeric value greater than 1" = (length(max_iters)==1 && is.numeric(max_iters) && max_iters>1))
  stopifnot("'seed' must be either NULL or a single numeric value" = ((length(seed)==1 && is.numeric(seed)) || is.null(seed)))
  if(is.null(prevalence_map[[1]]$likelihood) && is.null(c(amis_params[["delta"]], amis_params[["sigma"]], amis_params[["breaks"]]))){
    stop("At least one of the inputs ('delta','sigma','breaks') must not be NULL if a likelihood function is not provided.")
  }
  if(!amis_params[["bayesian"]] && is.null(c(amis_params[["delta"]], amis_params[["sigma"]], amis_params[["breaks"]]))){
    stop("At least one of the inputs ('delta','sigma','breaks') must not be NULL if 'bayesian' is set to FALSE.")
  }
  if(!is.null(prevalence_map[[1]]$likelihood)){
    cat("Likelihood function was provided and will be used to calculate likelihood terms. \n")
  }
  if(is.null(prevalence_map[[1]]$likelihood) && amis_params[["bayesian"]]){
    if(!is.null(amis_params[["breaks"]])){
      cat("Histogram method will be used in the estimation of the likelihood as 'breaks' was provided. \n")
    }else{
      if(!is.null(amis_params[["sigma"]])){
        cat("Gaussian kernel will be used in the estimation of the likelihood as 'sigma' was provided. \n")
      }else{
        cat("Uniform kernel will be used in the estimation of the likelihood. \n")
      }
    }
  }
  if(amis_params[["bayesian"]]){
    cat("Bayesian update of the weights will be implemented as 'bayesian' was set to TRUE. \n")
  }else{
    if(is.null(prevalence_map[[1]]$likelihood)){
      if(!is.null(amis_params[["breaks"]])){
        cat("Histogram method will be used in the empirical Radon-Nikodym derivative as 'breaks' was provided. \n")
      }else{
        if(!is.null(amis_params[["sigma"]])){
          cat("Gaussian kernel will be used in the empirical Radon-Nikodym derivative as 'sigma' was provided. \n")
        }else{
          cat("Uniform kernel will be used in the empirical Radon-Nikodym derivative. \n")
        }
      }
    }else{
      if(!is.null(amis_params[["breaks"]])){
        cat("Histogram method will be used in the denominator of the empirical Radon-Nikodym derivative as 'breaks' was provided. \n")
      }else{
        if(!is.null(amis_params[["sigma"]])){
          cat("Gaussian kernel will be used in the denominator of the empirical Radon-Nikodym derivative as 'sigma' was provided. \n")
        }else{
          cat("Uniform kernel will be used in the denominator of the empirical Radon-Nikodym derivative. \n")
        }
      }
    }
  }
  
  n_tims <- length(prevalence_map)
  locs_no_data <- NULL
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  for(l in 1:n_locs){
    num_valid_datapts <- 0L
    for(t in 1:n_tims){
      data_l_t <- prevalence_map[[t]]$data[l,]
      num_valid_datapts <- num_valid_datapts + sum(is.finite(data_l_t) & (data_l_t>=boundaries[1]) & (data_l_t<=boundaries[2]))
    }
    if(num_valid_datapts==0L){
      locs_no_data <- c(locs_no_data, l)
    }
  }
  if(!is.null(locs_no_data)){
    if(length(locs_no_data)==1L){
      message(paste0("Location ", shQuote(locs_no_data), " has no valid map samples. Its corresponding posterior samples will therefore be driven by the prior."))
    }else{
      message(paste0("Locations (", paste(shQuote(locs_no_data), collapse=","), ") have no valid map samples. Their corresponding posterior samples will therefore be driven by the prior."))
    }
  }
 
  return(list(locs_no_data=locs_no_data)) 
}


#' Compute likelihood for each additional simulation across timepoints
#'
#' Calls evaluate likelihood for each timepoint.
#' @param param A matrix containing the sampled parameter vectors.
#' @param prevalence_map A list with one entry for each timepoint.
#' Each entry must be a list containing objects \code{data} (an L x M matrix of data);
#' and optional function \code{likelihood} taking arguments \code{param} (model parameters used in the simulation), \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.
#' If a likelihood is not specified then it is assumed that
#' the data consist of samples from a geo-statistical model and empirical methods are used.
#' @param simulated_prevalences An n x timepoints matrix of prevalences simulated from the transmission model.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations (optional).
#' @param which_valid_sim_prev_iter List of T elements, where each one indicates which simulated prevalences are valid at the current iteration
#' @param which_valid_prev_map List showing which prevalence map samples are valid
#' @param log_norm_const_gaussian Normalising constant (in log scale) for the Gaussian kernel. It is only used if Gaussian kernel is used.
#' @return A larger array with the likelihoods of the new simulations joined to the existing array \code{likelihoods}.
compute_likelihood <- function(param,prevalence_map,simulated_prevalences,amis_params,likelihoods=NULL,
                               which_valid_sim_prev_iter,which_valid_prev_map,log_norm_const_gaussian) {
  n_tims <- length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  n_sims <- dim(simulated_prevalences)[1]
  lik <- array(NA, c(n_tims,n_locs,n_sims)) # this way around to avoid abind -- different to elsewhere
  for (t in 1:n_tims) {
    lik[t,,] <- evaluate_likelihood(param,prevalence_map[[t]],simulated_prevalences[,t],amis_params,
                                    which_valid_sim_prev_iter[[t]],which_valid_prev_map[[t]],log_norm_const_gaussian[t,,]) 
  }
  if (!is.null(likelihoods)) {
    lik <- array(c(likelihoods,lik), c(n_tims,n_locs,dim(likelihoods)[3]+n_sims))
  }
  return(lik)
}

#' Evaluate likelihood for each additional simulation for a single timepoint
#'
#' Implements analytical likelihoods if a likelihood function is available; otherwise histogram or empirical
#' likelihoods are generated based on samples from a geostatistical map.
#' @param param A matrix containing the sampled parameter vectors.
#' @param prevalence_map A list containing objects \code{data} (an L x M matrix of data);
#' and \code{likelihood} a function taking arguments \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.    
#' @param prev_sim A vector of simulated prevalences
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param which_valid_sim_prev_iter Vector showing which simulated prevalences are valid at the current iteration at time t.
#' @param which_valid_prev_map_t List showing which prevalence map samples are valid at time t.
#' @param log_norm_const_gaussian_t Normalising constant (in log scale) for the Gaussian kernel at time t. It is only used if Gaussian kernel is used.
#' @return A locations x simulations matrix of (log-)likelihoods.
evaluate_likelihood <- function(param,prevalence_map,prev_sim,amis_params, 
                                which_valid_sim_prev_iter,which_valid_prev_map_t,log_norm_const_gaussian_t) {

  logar <- amis_params[["log"]]
  # f <- matrix(NA,dim(prevalence_map$data)[1],length(prev_sim))
  
  if (!is.null(prevalence_map$likelihood)) {
    
    likelihood_fun <- prevalence_map$likelihood
    
    
    # Here, the user-defined 'likelihood_fun' returns an (M_l x nsims) matrix 
    # for a particular location l
    f <- f_user_defined_l(likelihood_fun,
                          param,
                          prevalence_map=prevalence_map$data, 
                          prev_sim=prev_sim, 
                          which_valid_sim_prev_iter=which_valid_sim_prev_iter, 
                          which_valid_prev_map_t=which_valid_prev_map_t, 
                          logar=logar)

    # # Here, the user-defined 'likelihood_fun' returns an M_l-length vector
    # # for a particular location l and a single simulated prevalence
    # f <- f_user_defined_l_r(likelihood_fun,
    #                         param,
    #                         prevalence_map=prevalence_map$data, 
    #                         prev_sim=prev_sim, 
    #                         which_valid_sim_prev_iter=which_valid_sim_prev_iter, 
    #                         which_valid_prev_map_t=which_valid_prev_map_t, 
    #                         logar=logar)
    # 
    # # Here, the user-defined 'likelihood_fun' returns a single point:
    # # the likelihood of observing a simulated value r given a sample m at a particular location l
    # # This approach can be slow because the R function has to be called from Rcpp too many times
    # f <- f_user_defined_l_m_r(likelihood_fun,
    #                           param,
    #                           prevalence_map=prevalence_map$data,
    #                           prev_sim=prev_sim,
    #                           which_valid_sim_prev_iter=which_valid_sim_prev_iter,
    #                           which_valid_prev_map_t=which_valid_prev_map_t,
    #                           logar=logar)

    # # R code of previous version of the package
    # locs<-which(!is.na(prevalence_map$data[,1])) # if there is no data for a location, do not update weights.
    # f[locs,]<-t(prevalence_map$likelihood(param,prevalence_map$data[locs,,drop=FALSE],prev_sim,amis_params[["log"]])) # likelihood function must be vectorised.
  } else {
    boundaries <- amis_params[["boundaries"]]
    if (is.null(amis_params[["breaks"]])) {
      if(!is.null(amis_params[["sigma"]])){
        sd <- amis_params[["sigma"]]
        f <- f_estimator_Gaussian(prevalence_map=prevalence_map$data, 
                                  prev_sim=prev_sim, 
                                  sd=sd, 
                                  which_valid_sim_prev_iter=which_valid_sim_prev_iter, 
                                  which_valid_prev_map_t=which_valid_prev_map_t,
                                  log_norm_const_gaussian_t=log_norm_const_gaussian_t, 
                                  logar=logar)
      }else{
        delta <- amis_params[["delta"]]
        f <- f_estimator_uniform(prevalence_map=prevalence_map$data, 
                                 prev_sim=prev_sim, 
                                 delta=delta,
                                 which_valid_sim_prev_iter=which_valid_sim_prev_iter, 
                                 which_valid_prev_map_t=which_valid_prev_map_t,
                                 boundaries=boundaries, 
                                 logar=logar)
        # # R code of previous version of the package
        # for (i in 1:length(prev_sim)) {
        #   # f[,i]<-rowSums(abs(prevalence_map$data[locs,,drop=FALSE]-prev_sim[i])<=delta/2)/delta
        #   f[,i]<-rowSums(abs(prevalence_map$data[locs,,drop=FALSE]-prev_sim[i])<=delta/2)/(delta*M)
        # }
      }
    } else {
      breaks <- amis_params[["breaks"]]
      f <- f_estimator_histogram(prevalence_map=prevalence_map$data, 
                                 prev_sim=prev_sim, 
                                 breaks=breaks, 
                                 which_valid_prev_map_t=which_valid_prev_map_t,
                                 logar=logar)
      # # R code of previous version of the package                                      
      # L<-length(breaks)
      # lwr<-breaks[1:(L-1)]
      # upr<-breaks[2:L]
      # wdt<-upr-lwr
      # for (l in 1:L) {
      #   wh<-which(prev_sim>=lwr[l] & prev_sim<upr[l])
      #   if (length(wh)>0) {
      #     f[locs,wh]<-rowSums(prevalence_map$data[locs,,drop=FALSE]>=lwr[l] & prevalence_map$data[locs,,drop=FALSE]<upr[l])/(ncol(prevalence_map$data)*wdt[l])
      #   }
      # }
    }
    # if (amis_params[["log"]]) {f <- log(f)}   # all functions already return f in the requested scale
  }
  return(f)
}

#' Compute weight matrix across timepoints using appropriate method
#' 
#' Wrapper function to select appropriate method to calculate weight matrix.
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations.
#' @param simulated_prevalence An n_sims x n_tims matrix containing the simulated prevalence values for each of the
#'     parameter samples. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as the rows in \code{simulated_prevalence}.
#' @param locs_empirical List indicating, at each time, which locations are updated without using Bayesian method.
#' @param locs_bayesian List indicating, at each time, which locations are updated using Bayesian method.
#' @param bool_valid_sim_prev Matrix of n_tims columns, where column is a logical vector indicating which simulated prevalences are valid.
#' @param which_valid_sim_prev List indicating, at each time, which simulated prevalences are valid.
#' @param which_invalid_sim_prev List indicating, at each time, which simulated prevalences are invalid
#' @param which_valid_locs_prev_map List showing which locations have valid data at each time
#' @param locations_with_no_data Vector indicating which locations have no data at any time point
#' @return normalised weight matrix.
compute_weight_matrix <- function(likelihoods, simulated_prevalence, amis_params, first_weight, 
                                  locs_empirical, locs_bayesian,
                                  bool_valid_sim_prev, which_valid_sim_prev, which_invalid_sim_prev, 
                                  which_valid_locs_prev_map, locations_with_no_data) {

  n_tims <- dim(likelihoods)[1]
  n_locs <- dim(likelihoods)[2]
  n_sims <- dim(likelihoods)[3]
  weight_matrix <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs)

  for (t in 1:n_tims) {

    lik_mat <- t(array(likelihoods[t,,], dim=c(n_locs, n_sims)))
    
    # Update the weights by the latest likelihood (filtering)
    if (!amis_params[["bayesian"]]){
      
      # If this is the first timepoint where there is data for a location then do not use bayesian update
      # locs_empirical = which(locations_first_t == t)
      # locs_bayesian = which(locations_first_t < t)

      if(!is.null(locs_bayesian[[t]])){
        weight_matrix <- compute_weight_matrix_bayesian_Rcpp(lik_mat, amis_params, weight_matrix,
                                                             which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]],
                                                             locs_bayesian[[t]])
      }

      if (is.null(amis_params[["breaks"]])){
          if(is.null(amis_params[["sigma"]])){
            if(!is.null(locs_empirical[[t]])){
              weight_matrix <- compute_weight_matrix_empirical_uniform(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix,
                                                                       bool_valid_sim_prev[,t], 
                                                                       which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                       locs_empirical[[t]])
            }
          }else{
            if(!is.null(locs_empirical[[t]])){
              weight_matrix <- compute_weight_matrix_empirical_gauss(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix, 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_empirical[[t]])
            }
          }
      } else {
        if(!is.null(locs_empirical[[t]])){
          weight_matrix <- compute_weight_matrix_empirical_histogram(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_empirical[[t]])
        }
      }
    } else {
      if(!is.null(which_valid_locs_prev_map[[t]])){
        weight_matrix <- compute_weight_matrix_bayesian_Rcpp(lik_mat, amis_params, weight_matrix,
                                                          which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                          which_valid_locs_prev_map[[t]])
      }
    }
    
    if(length(locations_with_no_data)>0 && length(which_invalid_sim_prev[[t]])>0){
      weight_inval_prev <- ifelse(amis_params[["log"]], -Inf, 0)
      weight_matrix[which_invalid_sim_prev[[t]]+1L, locations_with_no_data] <- weight_inval_prev
    }
    
  }
  
  # renormalise weights
  if (amis_params[["log"]]) {
    M<-apply(weight_matrix,2,max)
    wh<-which(M>-Inf)
    weight_matrix[,wh]<-weight_matrix[,wh]-rep(M[wh]+log(colSums(exp(weight_matrix[,wh,drop=FALSE]-rep(M[wh],each=n_sims)))),each=n_sims)
  } else {
    S<-colSums(weight_matrix)
    wh<-which(S>0)
    weight_matrix[,wh]<-weight_matrix[,wh]/rep(S[wh],each=n_sims)
  }
  return(weight_matrix)
}

#' Compute weight matrix using empirical Radon-Nikodym derivative using Uniform kernel (old code)
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based on the empirical Radon-Nikodym derivative, taking into account
#' geostatistical prevalence data for the specific location and the prevalence values
#' computed from the transmission model for the specific parameter sample.
#'
#' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
#' NB: transpose of slice of array.
#' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#' @param weight_matrix An n_sims x n_locs matrix containing the current values of the weights.
#' @noRd
#' @return An updated weight matrix.
compute_weight_matrix_empirical <- function(likelihoods, prev_sim, amis_params, weight_matrix) {
  delta<-amis_params[["delta"]]
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  new_weights<-weight_matrix
  for (i in 1:length(prev_sim)) {
    wh<-which(abs(prev_sim-prev_sim[i])<=delta/2)
    g_terms<-weight_matrix[wh,locs,drop=FALSE]
    if (amis_params[["log"]]) {
      M<-apply(g_terms,2,max)
      non_zero_locs<-locs[which(M>-Inf)]
      M<-M[which(M>-Inf)]
      g_terms = g_terms[,which(M>-Inf),drop=FALSE]
      new_weights[i,non_zero_locs]<-weight_matrix[i,non_zero_locs]+likelihoods[i,non_zero_locs]-M-log(colSums(exp(g_terms-rep(M,each=length(wh)))))+log(delta)
    } else {
      g<-colSums(g_terms)
      non_zero_locs<-locs[which(g>0)]
      g<-g[which(g>0)]
      new_weights[i,non_zero_locs]<-weight_matrix[i,non_zero_locs]*likelihoods[i,non_zero_locs]/g*delta
    }
    zero_locs<-setdiff(locs,non_zero_locs)
    if (length(zero_locs)>0) {new_weights[i,zero_locs]<-ifelse(amis_params[["log"]],-Inf,0)}
  }
  return(new_weights)
}

#' Compute weight matrix using empirical Radon-Nikodym derivative (with fixed breaks) (old code)
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based on the empirical Radon-Nikodym derivative, taking into account
#' geostatistical prevalence data for the specific location and the prevalence value
#' computed from the transmission model for the specific parameter sample.
#'
#' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
#' NB: transpose of slice of array.
#' @param prev_sim A vector containing the simulated prevalence value for each
#'     parameter sample. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#' @param weight_matrix A matrix containing the current values of the weights.
#' @noRd
#' @return An updated weight matrix.
compute_weight_matrix_histogram<-function(likelihoods, prev_sim, amis_params, weight_matrix) {
  breaks<-amis_params[["breaks"]] # NB top entry in breaks must be strictly larger than the largest possible prevalence.
  L<-length(breaks)
  lwr<-breaks[1:(L-1)]
  upr<-breaks[2:L]
  wdt<-upr-lwr
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  new_weights<-weight_matrix
  for (l in 1:L) {
    wh<-which(prev_sim>=lwr[l] & prev_sim<upr[l])
    if (length(wh)>0) {
      g_terms<-weight_matrix[wh,locs,drop=FALSE]
      if (amis_params[["log"]]) {
        M<-apply(g_terms,2,max)
        non_zero_locs<-locs[which(M>-Inf)]
        M<-M[which(M>-Inf)]
        g_terms = g_terms[,which(M>-Inf),drop=FALSE]   # Revise this
        new_weights[wh,non_zero_locs]<-weight_matrix[wh,non_zero_locs]+likelihoods[wh,non_zero_locs]-rep(M+log(colSums(exp(g_terms-M))),each=length(wh))+log(wdt[l])
      } else {
        g<-colSums(g_terms)
        non_zero_locs<-locs[which(g>0)]
        g<-g[which(g>0)]
        new_weights[wh,non_zero_locs]<-weight_matrix[wh,non_zero_locs]*likelihoods[wh,non_zero_locs]/rep(g,each=length(wh))*wdt[l]
      }
      zero_locs<-setdiff(locs,non_zero_locs)
      if (length(zero_locs)>0) {
        cat("Non-zero IUs",zero_locs,"\n")
        new_weights[wh,zero_locs]<-ifelse(amis_params[["log"]],-Inf,0)
      }
    }
  }
  wh<-which(prev_sim<min(lwr) | prev_sim>=max(upr))
  new_weights[wh,]<-ifelse(amis_params[["log"]],-Inf,0)
  return(new_weights)
}

#' Compute weight matrix without using empirical Radon-Nikodym derivative (old code)
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based only on the previous weight (prior) times the likelihood.
#'
#' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
#' NB: transpose of slice of array.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#' @param weight_matrix A matrix containing the current values of the weights.
#' @noRd
#' @return An updated weight matrix.
compute_weight_matrix_bayesian <- function(likelihoods, amis_params, weight_matrix) {
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  if (amis_params[["log"]]) {
    weight_matrix[,locs]<-weight_matrix[,locs]+likelihoods[,locs]
  } else {
    weight_matrix[,locs]<-weight_matrix[,locs]*likelihoods[,locs]
  }
  return(weight_matrix)
}

#' Compute the current effective sample size
#'
#' This function returns the effective sample size (ESS) for each location. 
#'  For each column of the weight matrix \code{weight_mat}, the ESS component is computed as
#' \deqn{(\sum_{i=1}^{N} w_{i}^2 )^{-1}}
#' where \eqn{N} is the number of sampled parameter values for each location.
#'
#' @param weight_mat The weight matrix. A N x L matrix with L the number of locations
#'     and N the number of sampled parameter values.
#' @param log logical indicating if the weights are on the log scale. 
#' @return A vector containing the ESS value for each location.
#'
#' @seealso \code{\link{compute_weight_matrix}}
calculate_ess <- function(weight_mat,log) {
  ess<-rep(0,dim(weight_mat)[2])
  if (log) {
    M<-apply(weight_mat,2,max)
    wh<-which(M>-Inf)
    M<-M[wh]
    ess[wh]<-exp(-2*M)*colSums(exp(2*(weight_mat[,wh,drop=FALSE]-rep(M,each=dim(weight_mat)[1]))))^(-1)
  } else {
    S<-colSums(weight_mat^2)
    wh<-which(S>0)
    ess[wh]<-1/S[wh]
  }
  return(ess)
}

#' Calculate sum of weight matrix for active locations
#'
#' This function sums the rows of the weight matrix \code{weight_matrix} for which
#' the effective sample size ESS is below a target size \code{target_size}.
#'
#' @param weight_matrix The weight_matrix as returned by
#'     \link{compute_weight_matrix}
#' @param ess The effective sample size vector as returned by
#'     \link{calculate_ess}
#' @param target_size A number representing the target size for the sample.
#' @param log A logical indicating if the weights are logged.
#' @return Vector containing the row sums of the active columns of the weight matrix.
update_according_to_ess_value <- function(weight_matrix, ess, target_size,log) {
  active_cols <- which(ess < target_size)
  if (log) {
    M<-apply(weight_matrix[,active_cols,drop=FALSE],1,max)
    wh<-which(M==-Inf)
    M[wh]<-0
    return(M+log(rowSums(exp(weight_matrix[,active_cols,drop=FALSE]-M))))
  } else {
    return(rowSums(weight_matrix[,active_cols,drop=FALSE]))
  }
}
#' Systematic resampling function
#' 
#' Implement systematic resampling to reduce variance in weighted particle selection
#' @param nsamples number of samples to draw
#' @param weights vector of length equal to the number of particles, containing their weights
#' @param log logical indicating if weights are log-weights
#' @return vector of indices of the sampled particles 
systematic_sample <- function(nsamples,weights,log=F) {
  if (log) {
    M<-max(weights)
    log_sum_weights<-M+log(sum(exp(weights-M)))
    cum <- cumsum(exp(weights-log_sum_weights))
  } else {
    cum <- cumsum(weights)/sum(weights) # cumulative sum of normalised weights
  }
  u <- stats::runif(1)/nsamples+0:(nsamples-1)/nsamples
  return(1+matrix(rep(u,length(weights))>rep(cum,each=nsamples),nsamples,length(weights))%*%matrix(1,length(weights),1))
}
#' Fit mixture to weighted sample
#' 
#' Weights are implemented by using systematic resampling to obtain an unweighted set of parameters.
#' An unweighted mixture is then fitted using \code{fit_mixture}.
#' @param parameters An N x d matrix containing the sampled values for the d parameters.
#' @param nsamples The number of parameter to resample as data to fit the mixture to.
#' @param weights A vector of weights with length N.
#' @param log logical indicating if weights are logged.
#' @return A list of the mixture components (see function \code{\link{fit_mixture}})
#'     \describe{
#'       \item{\code{probs}}{The mixture weights}
#'       \item{\code{Mean}}{The means of the components}
#'       \item{\code{Sigma}}{The covariance matrices of the components}
#'       \item{\code{G}}{Number of components}
#'       \item{\code{BIC}}{BIC of fitted mixture}
#'       \item{\code{ModelName}}{Model name from package mclust}
#'     }
#'
#' @seealso \code{\link{fit_mixture}}
weighted_mixture <- function(parameters, nsamples, weights, log=F) {
  sampled_idx <- systematic_sample(nsamples,weights,log)
  if (length(unique(sampled_idx))==1) {warning("Only one particle with sufficient weight. Will result in a non-invertible covariance matrix for the mixture. \n")}
  return(fit_mixture(parameters[sampled_idx,,drop=FALSE]))
}
#' Sample new parameters
#'
#' This function generates \code{nsamples} new model parameter values according the
#' t-distribution with \code{df} degrees of freedom and the mixture components \code{mixture}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{weighted_mixture}}
#' @param nsamples A number of new parameters to sample (integer)
#' @param df The degrees of freedom for the t-distributed proposal distribution.
#' @param prior list containing the functions \code{rprior} and \code{dprior}
#' @param log A logical indicating if densities 
#' @return A list containing \code{params}, an \code{nsamples} x d matrix containing the sampled parameter values and
#' \code{prior_density}, the corresponding vector of prior densities.
#'
#' @seealso \code{\link{fit_mixture}}
sample_new_parameters <- function(mixture, nsamples, df, prior, log) {
  prior_density<-rep(NA,nsamples) 
  i<-1
  while (i <= nsamples) {
    compo <- sample(1:mixture$G, 1, prob = mixture$probs)
    proposal <- mnormt::rmt(1,mean=mixture$Mean[,compo],S=mixture$Sigma[,,compo],df=df)
    density_of_proposal <- prior$dprior(proposal,log=log)
    if (all(!is.na(proposal)) && ((log && density_of_proposal>-Inf) || (!log && density_of_proposal>0))) {
      if (i==1) {params<-matrix(NA,nsamples,length(proposal))}
      params[i,]<-proposal
      prior_density[i]<-density_of_proposal
      i<-i+1
    }
  }
  return(list(params=params,prior_density=prior_density))
}

#' Update the components of the mixture
#'
#' This function updates the mixture \code{components} according to
#' the current mixture \code{mixture} generated at iteration \code{iter}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{fit_mixture}}
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{G}}{A numeric vector containing the number of components from each AMIS iteration}
#'       \item{\code{Sigma}}{A list of covariance matrices for each component}
#'       \item{\code{Mean}}{A list of means for each component}
#'       \item{\code{probs}}{A list probability weights for each component}
#'     }
#' @param iter The current iteration index (integer)
#' @return The updated \code{components} list
#'
#' @seealso \code{\link{weighted_mixture}}, \code{\link{fit_mixture}}
update_mixture_components <- function(mixture, components, iter) {
  components$G[iter] <- mixture$G
  G_previous <- sum(components$G[1:(iter - 1)]) # Number of pre-existing components
  for (i in 1:mixture$G) {
    components$Sigma[[i + G_previous]] <- mixture$Sigma[, , i]
    components$Mean[[i + G_previous]] <- mixture$Mean[,i]
    components$probs[[i + G_previous]] <- mixture$probs[i] ### scale by number of points if nsamples varies by iteration
  }
  return(components)
}

#' Compute the prior/proposal ratio
#'
#' This function returns the ratio between the prior and proposal distribution
#' for each sampled parameter value (i.e. each row in \code{param}).
#' This function returns the first weight
#' See step (4) of the AMIS algorithm in
#' Integrating geostatistical maps and infectious disease transmission models 
#' using adaptive multiple importance sampling.
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' Ann. Appl. Stat. 15 (4) 1980 - 1998, December 2021.
#' DOI: https://doi.org/10.1214/21-AOAS1486
#'
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{G}}{A numeric vector containing the number of components from each AMIS iteration}
#'       \item{\code{Sigma}}{A list of covariace matrices from each component}
#'       \item{\code{Mean}}{A list of means from each component}
#'       \item{\code{probs}}{A list of probability weights for each component (unnormalised)}
#'     }
#' @param param A matrix containing the sampled parameter vectors.
#' @param prior_density Vector containing the prior density of each sampled parameter vector.
#' @param df The degrees of freedom for the t-distributed proposal distribution.
#' @param log A logical indicating whether to work on the log scale.
#' @return A vector containing the prior/proposal ratio for each row in \code{param}
compute_prior_proposal_ratio <- function(components, param, prior_density, df, log) {
  probs <- components$probs # /sum(unlist(components$probs)) # to normalise?
  Sigma <- components$Sigma
  Mean <- components$Mean
  G <- sum(components$G)
  q_terms<-matrix(NA,nrow(param),G)
  for (g in 1:G) {
    if (log) {
      q_terms[,g]<-log(probs[[g]])+mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=T)
    } else {
      q_terms[,g]<-probs[[g]]*mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=F)
    }
  }
  if (log) {
    M<-pmax(apply(q_terms,1,max),prior_density)
    return(prior_density - M - log(rowSums(exp(q_terms-M))+exp(prior_density-M)))
  } else {
    return(prior_density/(rowSums(q_terms)+prior_density)) 
  }
}

#' Compute the model evidence
#'
#' This function returns a Monte Carlo estimation of the model evidence 
#' for a given set of unnormalised parameter weights.
#'
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. 
#' @return estimate of the model evidence.
#' @noRd
compute_model_evidence <- function(likelihoods, amis_params, first_weight){
  n_tims <- dim(likelihoods)[1]
  n_locs <- dim(likelihoods)[2]
  n_sims <- dim(likelihoods)[3]
  weight_matrix <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs)
  # ## If n_locs = 1, likelihood matrix for given timepoint is an atomic
  # ## vector and doesn't need to be transposed.
  # if (n_locs == 1) {
  #   lik_matrix <- function(l) as.matrix(l)
  # } else {
  #   lik_matrix <- function(l) t(l)
  # }
  for (t in 1:n_tims) {
    # lik_mat <- lik_matrix(likelihoods[t,,])
    lik_mat <- t(array(likelihoods[t,,], dim=c(n_locs, n_sims)))
    
    # Update the weights by the latest likelihood (filtering)
    weight_matrix <- compute_weight_matrix_bayesian(lik_mat,amis_params,weight_matrix)
  }

  if(amis_params[['log']] == T){ 
    weight_matrix = exp(weight_matrix)
  }
  model_evidence = mean(weight_matrix)
  model_evidence_var = 1/length(weight_matrix)^2 * model_evidence^2 * sum((weight_matrix-1)^2)
  return(c(model_evidence = model_evidence,variance = model_evidence_var))
}

