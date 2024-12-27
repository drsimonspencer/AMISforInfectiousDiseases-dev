#' @useDynLib AMISforInfectiousDiseases
#' @importFrom Rcpp sourceCpp
NULL

#' Produce list containing the default AMIS parameters
#' 
#' For description of AMIS parameters, see argument \code{amis_params} in \code{\link{amis}()}.
#' @return List containing the default AMIS parameters.
#' @export
default_amis_params <- function() {
  amis_params <- list(n_samples=500, target_ess=500, max_iters=12,
                      boundaries=c(0,1), boundaries_param=NULL, 
                      log=TRUE, delete_induced_prior=FALSE, mixture_samples=1000, df=3, q=0,
                      delta=0.01, sigma=NULL, breaks=NULL)
  return(amis_params)
}

#' Check inputs of \code{\link{amis}()} function
#' 
#' Check whether all the inputs of \code{\link{amis}()} function are as expected.
#' @inheritParams amis
#' @noRd
# #' @export
check_inputs <- function(prevalence_map, transmission_model, prior, amis_params, seed, output_dir, initial_amis_vals) {

  if (!is.null(output_dir)){
    if(!is.character(output_dir)){
      stop("'output_dir' must be either NULL or a character string.")
    }  
    if(!dir.exists(output_dir)){dir.create(output_dir)}
    message("Outputs will be saved in the user-specified directory after each iteration (this will use data storage space).")
  }
  
  if (!is.null(initial_amis_vals)){
    if(!inherits(initial_amis_vals, "amis")){
      stop("'initial_amis_vals' must be either NULL or an object of class 'amis'.")
    }
  }
  
  if(is.data.frame(prevalence_map) || is.matrix(prevalence_map)){
    n_tims <- 1
  }else if(is.list(prevalence_map)) {
    n_tims <- length(prevalence_map)
  }else{
    stop(("prevalence_map must be either \n - a matrix or data frame of size #locations by #samples (for one timepoint); or \n - a list with n_tims timepoints, each one with a matrix named 'data'."))
  }
  
  if(is.data.frame(prevalence_map)){
    prevalence_map <- as.matrix(prevalence_map)
  }
  
  if (is.list(prevalence_map) && !is.data.frame(prevalence_map)) {
    if(!all(sapply(prevalence_map, function(i) inherits(i, "list"))) || 
       !all(sapply(prevalence_map, function(i) "data"%in%names(i))) ||
       !all(sapply(prevalence_map, function(i) all(names(i)%in%c("data","likelihood"))))
       ) {
      stop("Since prevalence_map was provided as a list, each of its elements must itself be a list with an object called 'data' (and optionally an object called 'likelihood').")
    }

    if(any( sapply(prevalence_map, function(i) "likelihood"%in%names(i)) )){
      if(!all( sapply(prevalence_map, function(i) "likelihood"%in%names(i)) )){
        stop("If 'likelihood' is to be used, it must be supplied for each time point in 'prevalence_map'.")  
      }
      for(i in 1:n_tims){
        stopifnot("'likelihood' must be a function." = is.function(prevalence_map[[i]]$likelihood))
      }
    }
  }
  
  # if matrix/data.frame, convert it to list
  if (is.matrix(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}

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
  breaks <- amis_params$breaks
  n_samples <- amis_params$n_samples
  mixture_samples <- amis_params$mixture_samples
  df <- amis_params$df
  target_ess <- amis_params$target_ess
  max_iters <- amis_params$max_iters
  delete_induced_prior <- amis_params$delete_induced_prior
  boundaries <- amis_params$boundaries
  boundaries_param <- amis_params$boundaries_param
  boundaries <- as.numeric(boundaries)
  q <- amis_params$q
  if(!(is.numeric(q)&&(q>=0)&&(q<=1))){
    stop("Parameter 'q' must be in [0,1].")
  }
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
                (ncol(boundaries_param)==2) && (is.null(initial_amis_vals) && nrow(boundaries_param)==ncol(prior$rprior(1))))
  }
  stopifnot("'delta' must be either NULL or a single positive numeric value" = ((length(delta)==1 && is.numeric(delta) && delta>0) || is.null(delta)))
  stopifnot("'sigma' must be either NULL or a single positive numeric value" = ((length(sigma)==1 && is.numeric(sigma) && sigma>0) || is.null(sigma)))
  stopifnot("'n_samples' must be a single numeric value greater than 1" = (length(n_samples)==1 && is.numeric(n_samples) && n_samples>1))
  stopifnot("'mixture_samples' must be a single numeric value" = (length(mixture_samples)==1 && is.numeric(mixture_samples)))
  stopifnot("'df' must be a single numeric value greater than 0" = (length(df)==1 && is.numeric(df) && df>0))
  stopifnot("'target_ess' must be a single numeric value greater than 0" = (length(target_ess)==1 && is.numeric(target_ess) && target_ess>0))
  stopifnot("'max_iters' must be a single numeric value greater than 1" = (length(max_iters)==1 && is.numeric(max_iters) && max_iters>1))
  stopifnot("'seed' must be either NULL or a single numeric value" = ((length(seed)==1 && is.numeric(seed)) || is.null(seed)))
  if(is.null(prevalence_map[[1]]$likelihood) && is.null(c(delta, sigma, breaks))){
    stop("At least one of the inputs ('delta','sigma','breaks') must not be NULL if a likelihood function is not provided.")
  }
  if(delete_induced_prior && is.null(c(delta, sigma, breaks))){
    stop("At least one of the inputs ('delta','sigma','breaks') must not be NULL if 'delete_induced_prior' is set to FALSE.")
  }
  
  mes <- NULL
  mes_ <- NULL
  if(is.null(prevalence_map[[1]]$likelihood) && !delete_induced_prior){
    if(!is.null(breaks)){
      mes_ <- "- Histogram method will be used in the estimation of the likelihood as 'breaks' was provided."
    }else{
      if(!is.null(sigma)){
        mes_ <- "- Gaussian kernel will be used in the estimation of the likelihood as 'sigma' was provided."
      }else{
        mes_ <- "- Uniform kernel will be used in the estimation of the likelihood."
      }
    }
  }
  mes <- c(mes, mes_)
  if(!delete_induced_prior){
    mes_ <- "- Induced prior will not be deleted in the update of the weights."
  }else{
    if(is.null(prevalence_map[[1]]$likelihood)){
      if(!is.null(breaks)){
        mes_ <- "- Histogram method will be used to estimate the likelihood and induced prior densities as 'breaks' was provided."
      }else{
        if(!is.null(sigma)){
          mes_ <- "- Gaussian kernel will be used to estimate the likelihood and induced prior densities as 'sigma' was provided."
        }else{
          mes_ <- "- Uniform kernel will be used to estimate the likelihood and induced prior densities."
        }
      }
    }else{
      mes_0 <- "- A likelihood function was provided and will be used to calculate likelihood terms."
      if(!is.null(breaks)){
        mes_ <- "- Histogram method will be used to estimate the induced prior density as 'breaks' was provided."
      }else{
        if(!is.null(sigma)){
          mes_ <- "- Gaussian kernel will be used to estimate the induced prior density as 'sigma' was provided."
        }else{
          mes_ <- "- Uniform kernel will be used to estimate the induced prior density."
        }
      }
      mes_ <- c(mes_0, mes_)
    }
  }
  mes <- c(mes, mes_)
  
  locs_no_data <- NULL
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  likelihood_approach <- ifelse(is.null(prevalence_map[[1]]$likelihood), "nonparametric", "parametric")
  if(likelihood_approach=="nonparametric"){
    for(l in 1:n_locs){
      num_valid_datapts <- 0L
      for(t in 1:n_tims){
        data_l_t <- unlist(prevalence_map[[t]]$data[l,])
        num_valid_datapts <- num_valid_datapts + sum(is.finite(data_l_t) & (data_l_t>=boundaries[1]) & (data_l_t<=boundaries[2]))
      }
      if(num_valid_datapts==0L){
        locs_no_data <- c(locs_no_data, l)
      }
    }
    if(!is.null(locs_no_data)){
      if(length(locs_no_data)==1L){
        message(paste0("Location ", shQuote(locs_no_data), " has no valid map samples. Its corresponding posterior samples will therefore be driven by the prior."))
      }else if(length(locs_no_data)>1L && length(locs_no_data)<=10L){
        message(paste0("Locations (", paste(shQuote(locs_no_data), collapse=","), ") have no valid map samples. Their corresponding posterior samples will therefore be driven by the prior."))
      }else{
        message(paste0(length(locs_no_data), " locations have no valid map samples. Their corresponding posterior samples will therefore be driven by the prior."))
      }
    }
  }
  return(list(locs_no_data=locs_no_data, messages=mes)) 
}

#' Print error message if all particles are assigned weight zero for locations in the active set
#' 
#' @inheritParams amis
#' @param likelihood_approach A string saying the likelihood approach: "parametric" or "nonparametric".
#' @param nonparametric_method A string specifying which nonparametric method are used. Default to \code{NULL}.
#' @noRd
# #' @export
check_zero_weight_for_all_particles <- function(amis_params, mean_weights, likelihood_approach, nonparametric_method=NULL){
  if ((amis_params[["log"]] && max(mean_weights)==-Inf) || (!amis_params[["log"]] && max(mean_weights)==0)) {
    if(likelihood_approach=="nonparametric"){
      prefix <- "No weight on any particles for locations in the active set. "
      if(nonparametric_method=="histogram"){
        errorMessage <- paste0(prefix, "Try to use breakpoints more distant to each other.")
      }
      if(nonparametric_method=="gaussian"){
        errorMessage <- paste0(prefix, "Try to use larger sigma.")
      }
      if(nonparametric_method=="uniform"){
        errorMessage <- paste0(prefix, "Try to use larger delta.")
      }
      stop(errorMessage)
    }else{
      stop("No weight on any particles for locations in the active set. Check whether the prevalence map data are correct for all locations.")
    }
  }
}

#' Compute likelihood for each additional simulation across timepoints
#'
#' Calls evaluate likelihood for each timepoint.
#' @param prevalence_map A list with one entry for each timepoint.
#' Each entry must be a list containing objects \code{data} (an L x M matrix of data);
#' and optional function \code{likelihood} (see \code{\link{amis}()})
#' @param simulated_prevalences An n x timepoints matrix of prevalences simulated from the transmission model.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}()}.
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations (optional).
#' @param which_valid_sim_prev_iter List of T elements, where each one indicates which simulated prevalences are valid at the current iteration
#' @param which_valid_prev_map List showing which prevalence map samples are valid
#' @param log_norm_const_gaussian Normalising constant (in log scale) for the Gaussian kernel. It is only used if Gaussian kernel is used.
#' @return A larger array with the likelihoods of the new simulations joined to the existing array \code{likelihoods}.
#' @noRd
# #' @export
compute_likelihood <- function(prevalence_map,simulated_prevalences,amis_params,likelihoods=NULL,
                               which_valid_sim_prev_iter,which_valid_prev_map,log_norm_const_gaussian) {
  n_tims <- length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  n_sims <- dim(simulated_prevalences)[1]
  lik <- array(NA, c(n_tims,n_locs,n_sims)) # this way around to avoid abind -- different to elsewhere
  for (t in 1:n_tims) {
    if(!is.null(amis_params[["sigma"]])){
      log_norm_const_gaussian_t <- log_norm_const_gaussian[t,,]
      if(n_locs==1){
        log_norm_const_gaussian_t <- t(log_norm_const_gaussian_t)  
      }
    }
    lik[t,,] <- evaluate_likelihood(prevalence_map[[t]],simulated_prevalences[,t],amis_params,
                                    which_valid_sim_prev_iter[[t]],which_valid_prev_map[[t]],log_norm_const_gaussian_t) 
  }
  if (!is.null(likelihoods)) {
    lik <- array(c(likelihoods,lik), c(n_tims,n_locs,dim(likelihoods)[3]+n_sims))
  }
  return(lik)
}

#' Evaluate likelihood for each additional simulation for a single timepoint
#'
#' Implements a parametric likelihood function if it is supplied by the user; otherwise it implements a 
#' nonparametric method (histogram or kernel density smoothing) based on samples from a geostatistical map.
#' @param prevalence_map A list containing objects \code{data} (an L x M matrix of data);
#' and \code{likelihood} a function taking arguments \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.    
#' @param prev_sim A vector of simulated prevalences
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}()}.
#' @param which_valid_sim_prev_iter Vector showing which simulated prevalences are valid at the current iteration at time t.
#' @param which_valid_prev_map_t List showing which prevalence map samples are valid at time t.
#' @param log_norm_const_gaussian_t Normalising constant (in log scale) for the Gaussian kernel at time t. It is only used if Gaussian kernel is used.
#' @return A locations x simulations matrix of (log-)likelihoods.
#' @noRd
# #' @export
evaluate_likelihood <- function(prevalence_map,prev_sim,amis_params, 
                                which_valid_sim_prev_iter,which_valid_prev_map_t,log_norm_const_gaussian_t) {

  logar <- amis_params[["log"]]
  # f <- matrix(NA,dim(prevalence_map$data)[1],length(prev_sim))
  
  if (!is.null(prevalence_map$likelihood)) {
    
    likelihood_fun <- prevalence_map$likelihood
    
    f <- f_user_defined(likelihood_fun=likelihood_fun, 
                        prevalence_map=prevalence_map$data, 
                        which_valid_prev_map_t=which_valid_prev_map_t,
                        prev_sim=prev_sim, 
                        which_valid_sim_prev_iter=which_valid_sim_prev_iter, 
                        logar=logar)

    # # R code of previous version of the package (note that 'param' had to be an argument of evaluate_likelihood().)
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
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}()}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as the rows in \code{simulated_prevalence}.
#' @param locs_with_g List indicating, at each time, which locations are updated using induced prior.
#' @param locs_without_g List indicating, at each time, which locations are updated without using induced prior.
#' @param bool_valid_sim_prev Matrix of n_tims columns, where column is a logical vector indicating which simulated prevalences are valid.
#' @param which_valid_sim_prev List indicating, at each time, which simulated prevalences are valid.
#' @param which_invalid_sim_prev List indicating, at each time, which simulated prevalences are invalid
#' @param which_valid_locs_prev_map List showing which locations have valid data at each time
#' @param locations_with_no_data Vector indicating which locations have no data at any time point
#' @return normalised weight matrix.
#' @noRd
# #' @export
compute_weight_matrix <- function(likelihoods, simulated_prevalence, amis_params, first_weight, 
                                  locs_with_g, locs_without_g,
                                  bool_valid_sim_prev, which_valid_sim_prev, which_invalid_sim_prev, 
                                  which_valid_locs_prev_map, locations_with_no_data) {

  n_tims <- dim(likelihoods)[1]
  n_locs <- dim(likelihoods)[2]
  n_sims <- dim(likelihoods)[3]
  weight_matrix <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs)

  for (t in 1:n_tims) {

    lik_mat <- t(array(likelihoods[t,,], dim=c(n_locs, n_sims)))
    
    # Update the weights by the latest likelihood (filtering)
    if (amis_params[["delete_induced_prior"]]){
      
      # If this is the first timepoint where there is data for a location, then use induced prior
      # locs_with_g = which(locations_first_t == t)
      # locs_without_g = which(locations_first_t < t)

      if(!is.null(locs_without_g[[t]])){
        weight_matrix <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix,
                                                         which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]],
                                                         locs_without_g[[t]])
      }

      if (is.null(amis_params[["breaks"]])){
          if(is.null(amis_params[["sigma"]])){
            if(!is.null(locs_with_g[[t]])){
              weight_matrix <- compute_weight_matrix_empirical_uniform(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix,
                                                                       bool_valid_sim_prev[,t], 
                                                                       which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                       locs_with_g[[t]])
            }
          }else{
            if(!is.null(locs_with_g[[t]])){
              weight_matrix <- compute_weight_matrix_empirical_gauss(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix, 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
            }
          }
      } else {
        if(!is.null(locs_with_g[[t]])){
          weight_matrix <- compute_weight_matrix_empirical_histogram(lik_mat,simulated_prevalence[,t],amis_params,weight_matrix,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
        }
      }
    } else {
      if(!is.null(which_valid_locs_prev_map[[t]])){
        weight_matrix <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix,
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






#' Compute the current effective sample size
#'
#' This function returns the effective sample size (ESS) for each location. 
#'  For each column of the weight matrix \code{weight_mat}, the ESS component is computed as
#' \deqn{(\sum_{i=1}^{N} w_{i}^2 )^{-1}}
#' where \eqn{N} is the number of sampled parameter values for each location.
#'
#' @param weight_mat The weight matrix. A \eqn{N \times L} matrix with L the number of locations
#'     and N the number of sampled parameter values.
#' @param log logical indicating if the weights are on the log scale. 
#' @return A vector containing the ESS value for each location.
#'
#' @seealso \code{\link{compute_weight_matrix}()}.
#' @noRd
# #' @export
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
#'     \code{\link{compute_weight_matrix}()}
#' @param ess The effective sample size vector as returned by
#'     \code{\link{calculate_ess}()}
#' @param target_size A number representing the target size for the sample.
#' @param log A logical indicating if the weights are logged.
#' @param q Parameter (between 0 and 1) controlling how the weights are calculated for active locations. 
#' @return Vector containing the row sums of the active columns of the weight matrix.
#' @noRd
# #' @export
update_according_to_ess_value <- function(weight_matrix, ess, target_size, log, q) {
  n <- nrow(weight_matrix)
  active_cols <- which(ess < target_size)
  if (log) {
    ## original version:
    # M<-apply(weight_matrix[,active_cols,drop=FALSE],1,max)
    # wh<-which(M==-Inf)
    # M[wh]<-0
    # out <- M+log(rowSums(exp(weight_matrix[,active_cols,drop=FALSE]-M)))
    
    ## re-scaling depending on q:
    diffs <- ((target_size - ess)^q)[active_cols]
    ratio <- log(diffs/sum(diffs))                                        # an |A|-length vector
    if(length(active_cols)==1){
      ratio <- as.matrix(rep(ratio, n))
    }else{
      ratio <- t(replicate(n, ratio))
    }
    # ratio is N by |A|

    new_weight_matrix <- weight_matrix[,active_cols,drop=FALSE] + ratio   # both in log-scale
    M <- apply(new_weight_matrix,1,max)
    wh <- which(M==-Inf)
    M[wh] <- 0
    out <- M + log(rowSums(exp(new_weight_matrix-M)))
  } else {
    ## original version:
    # out <- rowSums(weight_matrix[,active_cols,drop=FALSE])
    
    ## re-scaling depending on q:
    diffs <- ((target_size - ess)^q)[active_cols]
    ratio <- diffs/sum(diffs)
    if(length(active_cols)==1){
      ratio <- as.matrix(rep(ratio, n))
    }else{
      ratio <- t(replicate(n, ratio))
    }
    new_weight_matrix <- weight_matrix[,active_cols,drop=FALSE] * ratio
    out <- rowSums(new_weight_matrix)
  }
  return(out)
}
#' Systematic resampling function
#' 
#' Implement systematic resampling to reduce variance in weighted particle selection
#' @param n_samples number of samples to draw
#' @param weights vector of length equal to the number of particles, containing their weights
#' @param log logical indicating if weights are log-weights
#' @return vector of indices of the sampled particles 
#' @noRd
# #' @export
systematic_sample <- function(n_samples,weights,log=FALSE) {
  if (log) {
    M<-max(weights)
    log_sum_weights<-M+log(sum(exp(weights-M)))
    cum <- cumsum(exp(weights-log_sum_weights))
  } else {
    cum <- cumsum(weights)/sum(weights) # cumulative sum of normalised weights
  }
  u <- stats::runif(1)/n_samples+0:(n_samples-1)/n_samples
  return(1+matrix(rep(u,length(weights))>rep(cum,each=n_samples),n_samples,length(weights))%*%matrix(1,length(weights),1))
}
#' Fit mixture to weighted sample
#' 
#' Weights are implemented by using systematic resampling to obtain an unweighted set of parameters.
#' An unweighted mixture is then fitted using \code{\link{fit_mixture}()}.
#' @param parameters An N x d matrix containing the sampled values for the d parameters.
#' @param n_samples The number of parameter to resample as data to fit the mixture to.
#' @param weights A vector of weights with length N.
#' @param log logical indicating if weights are logged.
#' @return A list of the mixture components (see function \code{\link{fit_mixture}()})
#'     \describe{
#'       \item{\code{probs}}{The mixture weights}
#'       \item{\code{Mean}}{The means of the components}
#'       \item{\code{Sigma}}{The covariance matrices of the components}
#'       \item{\code{G}}{Number of components}
#'       \item{\code{BIC}}{BIC of fitted mixture}
#'       \item{\code{ModelName}}{Model name from package mclust}
#'     }
#'
#' @seealso \code{\link{fit_mixture}()}
#' @noRd
# #' @export
weighted_mixture <- function(parameters, n_samples, weights, log=FALSE) {
  sampled_idx <- systematic_sample(n_samples,weights,log)
  if (length(unique(sampled_idx))==1) {warning("Only one particle with sufficient weight. Will result in a non-invertible covariance matrix for the mixture.")}
  return(fit_mixture(parameters[sampled_idx,,drop=FALSE]))
}
#' Sample new parameters
#'
#' This function generates \code{"n_samples"} new model parameter values according the
#' t-distribution with \code{"df"} degrees of freedom and the mixture components \code{"mixture"}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{weighted_mixture}()}
#' @param n_samples A number of new parameters to sample (integer)
#' @param df The degrees of freedom for the t-distributed proposal distribution.
#' @param prior list containing the functions \code{rprior} and \code{dprior}
#' @param log A logical indicating if densities 
#' 
#' @return A list containing the sampled parameters:
#' \describe{
#'   \item{\code{params}}{An \code{n_samples} x d matrix containing the sampled parameter values}
#'   \item{\code{prior_density}}{The corresponding vector of prior densities}
#'   \item{\code{compon_proposal}}{A vector indicating the mixture component each parameter value was simulated from}
#' }
#'
#' @seealso \code{\link{fit_mixture}()}
#' @noRd
# #' @export
sample_new_parameters <- function(mixture, n_samples, df, prior, log) {
  prior_density <- rep(NA,n_samples)
  compon_proposal <- rep(NA,n_samples)
  i<-1
  while (i <= n_samples) {
    compo <- sample(1:mixture$G, 1, prob = mixture$probs)
    compon_proposal[i] <- compo
    proposal <- mnormt::rmt(1,mean=mixture$Mean[,compo],S=mixture$Sigma[,,compo],df=df)
    density_of_proposal <- prior$dprior(proposal,log=log)
    if (all(!is.na(proposal)) && ((log && density_of_proposal>-Inf) || (!log && density_of_proposal>0))) {
      if (i==1) {params<-matrix(NA,n_samples,length(proposal))}
      params[i,]<-proposal
      prior_density[i]<-density_of_proposal
      i<-i+1
    }
  }
  return(list(params=params,prior_density=prior_density,compon_proposal=compon_proposal))
}

#' Update the components of the mixture
#'
#' This function updates the mixture \code{"components"} according to
#' the current mixture \code{"mixture"} generated at iteration \code{"iter"}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{fit_mixture}()}
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
#' @seealso \code{\link{weighted_mixture}()}, \code{\link{fit_mixture}()}
#' @noRd
# #' @export
update_mixture_components <- function(mixture, components, iter) {
  components$G[iter] <- mixture$G
  G_previous <- sum(components$G[1:(iter - 1)]) # Number of pre-existing components
  for (i in 1:mixture$G) {
    components$Sigma[[i + G_previous]] <- mixture$Sigma[, , i]
    components$Mean[[i + G_previous]] <- mixture$Mean[,i]
    components$probs[[i + G_previous]] <- mixture$probs[i] ### scale by number of points if n_samples varies by iteration
  }
  return(components)
}

#' Update the \code{Mclust} object used for plotting mixture components
#' 
#' The object is used in \code{\link{plot_mixture_components}()}.
#' 
#' @param mixture A list of the mixture components returned by \code{\link{weighted_mixture}()}
#' @param sampled_params List of parameter values sampled from the mixture. It is returned by \code{\link{sample_new_parameters}()}
#' @return An updated list of class \code{Mclust}
#' @noRd
# #' @export
update_Mclust_object <- function(mixture, sampled_params){
  if(!inherits(mixture$clustering, "Mclust")){
    stop("'mixture' must have an element called 'clustering' of class 'Mclust'.")
  }
  G <- mixture$G
  mixture$clustering$G <- G
  mixture$clustering$parameters$pro <- mixture$probs
  mixture$clustering$parameters$mean <- mixture$Mean
  mixture$clustering$parameters$variance$G <- G
  d <- mixture$clustering$d
  mixture$clustering$parameters$variance$sigma <- mixture$Sigma
  cholsigma_list <- lapply(1:G, function(g) chol(mixture$Sigma[,,g]))
  mixture$clustering$parameters$variance$cholsigma <- array(as.numeric(unlist(cholsigma_list)), dim=c(d, d, G))
  mixture$clustering$data_proposed <- sampled_params$params
  mixture$clustering$compon_proposal <- sampled_params$compon_proposal
  return(mixture$clustering)
}

#' Compute the prior/proposal ratio
#'
#' This function returns the ratio between the prior and proposal distribution
#' for each sampled parameter value (i.e. each row in \code{param}).
#' This function returns the first weight
#' See step (4) of the AMIS algorithm in 
#' Retkute, R., Touloupou, P., Basanez, M. G., Hollingsworth, T. D., 
#' Spencer, S. E. (2021). \emph{Integrating geostatistical maps and infectious disease 
#' transmission models using adaptive multiple importance sampling.} 
#' The Annals of Applied Statistics, 15(4), 1980-1998. 
#' DOI: \url{https://doi.org/10.1214/21-AOAS1486}.
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
#' @noRd
# #' @export
compute_prior_proposal_ratio <- function(components, param, prior_density, df, log) {
  probs <- components$probs # /sum(unlist(components$probs)) # to normalise?
  Sigma <- components$Sigma
  Mean <- components$Mean
  G <- sum(components$G)
  q_terms<-matrix(NA,nrow(param),G)
  for (g in 1:G) {
    if (log) {
      q_terms[,g]<-log(probs[[g]])+mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=TRUE)
    } else {
      q_terms[,g]<-probs[[g]]*mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=FALSE)
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
#' @param simulated_prevalence An n_sims x n_tims matrix containing the simulated prevalence values for each of the
#'     parameter samples. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}()}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as the rows in \code{simulated_prevalences}.
#' @param locs_with_g List indicating, at each time, which locations are updated using induced prior.
#' @param locs_without_g List indicating, at each time, which locations are updated without using induced prior.
#' @param bool_valid_sim_prev Matrix of n_tims columns, where column is a logical vector indicating which simulated prevalences are valid.
#' @param which_valid_sim_prev List indicating, at each time, which simulated prevalences are valid.
#' @param which_invalid_sim_prev List indicating, at each time, which simulated prevalences are invalid
#' @param which_valid_locs_prev_map List showing which locations have valid data at each time
#' @param locations_with_no_data Vector indicating which locations have no data at any time point
#' @return A list containing an estimate of the log model evidence and corresponding log variance of this estimate for both the full likelihood model 
#'     (product over all locations), and for each location individually.
#' @noRd
# #' @export
compute_model_evidence <- function(likelihoods, simulated_prevalences, 
                                   amis_params, first_weight,
                                   locs_with_g, locs_without_g,
                                   bool_valid_sim_prev, which_valid_sim_prev, 
                                   which_invalid_sim_prev, which_valid_locs_prev_map, 
                                   locations_with_no_data){
  
  n_tims <- dim(likelihoods)[1]
  n_locs <- dim(likelihoods)[2]
  n_sims <- dim(likelihoods)[3]
  weight_matrix <- matrix(rep(rep(1-amis_params[["log"]], amis_params[["n_samples"]]),n_locs), nrow = n_sims, ncol = n_locs) 
  weight_matrix_loc <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs) 
  # ## If n_locs = 1, likelihood matrix for given timepoint is an atomic
  # ## vector and doesn't need to be transposed.
  # if (n_locs == 1) {
  #   lik_matrix <- function(l) as.matrix(l)
  # } else {
  #   lik_matrix <- function(l) t(l)
  # }
  
  for (t in 1:n_tims) {
    
    lik_mat <- t(array(likelihoods[t,,], dim=c(n_locs, n_sims)))
    
    # Update the weights by the latest likelihood (filtering)
    if (amis_params[["delete_induced_prior"]]){
      
      # If this is the first timepoint where there is data for a location, then use induced prior
      # locs_with_g = which(locations_first_t == t)
      # locs_without_g = which(locations_first_t < t)
      
      if(!is.null(locs_without_g[[t]])){
        weight_matrix <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix,
                                                         which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]],
                                                         locs_without_g[[t]])
        weight_matrix_loc <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix_loc,
                                                         which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]],
                                                         locs_without_g[[t]])      }
      
      if (is.null(amis_params[["breaks"]])){
        if(is.null(amis_params[["sigma"]])){
          if(!is.null(locs_with_g[[t]])){
            weight_matrix <- compute_weight_matrix_empirical_uniform(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
            weight_matrix_loc <- compute_weight_matrix_empirical_uniform(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix_loc,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
          }
        }else{
          if(!is.null(locs_with_g[[t]])){
            weight_matrix <- compute_weight_matrix_empirical_gauss(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix, 
                                                                   which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                   locs_with_g[[t]])
            weight_matrix_loc <- compute_weight_matrix_empirical_gauss(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix_loc, 
                                                                   which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                   locs_with_g[[t]])
          }
        }
      } else {
        if(!is.null(locs_with_g[[t]])){
          weight_matrix <- compute_weight_matrix_empirical_histogram(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
          weight_matrix_loc <- compute_weight_matrix_empirical_histogram(lik_mat,simulated_prevalences[,t],amis_params,weight_matrix_loc,
                                                                     bool_valid_sim_prev[,t], 
                                                                     which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                                     locs_with_g[[t]])
        }
      }
    } else {
      if(!is.null(which_valid_locs_prev_map[[t]])){
        weight_matrix <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix,
                                                         which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                         which_valid_locs_prev_map[[t]])
        weight_matrix_loc <- compute_weight_matrix_without_g(lik_mat, amis_params, weight_matrix_loc,
                                                         which_valid_sim_prev[[t]], which_invalid_sim_prev[[t]], 
                                                         which_valid_locs_prev_map[[t]])
      }
    }
    
    if(length(locations_with_no_data)>0 && length(which_invalid_sim_prev[[t]])>0){
      weight_inval_prev <- ifelse(amis_params[["log"]], -Inf, 0)
      weight_matrix[which_invalid_sim_prev[[t]]+1L, locations_with_no_data] <- weight_inval_prev
      weight_matrix_loc[which_invalid_sim_prev[[t]]+1L, locations_with_no_data] <- weight_inval_prev
    }
  }
  
  # Model evidence of full model
  joint_log_posterior <- rowSums(weight_matrix) + first_weight
  M <- max(joint_log_posterior)
  log_model_evidence <-  -log(n_sims) + M + log(sum(exp(joint_log_posterior - M))) 
  M_var <- max((2*joint_log_posterior),(log(2)+joint_log_posterior),0)
  log_model_evidence_var <- -2*log(n_sims) + M_var + log(sum(exp(2*joint_log_posterior - M_var),(-2*exp(joint_log_posterior - M_var)), n_sims*exp(-M_var)))
  
  # Model evidence for each location subset
  joint_log_posterior_loc <- weight_matrix_loc
  M_loc <- apply(joint_log_posterior_loc,2,max)
  log_model_evidence_loc <- sapply(1:n_locs, function(v) {
    -log(n_sims) + M_loc[v] + log(sum(exp(joint_log_posterior_loc[,v] - M_loc[v])))
  })
  log_model_evidence_var_loc <- sapply(1:n_locs, function(v) {
    -2*log(n_sims) + 2*M_loc[v] + log(sum(exp(2*joint_log_posterior_loc[,v] - 2*M_loc[v]),(-2*exp(joint_log_posterior_loc[,v] - 2*M_loc[v])), n_sims*exp(-2*M_loc[v])))
  })
  
  return(list(evidence = cbind(log_model_evidence = log_model_evidence, log_variance = log_model_evidence_var), 
              joint_log_posterior = joint_log_posterior,
              evidence_by_location = cbind(log_model_evidence = log_model_evidence_loc, log_variance = log_model_evidence_var_loc), 
              joint_log_posterior_by_location = joint_log_posterior_loc))
}

