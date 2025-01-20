#' Run the AMIS algorithm to fit a transmission model to a map
#' 
#' This implements the AMIS algorithm as described in \cite{Retkute et al. (2021)}. 
#' Each iteration of the algorithm produces a set of parameters from a proposal 
#' distribution (the prior in the first iteration). For each parameter set, a simulation 
#' is run from the transmission model. Then, each preceding simulation is weighted at 
#' each location according to the distribution of prevalences (or likelihood function) 
#' at that location. A Gaussian mixture model is then fitted to the parameter samples 
#' with weights averaged over the active locations (ie locations that have yet to reach 
#' the effective sample size target). This Gaussian mixture informs the proposal for the 
#' next iteration. The algorithm continues until every location has reached the ESS target, 
#' or the maximum number of iterations is reached.
#'
#' @references Retkute, R., Touloupou, P., Basanez, M. G., Hollingsworth, T. D., 
#' Spencer, S. E. (2021). \emph{Integrating geostatistical maps and infectious disease 
#' transmission models using adaptive multiple importance sampling.} 
#' The Annals of Applied Statistics, 15(4), 1980-1998. \doi{10.1214/21-AOAS1486}.
#' @param prevalence_map For a single timepoint, \code{"prevalence_map"} can be an \eqn{L \times M} matrix or data frame containing samples from a geostatistical model, 
#' where \eqn{L} is the number of locations and \eqn{M} the number of samples per location.
#' \cr \cr If there are multiple timepoints and/or a parametric \code{likelihood} function is to be used, \code{"prevalence_map"} must be a list with \eqn{T} elements, one for each timepoint \eqn{t=1,\dots,T}.
#' Each element must itself be a list with the following objects: 
#' \describe{
#' \item{\code{data}}{An \eqn{L \times M} matrix as above}
#' \item{\code{likelihood}}{(optional) A function taking arguments:
#' \itemize{
#'    \item \code{data}: A vector of length \eqn{M_l}, where \eqn{M_l} 
#' is the number of samples from a geostatistical model for location \eqn{l} or 
#' the number of likelihood parameters;
#'    \item \code{sim_prev}: A numeric value for a prevalence simulated from the transmission model;
#'    \item \code{log}: Logical indicating if calculations are to be performed on log scale 
#'    (specified in \code{"amis_params"}, see below).
#' }
#' The function \code{likelihood} must return a numeric value representing the (log-)likelihood of 
#' observing a simulated prevalence given the data from a particular location.
#' }
#' }
#' The location names are inherited from \code{rownames(prevalence_map)} 
#' if \code{"prevalence_map"} is a matrix, and 
#' from \code{rownames(prevalence_map[[1]]$data)} if \code{"prevalence_map"} is a list.
#' \cr \cr If \code{likelihood} is not specified, then it is assumed that the data consist of 
#' samples from a geostatistical model and a nonparametric method is used. The nonparametric 
#' method to be used is specified in \code{"amis_params"} using the options 
#' \code{breaks}, \code{delta}, or \code{sigma} (see \code{"amis_params"}).
#' \cr \cr
#' @param transmission_model A function taking arguments:
#' \itemize{
#'    \item \code{seeds}: a vector of \eqn{n} seeds;
#'    \item \code{params}: an \eqn{n \times d} matrix of parameter vectors;
#'    \item \code{n_tims}: number of time points.
#'  }
#' This function must return an \eqn{n \times T} \bold{matrix} of prevalences 
#' (it must be a matrix even when \eqn{T=1}). The vector \code{seeds} will be the 
#' vector of indexes of the simulated samples. If \code{n_samples} new samples are 
#' drawn within each iteration of the AMIS algorithm, then the vector \code{seeds} will be 
#' \code{1:n_samples} at the first iteration, \code{(n_samples+1):(2*n_samples)} at 
#' the second iteration, and so on.
#' @param prior A list containing the functions \code{dprior} and \code{rprior} 
#' (density and random number generator, respectively).
#' The two arguments of \code{dprior} must be:
#' \itemize{
#'   \item a \eqn{d}-length vector of transmission model parameters; and
#'   \item a logical \code{log} to indicate whether to calculate log-density or not. 
#' }
#' The only argument of \code{rprior} must be a single integer \eqn{n} that determines the number of samples to draw. 
#' \code{rprior} must produce an \eqn{n \times d} \bold{matrix} of parameters even when \eqn{d=1}.
#' Parameter names are inherited from the \code{colnames} of the output of \code{rprior}.
#' @param amis_params A list containing control parameters for the AMIS algorithm 
#' (\code{default_amis_params()} returns the default values):
#' \describe{
#' \item{\code{n_samples}}{Number of new samples drawn within each AMIS iteration. Default to \code{500}.}
#' \item{\code{target_ess}}{Target effective sample size. Default to \code{500}.}
#' \item{\code{max_iters}}{Maximum number of AMIS iterations. Default to \code{12}.}
#' \item{\code{boundaries}}{A vector of length two with the left and right boundaries for prevalences. 
#' Default to \code{c(0,1)}. If, for instance, left boundary is zero and there is no right boundary, 
#' set \code{boundaries = c(0,Inf)}.}
#' \item{\code{boundaries_param}}{If specified, it should be a \eqn{d \times 2} matrix 
#' with the lower and upper boundaries for the \eqn{d} transmission model parameters. Default to \code{NULL}.}
#' \item{\code{log}}{Logical indicating if calculations are to be performed on log scale. Default to \code{TRUE}.}
#' \item{\code{delete_induced_prior}}{Logical indicating whether the induced prior density is to be deleted in the update of weights. Default to \code{FALSE}.}
#' \item{\code{mixture_samples}}{Number of samples used to represent the weighted parameters in the mixture fitting.}
#' \item{\code{df}}{Degrees of freedom in the \eqn{t}-distributions, used to yield a heavy tailed proposal. Default to \code{3}.}
#' \item{\code{q}}{Parameter (between 0 and 1) controlling how the weights are calculated for active locations. 
#' Default to \code{0}. See Details below.}
#' \item{\code{delta}}{Optional smoothing parameter if uniform kernel (default) is used. Default to \code{0.01}.}
#' \item{\code{sigma}}{Optional smoothing parameter if Gaussian kernel is used. Default to \code{NULL}.}
#' \item{\code{breaks}}{Optional vector specifying the breaks for the histogram. Default to \code{NULL}.
#' For finite \code{boundaries}, the first and last entries of \code{breaks} must be 
#' equal to the left and right boundaries, respectively.
#' For non-finite \code{boundaries}, ensure that the range of \code{breaks} includes any possible prevalence value.}
#' }
#' Uniform kernel is the default method for the density estimator of the likelihood. 
#' If \code{sigma} is provided, then Gaussian kernel will be used instead. 
#' If \code{breaks} is provided, then histogram-based method will be the nonparametric method being used. 
#' Note that if \code{likelihood} is provided in \code{prevalence_map}, then a parametric method will be implemented.
#' @param seed Optional single value interpreted as an integer. 
#' It is the seed for the random number generator for the AMIS algorithm. This is not the same as
#' the \code{seeds} argument passed to \code{"transmission_model"}.
#' @param output_dir A string specifying the local directory where to save outputs 
#' after each iteration of the algorithm. At the end of the string, 
#' use the correct path separator for your machine's operating system. 
#' If the directory is specified, the outputs will be saved in a file called `amis_output.rds'. 
#' Default to \code{NULL} (i.e. outputs are not saved in a local directory).
#' @param initial_amis_vals Optional list of intermittent outputs from a 
#' previous run (where at least one iteration was successful). These outputs can 
#' be saved by specifying the directory \code{"output_dir"}. 
#' @return A list of class \code{amis}. If the algorithm completed \eqn{I} iterations, 
#' it simulated a total of \eqn{N = I \times} \code{n_samples}, and therefore the list returned by \code{amis()} will contain:
#' \describe{
#' \item{\code{seeds}}{An \eqn{N}-length vector with the simulation seeds that were used.}
#' \item{\code{param}}{An \eqn{N \times d} matrix with the \eqn{d}-dimensional transmission model parameters simulated by the algorithm.}
#' \item{\code{simulated_prevalences}}{An \eqn{N \times T} matrix with the simulated prevalences, where \eqn{T} is the number of timepoints.}
#' \item{\code{weight_matrix}}{An \eqn{N \times L}, where \eqn{L} is the number of locations.}
#' \item{\code{likelihoods}}{A \eqn{T \times L \times N} array with the likelihood of observing a simulated prevalence in each location at each time.}
#' \item{\code{ess}}{An \eqn{L}-length vector with the final effective sample size (ESS) for each location.}
#' \item{\code{prevalence_map}}{List with the prevalence map supplied by the user.}
#' \item{\code{locations_with_no_data}}{Vector indicating which locations have no data at any time point.}
#' \item{\code{components}}{A list of the mixture components of all iterations, containing:
#'  \itemize{
#'    \item \code{G}: number of components in each iteration;
#'    \item \code{probs}: the mixture weights;
#'    \item \code{Mean}: the locations of the components;
#'    \item \code{Sigma}: the covariance matrices of the components.
#'  }
#' }
#' \item{\code{components_per_iteration}}{A list with the mixture components at each iteration. 
#' This object is used in \code{\link{plot_mixture_components}()}.}
#' \item{\code{ess_per_iteration}}{An \eqn{L \times I} matrix with with the ESS for each location after each iteration.}
#' \item{\code{prior_density}}{An \eqn{N}-length vector with the density function evaluated at the simulated parameter values.}
#' \item{\code{amis_params}}{List supplied by the user.}
#' \item{\code{evidence}}{A list containing an estimate of the log model evidence and corresponding log variance of this estimate for both the full likelihood model 
#'     (product over all locations), and for each location individually.}
#' }
#' @details The average weight of parameter vectors for the set of active locations at iteration \eqn{i} \eqn{\left(A_i\right)}
#' has weights determined by how far the effective sample size for location \eqn{l} \eqn{\left(\text{ESS}_l^i\right)} 
#' is from the target \eqn{\left(\text{ESS}^R\right)}:
#' \deqn{
#'  \bar{w}_j^i = 
#'  \frac{\sum_{l \in A_i}    \left(\text{ESS}^R-\text{ESS}_l^i\right)^q \hat{w}_{lj}^i   }{
#'  \sum_{l \in A_i} \left(\text{ESS}^R-\text{ESS}_l^i\right)^q} , \qquad q \in [0,1].
#'  }
#' If \eqn{q=0} (default), the simple average of individual weights will be calculated. 
#' If \eqn{q>0}, more weight will be assigned to locations with low ESS.
#' @examples
#' # Define simple "transmission" model where prevalence equals first parameter
#' transmission_model_identity <- function(seeds, parameters, n_tims=1) {
#'   return(matrix(parameters[,1], ncol=1))
#' }
#' 
#' # Generate samples for prevalence map with 3 locations given by B(2,1), B(1,1)=Uniform, B(1,2). 
#' set.seed(123)
#' L <- 3    # Number of locations
#' M <- 500 # Number of map samples
#' prevalence_map <- matrix(NA, L, M)
#' for (l in 1:L) {
#'   prevalence_map[l,] <- rbeta(M, max(1,l-1), max(1,3-l))
#' }
#' rownames(prevalence_map) <- c("Here","There","Somewhere else")
#' 
#' # Define 2D exponential prior
#' rprior <- function(n) {
#'   params <- matrix(NA, n, 2)
#'   colnames(params) <- c("a","b")
#'   params[,1] <- rexp(n)
#'   params[,2] <- rexp(n)
#'   return(params)
#' }
#' dprior <- function(x, log=FALSE) {
#'   if (log) {
#'     return(sum(dexp(x, log=TRUE)))
#'   } else {
#'     return(prod(dexp(x)))
#'   }
#' }
#' prior <- list(rprior=rprior,dprior=dprior)
#' 
#' # Run AMIS with default control parameters
#' amis_params <- default_amis_params()
#' output <- amis(prevalence_map, transmission_model_identity, prior, amis_params, seed=1)
#' 
#' print(output)
#' summary(output)
#' 
#' original_par <- par(no.readonly = TRUE)
#' par(cex.lab=1.5, cex.main=1.5, mar=c(5,4.5,4,2)+0.1)
#' 
#' par(mfrow=c(1,2))
#' plot_mixture_components(output, what = "uncertainty", cex=3)
#' plot_mixture_components(output, what = "density", nlevels=200)
#' 
#' par(mfrow=c(3,3))
#' plot(output, what = "a", type="hist", locations=1:L, breaks=100)
#' plot(output, what = "b", type="hist", locations=1:L, breaks=100)
#' plot(output, what = "prev", type="hist", locations=1:L, time=1, breaks=100)
#' 
#' par(mar=c(5,7.5,4,2)+0.1)
#' par(mfrow=c(1,3))
#' plot(output, what=c("a","b","prev"), type="CI", locations=1:L, ylab=NA,
#'      cex=3, lwd=3, measure_central="median", display_location_names=TRUE)
#' 
#' calculate_summaries(output, what="prev", locations=1:L, alpha=0.05)
#' 
#' # Generate new samples from the weighted posterior distributions
#' new_samples <- sample_parameters(output, n_samples = 200, locations = "Here")
#' head(new_samples)
#' plot_hist <- function(column_name){
#'   hist(new_samples[, column_name], xlab=column_name, main=paste("Histogram of", column_name))
#' }
#' par(mfrow=c(1,3))
#' plot_hist("a")
#' plot_hist("b")
#' plot_hist("prevalence")
#' 
#' par(original_par)
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params = default_amis_params(), 
                 seed = NULL, output_dir = NULL, initial_amis_vals = NULL) {

  if(!is.null(seed)){set.seed(seed)}
  
  # Checks
  checks <- check_inputs(prevalence_map, transmission_model, prior, amis_params, seed, output_dir, initial_amis_vals)
  locations_with_no_data <- checks$locs_no_data
  
  if(!is.null(initial_amis_vals)){
    message("Initialising algorithm from a previous run provided by the user.")
    amis_params <- initial_amis_vals$amis_params
    niter <- ncol(initial_amis_vals$ess_per_iteration) # number of completed iterations
    message(paste0("The previous run has completed ", niter, " iterations."))
    amis_params[["max_iters"]] <- amis_params[["max_iters"]] + niter
    message("'amis_params' used to generate the previous outputs will be used again.")
  }
  
  directory <- output_dir

  if(is.data.frame(prevalence_map)){
    prevalence_map <- as.matrix(prevalence_map)
  }
  
  save_output <- function(){
    allseeds <- 1:(niter * n_samples)
    if (is.null(rownames(prevalence_map[[1]]$data))) {
      iunames <- sapply(1:dim(weight_matrix)[2], function(idx) sprintf("iu%g", idx))
    } else {
      iunames <- rownames(prevalence_map[[1]]$data)
    }
    if (is.null(colnames(param))) {
      paramnames <- sapply(1:dim(param)[2], function(idx) sprintf("param%g", idx))
    } else {
      paramnames <- colnames(param)
    }
    if (is.null(colnames(simulated_prevalences))) {
      prevnames <- sapply(1:dim(simulated_prevalences)[2], function(idx) sprintf("prev%g", idx))
    } else {
      if(is.null(initial_amis_vals)){
        prevnames <- paste0("prev",colnames(simulated_prevalences))
      }else{
        prevnames <- colnames(simulated_prevalences)
      }
    }
    colnames(weight_matrix) <- iunames
    colnames(simulated_prevalences) <- prevnames
    names(ess) <- iunames
    rownames(ess_per_iteration) <- iunames
    colnames(ess_per_iteration) <- paste0("iter",1:niter)
    output <- list(seeds=allseeds,
                   param=param,
                   simulated_prevalences=simulated_prevalences, 
                   first_weight=first_weight,
                   weight_matrix=weight_matrix, 
                   likelihoods=likelihoods, 
                   ess=ess, 
                   prevalence_map=prevalence_map,
                   components=components, 
                   components_per_iteration=components_per_iteration,
                   ess_per_iteration=ess_per_iteration,
                   prior_density=prior_density,
                   amis_params=amis_params,
                   locs_with_g=locs_with_g, 
                   locs_without_g=locs_without_g,
                   bool_valid_sim_prev=bool_valid_sim_prev,
                   which_valid_sim_prev=which_valid_sim_prev, 
                   which_invalid_sim_prev=which_invalid_sim_prev,
                   which_valid_locs_prev_map=which_valid_locs_prev_map, 
                   locations_with_no_data=locations_with_no_data)
    class(output) <- 'amis'
    return(output)
  }
  
  # Formatting
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  likelihood_approach <- ifelse(is.null(prevalence_map[[1]]$likelihood), "nonparametric", "parametric")
  nonparametric_method <- NULL
  if(likelihood_approach=="nonparametric"){
    if(!is.null(amis_params[["breaks"]])){
      nonparametric_method <- "histogram"
    }else if(!is.null(amis_params[["sigma"]])){
      nonparametric_method <- "gaussian"
    }else{
      nonparametric_method <- "uniform"
    }
  }
  use_gaussian_kernel <- ifelse(is.null(amis_params[["breaks"]]) && !is.null(amis_params[["sigma"]]), TRUE, FALSE)
  boundaries <- amis_params[["boundaries"]]
  boundaries_param <- amis_params[["boundaries_param"]]
  n_samples <- amis_params[["n_samples"]]
  n_tims <- length(prevalence_map)
  n_locs <- nrow(prevalence_map[[1]]$data)
  M <- ncol(prevalence_map[[1]]$data)
  message("Data dimensions:")
  message(paste0("- Number of time points: ", n_tims))
  message(paste0("- Number of locations: ", n_locs))
  if(likelihood_approach=="nonparametric"){
    # Check which prevalence map samples are valid (non-NA, finite, and within boundaries)
    which_valid_prev_map <- get_which_valid_prev_map(prevalence_map, boundaries)
    message(paste0("- Number of map samples in each location: ", M))
  }else{
    which_valid_prev_map <- lapply(1:n_tims, function(t) lapply(1:n_locs, function(l) {which(!is.na(prevalence_map[[t]]$data[l,]))-1}))
  }
  message("----------------------- ")
  message("AMIS control parameters:")
  message(paste0("- Number of parameter vectors proposed at each iteration: ", n_samples))
  message(paste0("- Target effective sample size: ", amis_params[["target_ess"]]))
  message(paste0("- Maximum number of iterations: ", amis_params[["max_iters"]]))
  message("")
  message("Density estimation methods:")
  messages <- checks$messages
  if(!is.null(messages)){
    for(i in seq_along(messages)){
      message(messages[i])
    }
  }
  
  which_valid_locs_prev_map <- get_which_valid_locs_prev_map(which_valid_prev_map, n_tims, n_locs)
  # Determine at which time, for each location, denominator g will be calculated for
  locations_first_t <- get_locations_first_t(which_valid_locs_prev_map, n_tims, n_locs)
  locs_with_g <- get_locs_with_g(locations_first_t, n_tims)
  locs_without_g <- get_locs_without_g(locations_first_t, n_tims)
  # only include locations with data in locs_without_g and make NULL if empty
  locs_without_g = lapply(1:n_tims, function(t){
    wh <- which(locs_without_g[[t]]%in%which_valid_locs_prev_map[[t]])
    if(length(wh)>0){
      locs_without_g[[t]] <- locs_without_g[[t]][wh]
    }else{
      locs_without_g[[t]] <- NULL
    }
  })

  # calculate normalising constant for truncated Gaussian kernel
  log_norm_const_gaussian <- array(NA, c(1,1,1))  # for the Gaussian kernel case only, but needs to be declared
  if(use_gaussian_kernel){
    log_norm_const_gaussian <- calc_log_norm_const_gaussian(prevalence_map=prevalence_map, 
                                                            boundaries=boundaries, 
                                                            sd=amis_params[["sigma"]])
  }
  
  # Initialise
  iter <- 1
  components_per_iteration <- list()
  components_per_iteration[[1]] <- NA
  
  seeds <- function(iter) ((iter - 1) * n_samples + 1):(iter * n_samples)  # function to calculate the seeds for iteration iter.
  
  if(is.null(initial_amis_vals)){
    message("----------------------- ")
    message("AMIS iteration 1")
    message("Initialising algorithm by sampling the first set of parameters from the prior. ")
    # Sample first set of parameters from the prior
    param <- prior$rprior(n_samples)
    if(!is.matrix(param)) {stop("rprior function must produce a MATRIX of size #simulations by #parameters, even when #parameters is equal to 1.")}
    if(any(is.na(param))){warning("At least one sample from the prior was NA or NaN.")}
    if(ncol(param)==1 && amis_params[["delete_induced_prior"]]==TRUE) {warning("Currently running with amis_params[['delete_induced_prior']]=TRUE For models with only one parameter it is recommended to set amis_params[['delete_induced_prior']]=FALSE for prior to influence the weights calculation.")}
    # To avoid duplication, evaluate prior density now.
    prior_density <- sapply(1:n_samples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
    if(length(prior_density)!=n_samples) {stop("Output from dprior function must have length 1.")}
    if(any(is.na(prior_density))){warning("At least one prior density evaluation was NA or NaN.")}
    # Simulate from transmission model
    simulated_prevalences <- transmission_model(seeds = 1:n_samples, param, n_tims)
    if(!is.matrix(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, transmission_model function should produce a MATRIX of size #simulations by #timepoints, even when #timepoints is equal to 1.")}
    if(nrow(param) != nrow(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of rows in matrices from transmission_model and rprior functions must be equal (#simulations).")}
    if(length(prevalence_map) != ncol(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of timepoints in prevalence_map and the number of columns in output from transmission_model function must be equal to #timepoints.")}
    # Check validity of simulated prevalences
    bool_valid_sim_prev <- (simulated_prevalences>=boundaries[1]) & (simulated_prevalences<=boundaries[2]) & is.finite(simulated_prevalences)
    if(!is.null(boundaries_param)){
      bool_valid_sim_param <- rep(TRUE, n_samples)
      for(i_samp in 1:n_samples){
        bool_valid_sim_param[i_samp] <- all((param[i_samp,]>=boundaries_param[,1])&(param[i_samp,]<=boundaries_param[,2]))
      }
      prior_density[!bool_valid_sim_param] <- ifelse(amis_params[["log"]], -Inf, 0)
      bool_valid_sim_prev <- bool_valid_sim_prev & bool_valid_sim_param
    }
    which_valid_sim_prev <- lapply(1:n_tims, function(t) which(bool_valid_sim_prev[,t])-1L)
    which_invalid_sim_prev <- lapply(1:n_tims, function(t) which(!bool_valid_sim_prev[,t])-1L)
    # Evaluate likelihood
    likelihoods <- compute_likelihood(prevalence_map,simulated_prevalences,amis_params,
                                      likelihoods=NULL, which_valid_sim_prev,
                                      which_valid_prev_map,log_norm_const_gaussian)
    if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least 1 NaN value.")}
    # Update weight matrix
    first_weight = rep(1-amis_params[["log"]], n_samples)
    weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, amis_params,
      first_weight, locs_with_g, locs_without_g,bool_valid_sim_prev, 
      which_valid_sim_prev, which_invalid_sim_prev, which_valid_locs_prev_map, locations_with_no_data)
    if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value.")}

    ess <- calculate_ess(weight_matrix,amis_params[["log"]])

    message(paste0("  min ESS: ",round(min(ess)),", mean ESS: ",round(mean(ess)),", max ESS: ",round(max(ess))))
    message(paste0("  ",sum(ess<amis_params[["target_ess"]])," locations are below the target ESS."))
    # Make object to store the components of the AMIS mixtures of all iterations
    components <- list(
      G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
      Sigma = list(), # list of covariance matrices for each component
      Mean = list(), # list of means for each component
      probs = list() # probability of each component (unnormalised)
    )
    niter <- 1 # number of completed iterations
    ess_per_iteration <- NULL
    ess_per_iteration <- cbind(ess_per_iteration, ess)
    if (!is.null(directory)){
      res <- save_output()
      saveRDS(res, file = paste0(directory,"amis_output.rds"))
    }
  } else {
    # Check inputs of previous run provided by the user 
    check_initial_vals <- function(d){
      if(!is.null(initial_amis_vals[[d]])){
        initial_amis_vals[[d]]
      } else {
        if(d!="locations_with_no_data"){
          stop(paste0("Cannot find object 'initial_amis_vals[['",d,"']]'. To initialise from a previous run, use object saved by specifying 'output_dir'."))
        }
      }
    }

    locations_with_no_data <- check_initial_vals("locations_with_no_data")
    simulated_prevalences <- check_initial_vals("simulated_prevalences")
    likelihoods <- check_initial_vals("likelihoods")
    weight_matrix <- check_initial_vals("weight_matrix")
    ess <- check_initial_vals("ess")
    ess_per_iteration <- check_initial_vals("ess_per_iteration")
    components <- check_initial_vals("components")
    components_per_iteration <- check_initial_vals("components_per_iteration")
    param <- check_initial_vals("param")
    prior_density <- check_initial_vals("prior_density")
  }

  
  # Continue if target_ess not yet reached
  if (min(ess) >= amis_params[["target_ess"]]){
    message("----------------------- ")
    message("Algorithm finished after the first iteration with all locations below the target ESS.")
  }else{
    mixt_samples <- NULL
    mixt_samples_z <- NULL
    if(!is.null(initial_amis_vals)){
      next_iter <- niter + 1
    }else{
      next_iter <- 2
    }
    for (iter in next_iter:amis_params[["max_iters"]]) {
      message("----------------------- ")
      message("AMIS iteration ",iter)
      # Fit mixture cluster model and sample from it
      mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]], amis_params[["q"]])
      check_zero_weight_for_all_particles(amis_params, mean_weights, likelihood_approach, nonparametric_method)
      mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
      message("  A",mixture$G,"component mixture has been fitted.")
      components <- update_mixture_components(mixture, components, iter)
      new_params <- sample_new_parameters(mixture, n_samples, amis_params[["df"]], prior, amis_params[["log"]])
      components_per_iteration[[iter]] <- update_Mclust_object(mixture, new_params)
      # Check validity of sampled parameters
      if(any(is.na(new_params$params))){warning("At least one sample from the proposal after the first iteration of AMIS was NA or NaN.")}
      param <- rbind(param, new_params$params)
      if(!is.null(boundaries_param)){
        bool_valid_sim_param_iter <- rep(TRUE, n_samples)
        for(i_samp in 1:n_samples){
          bool_valid_sim_param_iter[i_samp] <- all((new_params$params[i_samp,]>=boundaries_param[,1])&(new_params$params[i_samp,]<=boundaries_param[,2]))
        }
        new_params$prior_density[!bool_valid_sim_param_iter] <- ifelse(amis_params[["log"]], -Inf, 0)
      }
      prior_density <- c(prior_density,new_params$prior_density)
      new_prevalences <- transmission_model(seeds(iter), new_params$params, n_tims)
      # Check validity of simulated prevalences
      bool_valid_sim_prev_iter <- (new_prevalences>=boundaries[1]) & (new_prevalences<=boundaries[2]) & is.finite(new_prevalences)
      if(!is.null(boundaries_param)){
        bool_valid_sim_prev_iter <- bool_valid_sim_prev_iter & bool_valid_sim_param_iter
      }
      bool_valid_sim_prev <- rbind(bool_valid_sim_prev, bool_valid_sim_prev_iter)
      which_valid_sim_prev_iter <- lapply(1:n_tims, function(t) which(bool_valid_sim_prev_iter[,t])-1L)
      which_valid_sim_prev <- lapply(1:n_tims, function(t) c(which_valid_sim_prev[[t]], which_valid_sim_prev_iter[[t]]+n_samples*(iter-1)))
      which_invalid_sim_prev_iter <- lapply(1:n_tims, function(t) which(!bool_valid_sim_prev_iter[,t])-1L)
      which_invalid_sim_prev <- lapply(1:n_tims, function(t) c(which_invalid_sim_prev[[t]], which_invalid_sim_prev_iter[[t]]+n_samples*(iter-1)))
      simulated_prevalences <- rbind(simulated_prevalences,new_prevalences)
      # Evaluate likelihood
      likelihoods <- compute_likelihood(prevalence_map,new_prevalences,amis_params, 
                                        likelihoods,which_valid_sim_prev_iter,
                                        which_valid_prev_map,log_norm_const_gaussian)
      if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least one NaN value.")}
      # Update weight matrix
      first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
      weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, 
                                             amis_params, first_weight,
                                             locs_with_g, locs_without_g,
                                             bool_valid_sim_prev, which_valid_sim_prev, 
                                             which_invalid_sim_prev, which_valid_locs_prev_map, 
                                             locations_with_no_data)
      if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value.")}
      # Calculate ESS
      ess <- calculate_ess(weight_matrix,amis_params[["log"]])
      ess_per_iteration <- cbind(ess_per_iteration, ess)
      message(paste0("  min ESS:", round(min(ess)),", mean ESS:", round(mean(ess)),", max ESS:", round(max(ess))))
      message(paste0("  ",sum(ess<amis_params[["target_ess"]])," locations are below the target ESS."))
      niter <- niter + 1
      if (!is.null(directory)){
        res <- save_output()
        saveRDS(res, file = paste0(directory,"amis_output.rds"))
      }
      if (min(ess) >= amis_params[["target_ess"]]) break
    }
  }
  message("-----------------------")

  if(niter == amis_params[["max_iters"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg, call. = FALSE)
  }
  
  # Save output
  output <- save_output()
  
  # Calculate model evidence
  model_evidence <- compute_model_evidence(likelihoods, simulated_prevalences, 
                                           amis_params, first_weight,
                                           locs_with_g, locs_without_g,
                                           bool_valid_sim_prev, which_valid_sim_prev, 
                                           which_invalid_sim_prev, which_valid_locs_prev_map, 
                                           locations_with_no_data)
  output$evidence <- model_evidence

  
  return(output)
}
