#' Run the AMIS algorithm to fit a transmission model to a map
#' 
#' For details of the algorithm, see 
#' \emph{Integrating geostatistical maps and infectious disease transmission models 
#' using adaptive multiple importance sampling.}
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer. 
#' Ann. Appl. Stat. 15 (4) 1980 - 1998, December 2021.
#' DOI: \url{https://doi.org/10.1214/21-AOAS1486}
#'
#' @param prevalence_map An \eqn{L \times M} matrix containing samples from the fitted prevalence map, 
#' where \eqn{L} is the number of locations and \eqn{M} the number of samples.
#' \cr \cr Alternatively, a list with \eqn{T} elements, one for each timepoint \eqn{t=1,\dots,T}.
#' Each entry must be a list with: 
#' \describe{
#' \item{\code{data}}{An \eqn{L \times M} matrix as above}
#' \item{\code{likelihood}}{(optional) A function taking arguments:
#' \itemize{
#'    \item \code{param}: An \eqn{n \times d} matrix with the \eqn{d}-dimensional 
#'    model parameter vectors, one for each of the \eqn{n} seeds used in the simulation;
#'    \item \code{data}: A vector with a given number \eqn{M_l} of valid samples for a particular location \eqn{l};
#'    \item \code{sim_prev}: An \eqn{n}-length vector with the prevalences simulated by the transmission model;
#'    \item \code{log}: Optional logical. If set to \code{TRUE}, the function returns the log-likelihoods.
#' }
#' The function \code{likelihood} must return an \eqn{M_l \times n} matrix which will be used to 
#' calculate the likelihood of observing the simulated prevalences at a particular location.
#' }
#' }
#' The location names are inherited from \code{rownames(prevalence_map)} 
#' if \code{prevalence_map} is a matrix, and 
#' from \code{rownames(prevalence_map[[1]]$data)} if \code{prevalence_map} is a list.
#' \cr \cr If \code{likelihood} is not specified, then it is assumed that the data consist of 
#' samples from a geo-statistical model and empirical methods are used.
#' \cr \cr
#' @param transmission_model A function taking arguments:
#' \itemize{
#'    \item \code{seeds}: a vector of \eqn{n} seeds;
#'    \item \code{params}: an \eqn{n \times d} matrix of parameter vectors;
#'    \item \code{n_tims}: number of time points.
#'  }
#' This function must return an \eqn{n \times T} \bold{matrix} of prevalences 
#' (it must be a matrix even when \eqn{T=1}).
#' @param prior A list containing the functions \code{dprior} and \code{rprior} (density and RNG, respectively).
#' The two arguments of \code{dprior} must be:
#' \itemize{
#'   \item a \eqn{d}-length vector of model parameters; and
#'   \item a logical \code{log} to indicate whether to calculate log-density or not. 
#' }
#' The only argument of \code{rprior} must be a single integer \eqn{n} that determines the number of samples to draw. 
#' \code{rprior} must produce an \eqn{n \times d} \bold{matrix} of parameters, even when \eqn{d=1}.
#' Parameter names are inherited from the \code{colnames} from the output of \code{rprior} if possible.
#' @param amis_params A list containing the control parameters for the AMIS algorithm (type default_amis_params() to see default values):
#' \describe{
#' \item{\code{n_samples}}{Number of new samples drawn within each AMIS iteration. Default to 500.}
#' \item{\code{target_ess}}{Target effective sample size. Default to 500.}
#' \item{\code{max_iters}}{Maximum number of AMIS iterations. Default to 12.}
#' \item{\code{boundaries}}{A vector of length two with the left and right boundaries for prevalences. 
#' Default to \code{c(0,1)}. If left boundary is zero and there is no right boundary, 
#' set \code{boundaries = c(0,Inf)}.}
#' \item{\code{boundaries_param}}{If specified, it should be a \eqn{d \times 2} matrix 
#' with the lower and upper boundaries for the \eqn{d} transmission model parameters. Default to NULL.}
#' \item{\code{log}}{Logical indicating if calculations are to be performed on log scale. Default to TRUE}
#' \item{\code{use_induced_prior}}{Logical indicating whether the induced prior density is to be used in the update of weights. Default to TRUE.}
#' \item{\code{mixture_samples}}{Number of samples used to represent the weighted parameters in the mixture fitting.}
#' \item{\code{df}}{Degrees of freedom in the \eqn{t}-distributions, used to yield a heavy tailed proposal. Default to 3.}
#' \item{\code{delta}}{Optional smoothing parameter if uniform kernel (default) is used. Default to 0.01.}
#' \item{\code{sigma}}{Optional smoothing parameter if Gaussian kernel is used. Default to NULL.}
#' \item{\code{breaks}}{Optional vector specifying the breaks for the histogram. Default to NULL.
#' For finite \code{boundaries}, the first and last entries of \code{breaks} must be 
#' equal to the left and right boundaries, respectively.
#' For non-finite \code{boundaries}, ensure that the range of \code{breaks} includes any possible prevalence value.}
#' }
#' Uniform kernel is the default method for the density estimator of the likelihood. 
#' If \code{sigma} is provided, then Gaussian kernel will be used instead. 
#' If \code{breaks} is provided, then histogram-based method will supersede all the other methods.
#' @param seed Optional seed for the random number generator.
#' @param output_dir A string specifying the local directory where to save outputs 
#' after each iteration of the algorithm. At the end of the string, 
#' use the correct path separator for your machine's operating system. 
#' If the directory is specified, the outputs will be saved in a file called 'amis_output.rds'. 
#' Default to NULL (i.e. outputs are not saved in a local directory).
#' @param initial_amis_vals Optional list of intermittent outputs from a 
#' previous run (where at least one iteration was successful). These outputs can 
#' be saved by specifying the directory `\code{output_dir}' before running \code{\link{amis}}. 
#' @return A list of class '\code{amis}' containing:
#' \describe{
#' \item{\code{seeds}}{Vector with the simulation seeds that were used.}
#' \item{\code{param}}{A matrix with \eqn{d} columns containing the model parameters.}
#' \item{\code{simulated_prevalences}}{A matrix with \eqn{T} columns containing the simulated prevalences at each time.}
#' \item{\code{weight_matrix}}{A matrix with \eqn{L} columns containing the weights for each location.}
#' \item{\code{likelihoods}}{An array with the likelihoods for each location at each time.}
#' \item{\code{ess}}{Vector with the final ESS for each location.}
#' \item{\code{prevalence_map}}{List with the prevalence map supplied by the user.}
#' \item{\code{locations_with_no_data}}{Vector indicating which locations have no data at any time point.}
#' \item{\code{components}}{A list of the mixture components of all iterations, containing:
#'  \itemize{
#'    \item \code{G}: number of components in each iteration;
#'    \item \code{probs}: the mixture weights;
#'    \item \code{Mean}: the means of the components;
#'    \item \code{Sigma}: the covariance matrices of the components.
#'  }
#' }
#' \item{\code{components_per_iteration}}{A list with the mixture components at each iteration. 
#' This object is used in \code{\link{plot_mixture_components}}.}
#' \item{\code{ess_per_iteration}}{Matrix with with \eqn{L} rows showing the ESS for each location after each iteration.}
#' \item{\code{prior_density}}{Vector with the density function evaluated at the simulated parameter values.}
#' \item{\code{last_simulation_seed}}{Last simulation seed that was used.}
#' \item{\code{amis_params}}{List supplied by the user.}
#' }
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params = default_amis_params(), 
                 seed = NULL, output_dir = NULL, initial_amis_vals = NULL) {

  if(is.null(initial_amis_vals)){
    amis_params <- initial_amis_vals$amis_params
    cat("'amis_params' used to generate intermittent outputs will be used again.\n")
  }
  
  directory <- output_dir

  # Checks
  checks <- check_inputs(prevalence_map, transmission_model, prior, amis_params, seed, output_dir)
  locations_with_no_data <- checks$locs_no_data
  
  save_output <- function(){
    if (is.null(initial_amis_vals)){
      allseeds <- 1:(niter * n_samples)
    } else {
      allseeds <- c(1:(initial_amis_vals$last_simulation_seed + (niter * n_samples)))
    }
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
      prevnames <- paste0("prev",colnames(simulated_prevalences))
    }
    colnames(weight_matrix) <- iunames
    colnames(simulated_prevalences) <- prevnames
    names(ess) <- iunames
    rownames(ess_per_iteration) <- iunames
    colnames(ess_per_iteration) <- paste0("iter",1:niter)
    output <- list(seeds=allseeds,
                   param=param,
                   simulated_prevalences=simulated_prevalences, 
                   weight_matrix=weight_matrix, 
                   likelihoods=likelihoods, 
                   ess=ess, 
                   prevalence_map=prevalence_map,
                   locations_with_no_data=locations_with_no_data,
                   components=components, 
                   components_per_iteration=components_per_iteration,
                   ess_per_iteration=ess_per_iteration,
                   prior_density=prior_density,
                   last_simulation_seed=max(allseeds), 
                   amis_params=amis_params)
    class(output) <- 'amis'
    return(output)
  }
  
  # Formatting
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  use_gaussian_kernel <- ifelse(is.null(amis_params[["breaks"]]) && !is.null(amis_params[["sigma"]]), TRUE, FALSE)
  boundaries <- amis_params[["boundaries"]]
  boundaries_param <- amis_params[["boundaries_param"]]
  n_samples <- amis_params[["n_samples"]]
  nparams <- ncol(prior$rprior(1))
  n_tims <- length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  # Check which prevalence map samples are valid (non-NA, finite, and within boundaries)
  which_valid_prev_map <- get_which_valid_prev_map(prevalence_map, boundaries)
  which_valid_locs_prev_map <- get_which_valid_locs_prev_map(which_valid_prev_map, n_tims, n_locs)
  # Determine at which time, for each location, denominator g will be calculated for
  locations_first_t <- get_locations_first_t(which_valid_locs_prev_map, n_tims, n_locs)
  locs_with_g <- get_locs_with_g(locations_first_t, n_tims)
  locs_without_g <- get_locs_without_g(locations_first_t, n_tims)
  # calculate normalising constant for truncated Gaussian kernel
  log_norm_const_gaussian <- array(NA, c(1,1,1))  # for the Gaussian kernel case only, but needs to be declared
  if(use_gaussian_kernel){
    log_norm_const_gaussian <- calc_log_norm_const_gaussian(prevalence_map=prevalence_map, 
                                                            boundaries=boundaries, 
                                                            sd=amis_params[["sigma"]])
  }
  
  # Initialise
  if(!is.null(seed)){set.seed(seed)}
  iter <- 1
  ess_per_iteration <- NULL
  components_per_iteration <- list()
  components_per_iteration[[1]] <- NA
  
  if(is.null(initial_amis_vals)){
    cat("----------------------- \n")
    cat("AMIS iteration 1\n")
    cat("Initialising algorithm by sampling the first set of parameters from the prior. \n")
    # Sample first set of parameters from the prior
    param <- prior$rprior(n_samples)
    if(!is.matrix(param)) {stop("rprior function must produce a MATRIX of size #simulations by #parameters, even when #parameters is equal to 1. \n")}
    if(any(is.na(param))){warning("At least one sample from the prior was NA or NaN. \n")}
    if(ncol(param)==1 && amis_params[["use_induced_prior"]]==T) {warning("Currently running with amis_params[['use_induced_prior']]=TRUE. For models with only one parameter it is recommended to set amis_params[['use_induced_prior']]=FALSE for prior to influence the weights calculation.\n")}
    # To avoid duplication, evaluate prior density now.
    prior_density <- sapply(1:n_samples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
    if(length(prior_density)!=n_samples) {stop("Output from dprior function must have length 1. \n")}
    if(any(is.na(prior_density))){warning("At least one prior density evaluation was NA or NaN. \n")}
    # Simulate from transmission model
    simulated_prevalences <- transmission_model(seeds = 1:n_samples, param, n_tims)
    if(!is.matrix(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, transmission_model function should produce a MATRIX of size #simulations by #timepoints, even when #timepoints is equal to 1. \n")}
    if(nrow(param) != nrow(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of rows in matrices from transmission_model and rprior functions must be equal (#simulations). \n")}
    if(length(prevalence_map) != ncol(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of timepoints in prevalence_map and the number of columns in output from transmission_model function must be equal to #timepoints. \n")}
    # Check validity of simulated prevalences
    bool_valid_sim_prev <- (simulated_prevalences>=boundaries[1]) & (simulated_prevalences<=boundaries[2]) & is.finite(simulated_prevalences)
    if(!is.null(boundaries_param)){
      bool_valid_sim_param <- rep(T, n_samples)
      for(i_samp in 1:n_samples){
        bool_valid_sim_param[i_samp] <- all((param[i_samp,]>=boundaries_param[,1])&(param[i_samp,]<=boundaries_param[,2]))
      }
      prior_density[!bool_valid_sim_param] <- ifelse(amis_params[["log"]], -Inf, 0)
      bool_valid_sim_prev <- bool_valid_sim_prev & bool_valid_sim_param
    }
    which_valid_sim_prev <- lapply(1:n_tims, function(t) which(bool_valid_sim_prev[,t])-1L)
    which_invalid_sim_prev <- lapply(1:n_tims, function(t) which(!bool_valid_sim_prev[,t])-1L)
    # Evaluate likelihood
    likelihoods <- compute_likelihood(param,prevalence_map,simulated_prevalences,amis_params,
                                      likelihoods=NULL, which_valid_sim_prev,
                                      which_valid_prev_map,log_norm_const_gaussian)
    if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least 1 NaN value. \n")}
    # Update weight matrix
    weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, amis_params,
      first_weight = rep(1-amis_params[["log"]], n_samples), locs_with_g, locs_without_g,
      bool_valid_sim_prev, which_valid_sim_prev, which_invalid_sim_prev, which_valid_locs_prev_map, 
      locations_with_no_data)
    if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value. \n")}

    ess <- calculate_ess(weight_matrix,amis_params[["log"]])

    cat(paste0("  min ESS: ",round(min(ess)),", mean ESS: ",round(mean(ess)),", max ESS: ",round(max(ess)),"\n"))
    cat(paste0("  ",sum(ess<amis_params[["target_ess"]])," locations are below the target ESS.\n"))
    # Make object to store the components of the AMIS mixtures of all iterations
    components <- list(
      G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
      Sigma = list(), # list of covariance matrices for each component
      Mean = list(), # list of means for each component
      probs = list() # probability of each component (unnormalised)
    )
    seeds <- function(iter) ((iter - 1) * n_samples + 1):(iter * n_samples)  # function to calculate the seeds for iteration iter.
    niter <- 1 # number of completed iterations 
    if (!is.null(directory)){
      res <- save_output()
      saveRDS(res, file = paste0(directory,"amis_output.rds"))
    }
  } else {
    cat("----------------------- \n")
    cat("AMIS iteration 1\n")
    cat("Initialising algorithm from a previous run provided by the user. \n")
    # Check inputs of previous run provided by the user 
    check_initial_vals = function(d){
      if(!is.null(initial_amis_vals[[d]])){
        initial_amis_vals[[d]]
      } else {
        stop(paste0("Cannot find object 'initial_amis_vals[['",d,"']]'. To initialise from a previous run, use object saved by specifying 'output_dir'.\n"))
      }
    }
    simulated_prevalences = check_initial_vals("simulated_prevalences")
    likelihoods <- check_initial_vals("likelihoods")
    weight_matrix <- check_initial_vals("weight_matrix")
    ess <- check_initial_vals("ess")
    components <- check_initial_vals("components")
    param <- check_initial_vals("param")
    prior_density <- check_initial_vals("prior_density")
    seeds <- function(iter) ((iter - 2) * n_samples + 1):((iter - 1) * n_samples) + initial_amis_vals$last_simulation_seed #function to calculate the seeds for iteration iter.
    niter <- 0 # number of completed iterations
  }
  
  ess_per_iteration <- cbind(ess_per_iteration, ess)
  # Define first_weight object in case target_ess reached in first iteration
  first_weight = rep(1-amis_params[["log"]], n_samples)
  # Continue if target_ess not yet reached
  if (min(ess) >= amis_params[["target_ess"]]){
    cat("----------------------- \n")
    cat("Algorithm finished after the first iteration with all locations below the target ESS. \n")
  }else{
    mixt_samples <- NULL
    mixt_samples_z <- NULL
    for (iter in 2:amis_params[["max_iters"]]) {
      cat("----------------------- \n")
      cat("AMIS iteration ",iter,"\n")
      # Fit mixture cluster model and sample from it
      mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]])
      if ((amis_params[["log"]] && max(mean_weights)==-Inf) || (!amis_params[["log"]] && max(mean_weights)==0)) {stop("No weight on any particles for locations in the active set.\n")}
      mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
      cat("  A",mixture$G,"component mixture has been fitted.\n")
      components <- update_mixture_components(mixture, components, iter)
      new_params <- sample_new_parameters(mixture, n_samples, amis_params[["df"]], prior, amis_params[["log"]])
      components_per_iteration[[iter]] <- update_Mclust_object(mixture, new_params)
      # Check validity of sampled parameters
      if(any(is.na(new_params$params))){warning("At least one sample from the proposal after the first iteration of AMIS was NA or NaN. \n")}
      param <- rbind(param, new_params$params)
      if(!is.null(boundaries_param)){
        bool_valid_sim_param_iter <- rep(T, n_samples)
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
      likelihoods <- compute_likelihood(new_params$params,prevalence_map,new_prevalences,amis_params, 
                                        likelihoods,which_valid_sim_prev_iter,
                                        which_valid_prev_map,log_norm_const_gaussian)
      if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least one NaN value. \n")}
      # Update weight matrix
      first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
      weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, 
                                             amis_params, first_weight,
                                             locs_with_g, locs_without_g,
                                             bool_valid_sim_prev, which_valid_sim_prev, 
                                             which_invalid_sim_prev, which_valid_locs_prev_map, 
                                             locations_with_no_data)
      if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value. \n")}
      # Calculate ESS
      ess <- calculate_ess(weight_matrix,amis_params[["log"]])
      ess_per_iteration <- cbind(ess_per_iteration, ess)
      cat(paste0("  min ESS:", round(min(ess)),", mean ESS:", round(mean(ess)),", max ESS:", round(max(ess)),"\n"))
      cat(paste0("  ",sum(ess<amis_params[["target_ess"]])," locations are below the target ESS.\n"))
      niter <- niter + 1
      if (!is.null(directory)){
        res <- save_output()
        saveRDS(res, file = paste0(directory,"amis_output.rds"))
      }
      if (min(ess) >= amis_params[["target_ess"]]) break
    }
  }
  cat("----------------------- \n")

  if(niter == amis_params[["max_iters"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg, call. = FALSE)
  }
  
  # Save output
  output <- save_output()
  
  # Calculate model evidence only if use_induced_prior==FALSE
  if (!amis_params[["use_induced_prior"]]){
    model_evidence <- NULL
    warning("model_evidence not calculated. Function compute_model_evidence() is under development.")
    # model_evidence <- compute_model_evidence(likelihoods, amis_params, first_weight)
    output$evidence <- model_evidence
  } else {
    output$evidence <- NULL
  }
  
  return(output)
}
