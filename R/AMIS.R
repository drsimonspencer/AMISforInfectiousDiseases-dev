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
#' Each entry must be a list containing objects 
#' \describe{
#' \item{\code{data}}{An \eqn{L \times M} matrix as above}
#' \item{\code{likelihood}}{A function taking arguments 
#' \itemize{
#'    \item \code{param}: An \eqn{n \times d} matrix with the \eqn{d}-dimensional 
#'    model parameters, one for each of the \eqn{n} seeds, used in the simulation.
#'    \item \code{data}: An \eqn{L \times M} matrix as above
#'    \item \code{prevalence}: An \eqn{n}-length vector with the outputs from the transmission model.
#'    \item \code{log}: Optional logical. If set to \code{TRUE}, the function returns the log-likelihoods.
#' }
#' The function \code{likelihood} is expected to return an \eqn{L \times n} matrix.
#' }
#' }
#' The location names are inherited from \code{rownames(prevalence_map)} 
#' if \code{prevalence_map} is a matrix, or 
#' from \code{rownames(prevalence_map[[1]]$data)} if \code{prevalence_map} is a list.
#' \cr \cr If \code{likelihood} is not specified, then it is assumed that the data consist of 
#' samples from a geo-statistical model and empirical methods are used.  
#' @param transmission_model A function taking a vector of \eqn{n} seeds and 
#' an \eqn{n \times d} matrix of parameter vectors as inputs and producing 
#' an \eqn{n \times T} \bold{matrix} of prevalences as output (it must be a matrix even when \eqn{T==1}).
#' @param prior A list containing the functions \code{dprior} and \code{rprior} (density and RNG, respectively).
#' The two arguments of \code{dprior} must be:
#' \itemize{
#'   \item a \eqn{d}-length vector of model parameters; and
#'   \item a logical \code{log} to indicate whether to calculate log-density or not. 
#' }
#' The only argument of \code{rprior} must be a single integer \eqn{n} that determines the number of samples to draw. 
#' \code{rprior} must produce an \eqn{n \times d} \bold{matrix} of parameters, even when \eqn{d==1}.
#' \cr \cr Parameter names are inherited from the \code{colnames} from the output of \code{rprior} if possible.
#' @param amis_params A list containing the control parameters for the AMIS algorithm:
#' \describe{
#' \item{\code{nsamples}}{Number of new samples drawn within each AMIS iteration. Default to 500.}
#' \item{\code{boundaries}}{A vector of length two with the left and right boundaries. 
#' Default to \code{c(0,1)}. If left boundary is zero and there is no right boundary, 
#' set \code{boundaries = c(0,Inf)}.}
#' \item{\code{RN}}{Logical indicating whether to use empirical Radon-Nikodym (RN) derivative or not. 
#' Default to TRUE. If set to FALSE, likelihood is not divided by the induced prior 
#' over the simulated prevalences when calculating weights.}
#' \item{\code{delta}}{Optional smoothing parameter in the empirical RN derivative if uniform 
#' kernel (default) is used. Default to 0.01. If any of \code{sigma} or \code{breaks} is supplied, then uniform kernel will not be used.}
#' \item{\code{sigma}}{Optional smoothing parameter in the empirical RN derivative 
#' if Gaussian kernel is used (i.e. Gaussian kernel will not used). Default to NULL. 
#' If \code{breaks} is supplied, Gaussian kernel will not be used.}
#' \item{\code{breaks}}{Optional vector specifying the breaks for the histogram. 
#' Default to NULL (i.e. histogram method will not used).
#' For finite \code{boundaries}, the first and last entries of \code{breaks} must be 
#' equal to the left and right boundaries, respectively.
#' For non-finite \code{boundaries}, ensure that the range of 'breaks' includes any possible prevalence value.}
#' \item{\code{mixture_samples}}{Number of samples used to represent the weighted parameters in the mixture fitting.}
#' \item{\code{df}}{Degrees of freedom in the \eqn{t}-distributions, used to yield a heavy tailed proposal. Default to 3.}
#' \item{\code{target_ess}}{Target effective sample size. Default to 500.}
#' \item{\code{log}}{Logical indicating if calculations are to be performed on log scale. Default to FALSE.}
#' \item{\code{max_iters}}{Maximum number of AMIS iterations.}
#' \item{\code{intermittent_output}}{Optional logical indicating whether to save output to 
#' \code{\link{amis_env}} environment at each iteration of the algorithm. Default to FALSE.}
#' }
#' @param seed Optional seed for the random number generator.
#' @param initial_amis_vals Optional list containing the object created by setting \code{intermittent_output=TRUE}. 
#' Necessary to initialise the algorithm from intermittent output from a previous run (where at least one iteration was successful!).
#' @return A dataframe with the simulation seeds, sampled parameters, simulated prevalences, and weight in each location.
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params, seed = NULL, initial_amis_vals = NULL) {

  # Checks
  check_inputs(prevalence_map, transmission_model, prior, amis_params, seed)
  
  if(amis_params[["intermittent_output"]]){
    message("Saving output after each iteration (this will increase memory usage). \n")
  }

  save_output <- function(){
    
    if (is.null(initial_amis_vals)){
      allseeds <- 1:(niter * nsamples)
    } else {
      allseeds <- c(1:(initial_amis_vals$amis_bits$last_simulation_seed + (niter * nsamples)))
    }
    
    ret <- data.frame(allseeds, param, simulated_prevalences, weight_matrix)
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
    colnames(ret) <- c("seeds",paramnames,prevnames,iunames)
    
    return(list(results=ret, amis_bits = list(simulated_prevalences=simulated_prevalences, 
                                              likelihoods=likelihoods, 
                                              weight_matrix=weight_matrix, 
                                              ess=ess, 
                                              components=components, 
                                              param=param, 
                                              prior_density=prior_density,
                                              last_simulation_seed=max(allseeds))))
  }
  
  # Formatting
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  use_gaussian_kernel <- ifelse(is.null(amis_params[["breaks"]]) && !is.null(amis_params[["sigma"]]), TRUE, FALSE)
  boundaries <- amis_params[["boundaries"]]
  nsamples <- amis_params[["nsamples"]]

  # Check which prevalence map samples are valid (non-NA and within boundaries); and 
  # calculate normalising constant for truncated Gaussian kernel
  which_valid_prev_map <- get_which_valid_prev_map(prevalence_map, boundaries)
  
  log_norm_const_gaussian <- array(NA, c(1,1,1))  # for the Gaussian kernel case only, but needs to be declared
  if(use_gaussian_kernel){
    log_norm_const_gaussian <- calc_log_norm_const_gaussian(prevalence_map=prevalence_map, 
                                                            boundaries=boundaries, 
                                                            sd=amis_params[["sigma"]])
  }
  
  n_tims <- length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  which_valid_locs_prev_map <- get_which_valid_locs_prev_map(which_valid_prev_map, n_tims, n_locs)
  locations_first_t <- get_locations_first_t(which_valid_locs_prev_map, n_tims, n_locs)
  locs_RN <- get_locs_RN(locations_first_t, n_tims)
  locs_nonRN <- get_locs_nonRN(locations_first_t, n_tims)

  # Initialise
  if(!is.null(seed)) set.seed(seed)
  iter <- 1
  if(is.null(initial_amis_vals)){
    cat("AMIS iteration 1\n")
    message("Initialising algorithm by sampling the first set of parameters from the prior. \n")
    # Sample first set of parameters from the prior
    param <- prior$rprior(nsamples)
    if(!is.matrix(param)) {stop("rprior function must produce a MATRIX of size #simulations by #parameters, even when #parameters is equal to 1. \n")}
    if(any(is.na(param))){warning("At least one sample from the prior was NA or NaN. \n")}
    if(ncol(param)==1 && amis_params[["RN"]]) {warning("Currently running with amis_params[['RN']]==TRUE. 
                                                          For models with only one parameter it is recommended to set 
                                                          amis_params[['RN']] = FALSE for prior to influence the weights calculation.\n")}
    # to avoid duplication, evaluate prior density now.
    prior_density<-sapply(1:nsamples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
    if(length(prior_density)!=nsamples) {stop("Output from dprior function must have length 1. \n")}
    if(any(is.na(prior_density))){warning("At least one prior density evaluation was NA or NaN. \n")}
    # Simulate from transmission model
    simulated_prevalences <- transmission_model(seeds = 1:nsamples, param)
    if(!is.matrix(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, transmission_model function should produce a MATRIX of size #simulations by #timepoints, even when #timepoints is equal to 1. \n")}
    if(nrow(param) != nrow(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of rows in matrices from transmission_model and rprior functions must be equal (#simulations). \n")}
    if(length(prevalence_map) != ncol(simulated_prevalences)) {warning("Unless specifying a bespoke likelihood function, number of timepoints in prevalence_map and the number of columns in output from transmission_model function must be equal to #timepoints. \n")}

    is_within_boundaries <- (simulated_prevalences>=boundaries[1]) & 
      (simulated_prevalences<=boundaries[2]) &
      is.finite(simulated_prevalences)
    sim_within_boundaries <- which(is_within_boundaries)-1L
    sim_outside_boundaries <- which(!is_within_boundaries)-1L
    if(length(sim_within_boundaries)<length(simulated_prevalences)){
      warning(paste0("At iteration ",iter, ", transmission_model produced ", 
             length(simulated_prevalences)-length(sim_within_boundaries), 
             " invalid samples which will not be used.")) # invalid means -Inf, Inf, NA, NaN and values outside of boundaries
    }
    
    # to avoid duplication, evaluate likelihood now.
    likelihoods <- compute_likelihood(param,prevalence_map,simulated_prevalences,amis_params,
                                      likelihoods=NULL, sim_within_boundaries,
                                      which_valid_prev_map,log_norm_const_gaussian)
    if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least 1 NaN value. \n")}
    # Determine first time each location appears in the data


    weight_matrix <- compute_weight_matrix(
      likelihoods,
      simulated_prevalences,
      amis_params,
      first_weight = rep(1-amis_params[["log"]], nsamples),
      locs_RN,
      locs_nonRN,
      is_within_boundaries,
      sim_within_boundaries, 
      sim_outside_boundaries, 
      which_valid_locs_prev_map)
    if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value. \n")}
    
    ess <- calculate_ess(weight_matrix,amis_params[["log"]])
    cat("  min ESS:",round(min(ess))," mean ESS:",round(mean(ess))," max ESS:",round(max(ess)),"\n")
    cat(" ",length(which(ess<amis_params[["target_ess"]])),"locations are below the target ESS.\n")
    # Make object to store the components of the AMIS mixture.
    components <- list(
      G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
      Sigma = list(), # list of covariance matrices for each component
      Mean = list(), # list of means for each component
      probs = list() # probability of each component (unnormalised)
    )
    seeds <- function(iter) ((iter - 1) * nsamples + 1):(iter * nsamples)  # function to calculate the seeds for iteration iter.
    niter <- 1 # number of completed iterations 
    if (amis_params[["intermittent_output"]]){
      res <- save_output()
      assign(x="intermittent_output",value=res,envir=amis_env)
    }
  } else {
    cat("AMIS iteration 1\n")
    message("Initialising algorithm from a previous run provided by the user. \n")
    check_initial_vals = function(d){
      if(!is.null(initial_amis_vals$amis_bits[[d]])){
        initial_amis_vals$amis_bits[[d]]
      } else {
        stop(paste0("Cannot find object 'initial_amis_vals$amis_bits[['",d,"']]'. To initialise from a previous run, use object saved by setting amis_params[['intermittent output']]=TRUE.\n"))
      }
    }
    
    simulated_prevalences = check_initial_vals("simulated_prevalences")
    
    
    
    ## Do we need to check for Inf and NAs in data and simulations here as well?
    

    likelihoods <- check_initial_vals("likelihoods")
    weight_matrix <- check_initial_vals("weight_matrix")
    ess <- check_initial_vals("ess")
    components <- check_initial_vals("components")
    param <- check_initial_vals("param")
    prior_density <- check_initial_vals("prior_density")
    
    
    
    
    # locations_first_t is based on data above. Is it different when user supplies likelihood?
    
    n_tims = dim(likelihoods)[1]
    n_locs = dim(likelihoods)[2]
    if (n_locs == 1) {
      lik_matrix <- function(l) as.matrix(l)
    } else {
      lik_matrix <- function(l) t(l)
    }
    
    locations_first_t = rep(NA,n_locs)
    for (j in 1:n_tims) {
      lik_mat = lik_matrix(likelihoods[j,,])
      for (l in 1:n_locs){
        locations_first_t[!is.na(lik_mat[1,]) & is.na(locations_first_t)] = j
      }
    }
    
    seeds <- function(iter) ((iter - 2) * nsamples + 1):((iter - 1) * nsamples) + initial_amis_vals$amis_bits$last_simulation_seed #function to calculate the seeds for iteration iter.
    
    
    niter <- 0 # number of completed iterations 
    
    
  }
  
  # Checking whether some locations have no data at any time point
  if(any(locations_first_t==-1L)){
    warning(paste0(sum(locations_first_t==-1L), " location(s) provided with no data. Using prior information to determine weights for these locations.\n"))
  }
  
  # Define first_weight object in case target_ess reached in first iteration
  first_weight = rep(1-amis_params[["log"]], nsamples)
  # Continue if target_ess not yet reached
  if (min(ess) < amis_params[["target_ess"]]){
    for (iter in 2:amis_params[["max_iters"]]) {
      cat("AMIS iteration ",iter,"\n")
      mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]])
      if ((amis_params[["log"]] && max(mean_weights)==-Inf) || (!amis_params[["log"]] && max(mean_weights)==0)) {stop("No weight on any particles for locations in the active set.\n")}
      mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
      cat("  A",mixture$G,"component mixture has been fitted.\n")
      components <- update_mixture_components(mixture, components, iter)
      new_params <- sample_new_parameters(mixture, nsamples, amis_params[["df"]], prior, amis_params[["log"]])
      if(any(is.na(new_params$params))){warning("At least one sample from the proposal after the first iteration of AMIS was NA or NaN. \n")}
      param <- rbind(param, new_params$params)
      prior_density <- c(prior_density,new_params$prior_density)
      new_prevalences <- transmission_model(seeds(iter), new_params$params)

      is_within_boundaries_iter <- (new_prevalences>=boundaries[1]) & 
        (new_prevalences<=boundaries[2]) & 
        is.finite(new_prevalences)
      is_within_boundaries <- c(is_within_boundaries, is_within_boundaries_iter)
      sim_within_boundaries_iter <- which(is_within_boundaries_iter)-1L
      sim_within_boundaries <- c(sim_within_boundaries, sim_within_boundaries_iter+nsamples*(iter-1))
      sim_outside_boundaries <- c(sim_outside_boundaries, which(!is_within_boundaries_iter)-1L+nsamples*(iter-1))
      if(sum(is_within_boundaries_iter)<length(new_prevalences)){
        warning(paste0("At iteration iter=",iter, ", transmission_model produced ", 
               length(new_prevalences)-sum(is_within_boundaries_iter), 
               " invalid samples which will not be used.")) # invalid means -Inf, Inf, NA, NaN and values outside of boundaries
      }
            
      simulated_prevalences <- rbind(simulated_prevalences,new_prevalences)
      likelihoods <- compute_likelihood(new_params$params,prevalence_map,new_prevalences,amis_params, 
                                        likelihoods,sim_within_boundaries_iter,
                                        which_valid_prev_map,log_norm_const_gaussian)
      if(any(is.nan(likelihoods))) {warning("Likelihood evaluation produced at least 1 NaN value. \n")}
      first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
      weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, 
                                             amis_params, first_weight,
                                             locs_RN, locs_nonRN,
                                             is_within_boundaries, sim_within_boundaries, 
                                             sim_outside_boundaries, which_valid_locs_prev_map) # RN derivative (shd take all amis_params)
      if(any(is.na(weight_matrix))) {warning("Weight matrix contains at least one NA or NaN value. \n")}
      
      ess <- calculate_ess(weight_matrix,amis_params[["log"]])
      cat("  min ESS:",round(min(ess))," mean ESS:",round(mean(ess))," max ESS:",round(max(ess)),"\n")
      cat(" ",length(which(ess<amis_params[["target_ess"]])),"locations are below the target ESS.\n")
      niter <- niter + 1
      if (amis_params[["intermittent_output"]]){
        res <- save_output()
        assign(x="intermittent_output",value=res,envir=amis_env)
      }
      if (min(ess) >= amis_params[["target_ess"]]) break
    }
  }else{
    cat("Algorithm finished after the first iteration with all locations below the target ESS. \n")
  }

  if(niter == amis_params[["max_iters"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg)
  }
  
  # Save output
  res <- save_output()
  
  # Calculate model evidence if RN derivative not being used
  if (amis_params[["RN"]]){
    return(list(sample=res$results, evidence=NULL))
  } else {
    
    model_evidence <- NULL
    warning("model_evidence not calculated. Function compute_model_evidence() is under development.")
    # model_evidence <- compute_model_evidence(likelihoods, amis_params, first_weight)
    
    return(list(sample=res$results, evidence=model_evidence))
  }
}
