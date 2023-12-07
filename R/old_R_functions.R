# These are original R functions which were replaced by Rcpp versions

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
compute_weight_matrix_no_induced_prior <- function(likelihoods, amis_params, weight_matrix) {
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  if (amis_params[["log"]]) {
    weight_matrix[,locs]<-weight_matrix[,locs]+likelihoods[,locs]
  } else {
    weight_matrix[,locs]<-weight_matrix[,locs]*likelihoods[,locs]
  }
  return(weight_matrix)
}
