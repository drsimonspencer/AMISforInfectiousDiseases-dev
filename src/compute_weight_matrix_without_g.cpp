#include <RcppArmadillo.h>
#include "functions.h"

//' @title  Compute weight matrix without using induced prior
//' @description Compute matrix describing the weights for each parameter sampled, for each
//' location. One row per sample, one column per location.  Each weight 
//' is computed based on the empirical Radon-Nikodym derivative, taking into account 
//' geostatistical prevalence data for the specific location and the prevalence values 
//' computed from the transmission model for the specific parameter sample.
//' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
//' NB: transpose of slice of array. 
//' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
//' @param weight_matrix An n_sims x n_locs matrix containing the current values of the weights.
//' @param which_valid_sim_prev Vector showing which simulated values are valid.
//' @param which_invalid_sim_prev Vector showing which simulated values are invalid.
//' @param locs Vector showing which locations have data.
//' @return An updated weight matrix.
//' @noRd
// [[Rcpp::export]]
arma::mat compute_weight_matrix_without_g(const arma::mat& likelihoods, 
                                           List amis_params,
                                           const arma::mat& weight_matrix,
                                           arma::uvec& which_valid_sim_prev,
                                           arma::uvec& which_invalid_sim_prev,
                                           arma::uvec& locs){

  bool logar = amis_params["log"];
  arma::mat new_weights = weight_matrix;

  if(which_valid_sim_prev.n_elem>0){
    arma::uvec i_row = arma::zeros<arma::uvec>(1L);
    for(auto & i : which_valid_sim_prev){
      i_row = i;
      if(logar){
        new_weights(i_row,locs) = weight_matrix(i_row,locs) + likelihoods(i_row,locs);
      }else{
        new_weights(i_row,locs) = weight_matrix(i_row,locs) % likelihoods(i_row,locs);
      }
    }
  }
  
  if(which_invalid_sim_prev.n_elem>0){
    if(logar){
      (new_weights(which_invalid_sim_prev,locs)).fill(-arma::datum::inf);
    }else{
      (new_weights(which_invalid_sim_prev,locs)).fill(0.0);
    }
  }
  
  return(new_weights);
}
