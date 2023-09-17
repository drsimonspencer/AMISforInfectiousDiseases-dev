#include <RcppArmadillo.h>
#include "functions.h"

//' @title  Compute weight matrix using empirical Radon-Nikodym derivative (Rcpp version) and Gaussian kernel
//' @description Compute matrix describing the weights for each parameter sampled, for each
//' location. One row per sample, one column per location.  Each weight 
//' is computed based on the empirical Radon-Nikodym derivative, taking into account 
//' geostatistical prevalence data for the specific location and the prevalence values 
//' computed from the transmission model for the specific parameter sample.
//' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
//' NB: transpose of slice of array. 
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
//' @param weight_matrix An n_sims x n_locs matrix containing the current values of the weights.
//' @param which_valid_sim_prev Vector showing which simulated values are valid.
//' @param which_invalid_sim_prev Vector showing which simulated values are invalid.
//' @param locs Vector showing which locations have data.
//' @return An updated weight matrix.
//' @noRd
// [[Rcpp::export]]
arma::mat compute_weight_matrix_empirical_gauss(const arma::mat& likelihoods, 
                                                const arma::vec& prev_sim, 
                                                List amis_params,
                                                const arma::mat& weight_matrix,
                                                arma::uvec& which_valid_sim_prev,
                                                arma::uvec& which_invalid_sim_prev,
                                                arma::uvec& locs){
  
  int n_sims = likelihoods.n_rows;
  // int n_locs = likelihoods.n_cols;
  bool logar = amis_params["log"];
  double sd = amis_params["sigma"];
  arma::vec boundaries = amis_params["boundaries"];
  double left_boundary = boundaries[0];
  double right_boundary = boundaries[1];
    
  arma::mat new_weights = weight_matrix;
  arma::vec norm_const = arma::zeros<arma::vec>(n_sims);
  norm_const.fill(arma::datum::inf);
  for(auto & r : which_valid_sim_prev){
      norm_const[r] = R::pnorm(right_boundary, prev_sim[r], sd, 1, 0) -
        R::pnorm(left_boundary, prev_sim[r], sd, 1, 0);
  }
  double c = 1.0/(sd*sqrt(2*M_PI));
  arma::vec kern = arma::zeros<arma::vec>(n_sims);
  arma::vec g_weight = arma::zeros<arma::vec>(n_sims);
  double sum_g_l = 0.0;
  double M = 0.0;
  double se = 0.0;
  
  if(which_valid_sim_prev.n_elem>0){
    arma::uvec i_row = arma::zeros<arma::uvec>(1L);
    arma::uvec l_col = arma::zeros<arma::uvec>(1L);
    for(auto & i : which_valid_sim_prev){
      i_row = i;
      for(auto & r : which_valid_sim_prev){
        kern[r] = exp(-0.5*pow((prev_sim[i]-prev_sim[r])/sd, 2))/norm_const[r];
      }
      kern = c*kern; 
// Here we have the same prior for each location. ----------------
      if(logar){
        g_weight = log(kern(which_valid_sim_prev)) + weight_matrix(which_valid_sim_prev,l_col);
        M = max(g_weight);
        if(!(M == -arma::datum::inf)){
          se = sum(exp(g_weight-M));
          new_weights(i_row,locs) = weight_matrix(i_row,locs) + likelihoods(i_row,locs) - M - log(se);
        }else{
          (new_weights(i_row,locs)).fill(-arma::datum::inf);
        }
      }else{
        sum_g_l = dot(kern(which_valid_sim_prev), weight_matrix(which_valid_sim_prev,l_col));
        if(sum_g_l>0){
          new_weights(i_row,locs) = weight_matrix(i_row,locs) % likelihoods(i_row,locs)/sum_g_l;
        }else{
          (new_weights(i_row,locs)).fill(0.0);
        }
      }
// --------------------------------------------------------------
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



// // If priors are different for each location, use this loop instead:
// for(auto & l : locs){
//   l_col = l;
//   if(logar){
//     g_weight = log(kern(which_valid_sim_prev)) + weight_matrix.submat(which_valid_sim_prev,l_col);
//     M = max(g_weight);
//     if(!(M == -arma::datum::inf)){
//       se = sum(exp(g_weight-M));
//       new_weights(i,l) = weight_matrix(i,l) + likelihoods(i,l) - M - log(se);
//     }else{
//       new_weights(i,l) = -arma::datum::inf;
//     }
//   }else{
//     sum_g_l = dot(kern(which_valid_sim_prev), weight_matrix.submat(which_valid_sim_prev,l_col));
//     if(sum_g_l>0){
//       new_weights(i,l) = weight_matrix(i,l)*likelihoods(i,l)/sum_g_l;
//     }else{
//       new_weights(i,l) = 0.0;
//     }
//   }
// }



// // arma::uvec locs = arma::find_finite(likelihoods.row(0L)); // this checks for NA, -Inf and Inf
// // NA does not exist in Rcpp armadillo, so I am using R_IsNA which only accepts doubles (not vec).
// for(int i=0; i<n_sims; i++){
//   not_na_lik.fill(0L);
//   for (int j=0; j<n_locs; j++) {
//     if (!R_IsNA(likelihoods(i,j))) {
//       not_na_lik[j] = 1L;
//     }
//   }
//   locs = arma::find(not_na_lik==1L);
// }

// arma::uvec not_na_lik = arma::zeros<arma::uvec>(n_locs);
// for (int l=0; l<n_locs; l++) {
//   arma::uvec valid_samples_t_l = valid_prev_map_t[l];
//   if(valid_samples_t_l.n_elem>0){
//     not_na_lik[l] = 1L;
//   }
// }
// arma::uvec locs = arma::find(not_na_lik==1L);

