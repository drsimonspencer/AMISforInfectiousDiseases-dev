#include <RcppArmadillo.h>
#include "functions.h"

//' @title  Compute weight matrix using empirical Radon-Nikodym derivative (Rcpp version)
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
//' @param bool_valid_sim_prev Logical vector showing which simulated values are valid.
//' @param which_valid_sim_prev Vector showing which simulated values are valid.
//' @param which_invalid_sim_prev Vector showing which simulated values are invalid.
//' @param locs Vector showing which locations have data.
//' @return An updated weight matrix.
//' @noRd
// [[Rcpp::export]]
arma::mat compute_weight_matrix_empirical_histogram(const arma::mat& likelihoods, 
                                                    const arma::vec& prev_sim, 
                                                    List amis_params,
                                                    const arma::mat& weight_matrix,
                                                    arma::uvec& bool_valid_sim_prev,
                                                    arma::uvec& which_valid_sim_prev,
                                                    arma::uvec& which_invalid_sim_prev,
                                                    arma::uvec& locs){

  int n_sims = likelihoods.n_rows;
  bool logar = amis_params["log"];
  arma::vec breaks = amis_params["breaks"];
  int B = breaks.n_elem;
  arma::uvec idx_lwr = arma::linspace<arma::uvec>(0,B-2,B-1);
  arma::uvec idx_upr = arma::linspace<arma::uvec>(1,B-1,B-1);
  arma::vec lwr = breaks(idx_lwr);
  arma::vec upr = breaks(idx_upr);
  arma::vec wdt = upr-lwr;
  double lwr_r;
  double upr_r;
  arma::ivec bin_prev_sim = arma::zeros<arma::ivec>(n_sims);

  arma::mat new_weights = weight_matrix;
  
  arma::vec kern = arma::zeros<arma::vec>(n_sims);
  arma::vec norm_const = arma::zeros<arma::vec>(n_sims);
  norm_const.fill(arma::datum::inf);
  
  if(which_valid_sim_prev.n_elem>0){

    arma::vec g_weight = arma::zeros<arma::vec>(n_sims);
    double sum_g_l = 0.0;
    double M = 0.0;
    double se = 0.0;
    arma::uvec i_row = arma::zeros<arma::uvec>(1L);
    arma::uvec l_col = arma::zeros<arma::uvec>(1L);

    // which bin each prev_sim[i] is in
    for (int b=0; b<(B-1L); b++) {
      lwr_r = lwr[b];
      upr_r = upr[b];
      arma::uvec which_in_bin_b = arma::find((prev_sim>=lwr_r) && (prev_sim<upr_r));
      if(which_in_bin_b.n_elem>0){
        (bin_prev_sim(which_in_bin_b)).fill(b);
      }
    }

    int b;
    for(auto & i : which_valid_sim_prev){
      i_row = i;
      kern.fill(0.0);
      b = bin_prev_sim[i];
      lwr_r = lwr[b];
      upr_r = upr[b];

      arma::uvec wh = arma::find(bool_valid_sim_prev && (prev_sim>=lwr_r) && (prev_sim<upr_r));
      if(wh.n_elem==0){
        Rcout << "-------- " << std::endl;
        Rcout << "A prevalence value within boundaries does not belong to any bin of the histogram: " << std::endl;
        arma::vec boundaries = amis_params["boundaries"];
        Rcout << "Boundaries: " << boundaries.t();
        Rcout << "Simulated prevalence value: " << prev_sim(i) << std::endl;
        Rcout << "First entry of breaks: " << breaks << std::endl;
        Rcout << "Last entry of breaks: " << max(breaks) << std::endl;
        Rcpp::stop("Prevalence value simulated from the transmission model is not assigned to any bin of the histogram. Ensure the range of 'breaks' includes any possible prevalence value.");
      }

      (kern(wh)).fill(1.0/wdt(b));

// Here we have the same prior for each location. ----------------
      if(logar){
        g_weight = log(kern(wh)) + weight_matrix(wh,l_col);
        M = max(g_weight);
        if(!(M == -arma::datum::inf)){
          se = sum(exp(g_weight-M));
          new_weights(i_row,locs) = weight_matrix(i_row,locs) + likelihoods(i_row,locs) - M - log(se);
        }else{
          (new_weights(i_row,locs)).fill(-arma::datum::inf);
        }
      }else{
        sum_g_l = dot(kern(wh), weight_matrix(wh,l_col));
        if(sum_g_l>0){
          new_weights(i_row,locs) = weight_matrix(i_row,locs) % likelihoods(i_row,locs)/sum_g_l;
        }else{
          (new_weights(i_row,locs)).fill(0.0);
        }
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
