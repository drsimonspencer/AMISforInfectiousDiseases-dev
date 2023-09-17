#include <RcppArmadillo.h>
#include "functions.h"

//' @title Empirical estimator for the likelihood using Uniform kernel
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param breaks vector specifying the breaks for the histogram.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
arma::mat f_estimator_histogram(arma::mat& prevalence_map, 
                                arma::vec& prev_sim, 
                                arma::vec& breaks, 
                                List& which_valid_prev_map_t, 
                                bool logar){
  int R = prev_sim.n_elem;
  int L = prevalence_map.n_rows;
  int B = breaks.n_elem;
  arma::uvec idx_lwr = arma::linspace<arma::uvec>(0,B-2,B-1);
  arma::uvec idx_upr = arma::linspace<arma::uvec>(1,B-1,B-1);
  arma::vec lwr = breaks(idx_lwr);
  arma::vec upr = breaks(idx_upr);
  arma::vec wdt = upr-lwr;
  arma::mat f = arma::zeros<arma::mat>(L, R);
  arma::rowvec prevalence_map_l_valid;
  double c; 
  int M_l;
  int count;
  double lwr_r;
  double upr_r;
  arma::uvec loc = arma::zeros<arma::uvec>(1L);
  for (int l=0; l<L; l++) {
    loc = l;
    arma::uvec valid_samples_t_l = which_valid_prev_map_t[l];
    M_l = valid_samples_t_l.n_elem;
    if(M_l>0L){
      prevalence_map_l_valid = prevalence_map(loc, valid_samples_t_l);
      c = 1.0/((double)M_l);
      for (int b=0; b<(B-1L); b++) {
        lwr_r = lwr[b];
        upr_r = upr[b];
        arma::uvec wh = arma::find((prev_sim>=lwr_r) && (prev_sim<upr_r));
        if(wh.n_elem>0L){
          count = 0L;
          for (int m=0; m<M_l; m++) {
            if((prevalence_map_l_valid[m]>=lwr_r)&&(prevalence_map_l_valid[m]<upr_r)){
              count++;
            }
          }
          (f(loc, wh)).fill(c*(double)count/(wdt[b]));
        }
      }
    }
  }
  if(logar){
    f = log(f);
  }
  return(f);
}

// for(auto & r : which_valid_sim_prev_iter){
//   bin = max(arma::find(prev_sim>=lwr));   // this is wrong
//   lwr_r = lwr[bin];
//   upr_r = upr[bin];
//   count = 0L;
//   for (int m=0; m<M_l; m++) {
//     if((prevalence_map_l_valid[m]>=lwr_r)&&(prevalence_map_l_valid[m]<upr_r)){
//       count++;
//     }
//   }
//   f(l,r) = c*(double)count/(wdt[bin]);
// }
