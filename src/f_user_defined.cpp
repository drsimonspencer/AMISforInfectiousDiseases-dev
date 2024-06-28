#include <RcppArmadillo.h>
#include "functions.h"

//' @title Calculates an (n_locs x n_sims) matrix of likelihood values using a parametric likelihood function
//' @param likelihood_fun User-defined parametric likelihood function.
//' @param prevalence_map An (n_locs x M) matrix, where M is the number of samples from a geostatistical model 
//' or the number of likelihood parameters.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return An (n_locs x n_sims) matrix with the parametric likelihood function evaluated at the simulated prevalences.
//' @noRd
// [[Rcpp::export]]
arma::mat f_user_defined(Rcpp::Function likelihood_fun, 
                         NumericMatrix prevalence_map,
                         List& which_valid_prev_map_t,
                         NumericVector prev_sim, 
                         NumericVector which_valid_sim_prev_iter,
                         bool logar){
 int R = prev_sim.length();
 int L = prevalence_map.nrow();
 int M = prevalence_map.ncol();
 NumericVector prevalence_map_l(M);
 NumericVector prevalence_map_l_valid;
 int M_l;
 int num_valid_prev_sim;
 arma::mat f = arma::zeros<arma::mat>(L, R);
 if(logar){
   f.fill(-arma::datum::inf);
 }
 for (int l=0; l<L; l++) {
   IntegerVector valid_samples_t_l = which_valid_prev_map_t[l];
   M_l = valid_samples_t_l.length();
   if(M_l>0L){
     prevalence_map_l = prevalence_map(l, _);
     prevalence_map_l_valid = prevalence_map_l[valid_samples_t_l];
     num_valid_prev_sim = which_valid_sim_prev_iter.length();
     if(num_valid_prev_sim>0L){
       for(auto & r : which_valid_sim_prev_iter){
         f(l,r) = Rcpp::as<double>(likelihood_fun(Rcpp::_["data"] = prevalence_map_l_valid,
                                                  Rcpp::_["sim_prev"] = prev_sim[r],
                                                  Rcpp::_["log"] = logar));
       }
     }
   }
 }
 return(f);
}
