#include <RcppArmadillo.h>
#include "functions.h"

//' @title Calculates likelihood matrix given user-defined vectorised likelihood function
//' @param likelihood_fun User-defined vectorised likelihood function
//' @param param An (n x d) matrix with the d-dimensional model parameters, 
//' one for each of the n seeds, used in the simulation.
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param sim_within_boundaries Vector showing which simulated values are within boundaries.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
NumericMatrix f_user_defined_vectorised(Rcpp::Function likelihood_fun, 
                                        NumericMatrix param,
                                        NumericMatrix prevalence_map, 
                                        NumericVector prev_sim, 
                                        arma::uvec& sim_within_boundaries, 
                                        List& which_valid_prev_map_t, 
                                        bool logar){
 
 int R = prev_sim.length();
 int L = prevalence_map.nrow();
 NumericMatrix f(L, R);
 NumericVector prevalence_map_l = prevalence_map(0,_);
 double logf;
 int M_l;
 double cmax;
 for (int l=0; l<L; l++) {
   IntegerVector valid_samples_t_l = which_valid_prev_map_t[l];
   M_l = valid_samples_t_l.length();
   if(M_l>0L){
     prevalence_map_l = prevalence_map(l, _);
     NumericVector prevalence_map_l_valid = prevalence_map_l[valid_samples_t_l];
     for(auto & r : sim_within_boundaries){
       NumericVector x = likelihood_fun(param(r,_), prevalence_map_l_valid, prev_sim[r], logar);
       if(logar){
         cmax = max(x);
         logf = cmax + log(sum(exp(x - cmax)));
         f(l,r) = logf;
       }else{
         x = log(x);
         cmax = max(x);
         logf = cmax + log(sum(exp(x - cmax)));
         f(l,r) = exp(logf);
       }
     }
   }
 }
 return(f);
}


//' @title Calculates likelihood matrix given user-defined likelihood function evaluated pointwise
//' @param likelihood_fun User-defined likelihood function evaluated pointwise
//' @param param An (n x d) matrix with the d-dimensional model parameters, 
//' one for each of the n seeds, used in the simulation.
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param sim_within_boundaries Vector showing which simulated values are within boundaries.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
arma::mat f_user_defined_pointwise(Rcpp::Function likelihood_fun, 
                                   NumericMatrix param,
                                   arma::mat& prevalence_map, 
                                   arma::vec& prev_sim, 
                                   arma::uvec& sim_within_boundaries,
                                   List& which_valid_prev_map_t, 
                                   bool logar){
  
  int R = prev_sim.n_elem;
  int L = prevalence_map.n_rows;
  arma::mat f = arma::zeros<arma::mat>(L, R);
  arma::rowvec prevalence_map_l = prevalence_map.row(0);
  double logf;
  double prev_sim_r;
  int m_i;
  int M_l;
  double cmax;
  for (int l=0; l<L; l++) {
    arma::uvec valid_samples_t_l = which_valid_prev_map_t[l];
    M_l = valid_samples_t_l.n_elem;
    if(M_l>0L){
      prevalence_map_l = prevalence_map.row(l);
      for(auto & r : sim_within_boundaries){
        prev_sim_r = prev_sim[r];
        arma::vec x = arma::zeros<arma::vec>(M_l);
        m_i = 0L;
        for(auto & m : valid_samples_t_l){
          x[m_i] = Rcpp::as<double>(likelihood_fun(param(r,_), prevalence_map_l[m], prev_sim_r, logar));
          m_i++;
        }
        if(logar){
          cmax = max(x);
          logf = cmax + log(sum(exp(x - cmax)));
          f(l,r) = logf;
        }else{
          x = log(x);
          cmax = max(x);
          logf = cmax + log(sum(exp(x - cmax)));
          f(l,r) = exp(logf);
        }
      }
    }
  }
  return(f);
}
