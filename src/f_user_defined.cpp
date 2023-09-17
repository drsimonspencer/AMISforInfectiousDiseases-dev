#include <RcppArmadillo.h>
#include "functions.h"

//' @title Calculates likelihood matrix given user-defined likelihood function (1st version)
//' @param likelihood_fun User-defined likelihood function evaluated at simulation r at location l 
//' (i.e., it is vectorised for the M valid samples and r valid simulations and 
//' returns an M x R matrix)
//' @param param An (n x d) matrix with the d-dimensional model parameters, 
//' one for each of the n seeds, used in the simulation.
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return An (n_locs x n_sims) matrix with the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
NumericMatrix f_user_defined_l(Rcpp::Function likelihood_fun, 
                               NumericMatrix param,
                               NumericMatrix prevalence_map, 
                               NumericVector prev_sim, 
                               IntegerVector which_valid_sim_prev_iter, 
                               List& which_valid_prev_map_t, 
                               bool logar){
 int R = prev_sim.length();
 int L = prevalence_map.nrow();
 NumericMatrix f(L, R);
 NumericVector prevalence_map_l = prevalence_map(0,_);
 double logf;
 int M_l;
 int R_l = which_valid_sim_prev_iter.length();
 int r_valid_idx;
 double cmax;
 if(R_l>0){
   NumericVector prev_sim_valid = prev_sim[which_valid_sim_prev_iter];
   for (int l=0; l<L; l++) {
     IntegerVector valid_samples_t_l = which_valid_prev_map_t[l];
     M_l = valid_samples_t_l.length();
     if(M_l>0L){
       prevalence_map_l = prevalence_map(l, _);
       NumericVector prevalence_map_l_valid = prevalence_map_l[valid_samples_t_l];
       NumericMatrix xMat = likelihood_fun(param, prevalence_map_l_valid, prev_sim_valid, logar);
       for (int r=0; r<R_l; r++) {
         r_valid_idx = which_valid_sim_prev_iter[r];
         NumericVector x = xMat(_,r);
         if(logar){
           cmax = max(x);
           logf = cmax + log(sum(exp(x - cmax)));
           f(l, r_valid_idx) = logf;
         }else{
           x = log(x);
           cmax = max(x);
           logf = cmax + log(sum(exp(x - cmax)));
           f(l, r_valid_idx) = exp(logf);
         }
       }
     }
   }
 }
 return(f);
}


//' @title Calculates likelihood matrix given user-defined likelihood function (2nd version)
//' @param likelihood_fun User-defined likelihood function evaluated at simulation r at location l 
//' (i.e., it is vectorised for the M valid samples and returns an M-length vector)
//' @param param An (n x d) matrix with the d-dimensional model parameters, 
//' one for each of the n seeds, used in the simulation.
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return An (n_locs x n_sims) matrix with the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
NumericMatrix f_user_defined_l_r(Rcpp::Function likelihood_fun, 
                                 NumericMatrix param,
                                 NumericMatrix prevalence_map, 
                                 NumericVector prev_sim, 
                                 arma::uvec& which_valid_sim_prev_iter, 
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
     for(auto & r : which_valid_sim_prev_iter){
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


//' @title Calculates likelihood matrix given user-defined likelihood function (3rd version)
//' @param likelihood_fun User-defined likelihood function evaluated at simulation r for sample m at location l 
//' @param param An (n x d) matrix with the d-dimensional model parameters, 
//' one for each of the n seeds, used in the simulation.
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param logar Logical indicating if the outputs should be in log-scale or not
//' @return An (n_locs x n_sims) matrix with the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
arma::mat f_user_defined_l_m_r(Rcpp::Function likelihood_fun, 
                               NumericMatrix param,
                               arma::mat& prevalence_map, 
                               arma::vec& prev_sim, 
                               arma::uvec& which_valid_sim_prev_iter,
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
      for(auto & r : which_valid_sim_prev_iter){
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
