#include <RcppArmadillo.h>
#include "functions.h"

//' @title Empirical estimator for the likelihood using Gaussian kernel
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param sd Bandwith value.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param log_norm_const_gaussian_t (n_locs x M) matrix showing the log normalising constant for the Gaussian kernels.
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
arma::mat f_estimator_Gaussian(arma::mat& prevalence_map, 
                               arma::vec& prev_sim, 
                               double sd, 
                               arma::uvec& which_valid_sim_prev_iter,
                               List& which_valid_prev_map_t,
                               arma::mat& log_norm_const_gaussian_t, 
                               bool logar){
  int R = prev_sim.n_elem;
  int L = prevalence_map.n_rows;
  arma::mat f = arma::zeros<arma::mat>(L, R);
  arma::rowvec prevalence_map_l = prevalence_map.row(0);
  double front = -0.5/pow(sd,2);
  double logf;
  double prev_sim_r;
  arma::uvec loc = arma::zeros<arma::uvec>(1L);
  int m_i;
  int M_l;
  double cmax;
  for (int l=0; l<L; l++) {
    loc = l;
    arma::uvec valid_samples_t_l = which_valid_prev_map_t[l];
    M_l = valid_samples_t_l.n_elem;
    if(M_l>0L){
      arma::vec log_norm_const_gaussian_t_valid = (log_norm_const_gaussian_t(loc, valid_samples_t_l)).t();
      prevalence_map_l = prevalence_map.row(l);
      for(auto & r : which_valid_sim_prev_iter){
        prev_sim_r = prev_sim[r];
        arma::vec x = arma::zeros<arma::vec>(M_l);
        m_i = 0L;
        for(auto & m : valid_samples_t_l){
          x[m_i] = front*pow((prev_sim_r-prevalence_map_l[m]), 2);
          m_i++;
        }
        x -= log_norm_const_gaussian_t_valid;
        cmax = max(x);
        logf = cmax + log(sum(exp(x - cmax)));
        if(logar){
          f(l,r) = logf;
        }else{
          f(l,r) = exp(logf);  
        }
      }
    }
  }
  return(f);
}

// // // Rcpp version without arma is actually slower:
// // [[Rcpp::export]]
// NumericMatrix f_estimator_Gaussian(NumericMatrix prevalence_map, 
//                                     NumericVector prev_sim, 
//                                     double sd, 
//                                     NumericVector which_valid_sim_prev_iter,
//                                     List& which_valid_prev_map_t,
//                                     NumericMatrix log_norm_const_gaussian_t, 
//                                     bool logar){
//   int R = prev_sim.length();
//   int L = prevalence_map.nrow();
//   int M = prevalence_map.ncol();
//   NumericMatrix f(L, R);
//   NumericVector prevalence_map_l(M);
//   int int2 = 2L;
//   double front = -0.5/pow(sd,int2);
//   double logf;
//   double prev_sim_r;
//   NumericVector log_norm_const_gaussian_t_l(M);
//   int m_i;
//   int M_l;
//   double cmax;
//   for (int l=0; l<L; l++) {
//     IntegerVector valid_samples_t_l = which_valid_prev_map_t[l];
//     M_l = valid_samples_t_l.length();
//     if(M_l>0L){
//       log_norm_const_gaussian_t_l = log_norm_const_gaussian_t(l,_);
//       NumericVector log_norm_const_gaussian_t_valid = log_norm_const_gaussian_t_l[valid_samples_t_l];
//       prevalence_map_l = prevalence_map(l,_);
//       for(auto & r : which_valid_sim_prev_iter){
//         prev_sim_r = prev_sim[r];
//         NumericVector x(M_l);
//         m_i = 0L;
//         for(auto & m : valid_samples_t_l){
//           x[m_i] = front*pow((prev_sim_r-prevalence_map_l[m]), int2);
//           m_i++;
//         }
//         // NumericVector prev_sim_r_vec = rep(prev_sim_r, M_l);
//         // NumericVector prevalence_map_l_valid = prevalence_map_l[valid_samples_t_l];
//         // NumericVector x = front*pow((prev_sim_r_vec - prevalence_map_l_valid), int2);
//         x = x - log_norm_const_gaussian_t_valid;
//         cmax = max(x);
//         logf = cmax + log(sum(exp(x - cmax)));
//         if(logar){
//           f(l,r) = logf;
//         }else{
//           f(l,r) = exp(logf);  
//         }
//       }
//     }
//   }
//   return(f);
// }


// Other ways to use norm_const_gauss (NOT in log scale, so log must not be taken in calc_log_norm_const_gaussian)

// arma::rowvec norm_const_gaussian_t_l = norm_const_gaussian_t.row(l);

// // NOT using LSE trick ------------------------------
// res = 0.0;
// for(auto & m : valid_samples_t_l){
//   res += exp(front*pow((prev_sim_r-prevalence_map_l[m]), 2))/norm_const_gaussian_t_l[m];
// }

// // using LSE trick - vectorised version -------------
// arma::vec x = front*pow((prev_sim_r-prevalence_map_l(valid_samples_t_l)), 2) - log(norm_const_gaussian_t_l(valid_samples_t_l));
// double cmax = max(x);
// res = exp(cmax + log(sum(exp(x - cmax))));

// // using LSE trick - non-vectorised version ---------
// arma::vec x = arma::zeros<arma::mat>(valid_samples_t_l.n_elem);
// int m_i = 0L;
// for(auto & m : valid_samples_t_l){
//   x[m_i] = front*pow((prev_sim_r-prevalence_map_l[m]), 2) - log(norm_const_gaussian_t_l[m]);
//   m_i++;
// }
// double cmax = max(x);
// res = exp(cmax + log(sum(exp(x - cmax))));


// // // // // // old versions

// arma::mat f_estimator_Gaussian(arma::mat prevalence_map, 
//                                arma::vec prev_sim, 
//                                double sd, 
//                                double left_boundary, 
//                                double right_boundary){
//   int R = prev_sim.n_elem;
//   int L = prevalence_map.n_rows;
//   int M = prevalence_map.n_cols;
//   arma::mat f = arma::zeros<arma::mat>(L, R);
//   arma::vec norm_const = arma::zeros<arma::vec>(M);
//   double c = 1.0/(sd*M*sqrt(2*M_PI));
//   for (int l=0; l<L; l++) {
//     for (int m=0; m<M; m++) {
//       norm_const[m] = R::pnorm(right_boundary, prevalence_map(l,m), sd, 1, 0) -
//         R::pnorm(left_boundary, prevalence_map(l,m), sd, 1, 0);
//     }
//     for (int r=0; r<R; r++) {
//       for (int m=0; m<M; m++) {
//         f(l,r) += exp(-0.5*pow((prev_sim[r]-prevalence_map(l,m))/sd, 2))/norm_const[m];
//       }
//     }
//   }
//   f = c*f;
//   return(f);
// }


// int M = prevalence_map.n_cols;
// if(which_valid_sim_prev_iter.n_elem>0){ // c++ indexing
//   which_valid_sim_prev_iter -= 1L;
// }

// int M_l = 1L;
// double q = 0.0;
// arma::uvec valid = arma::zeros<arma::uvec>(M);
// arma::uvec loc = arma::zeros<arma::uvec>(1L);

// for (int l=0; l<L; l++) {
//   loc = l;
//   valid.fill(0L);
//   for (int m=0; m<M; m++) {
//     q = prevalence_map(l,m);
//     if(!NumericVector::is_na(q) && (q>=left_boundary) && (q<=right_boundary)){
//       valid[m] = 1L;
//     }
//   }
//   arma::uvec idx = arma::find(valid == 1L);
//   M_l = idx.n_elem;
//   if(M_l>0L){
//     prevalence_map_l_valid = prevalence_map(loc, idx);
//     arma::vec norm_const = arma::zeros<arma::vec>(M_l);
//     for (int mv=0; mv<M_l; mv++) {
//       norm_const[mv] = R::pnorm(right_boundary, prevalence_map_l_valid[mv], sd, 1, 0) -
//         R::pnorm(left_boundary, prevalence_map_l_valid[mv], sd, 1, 0);
//     }
//     norm_const *= sd*(double)M_l*sqrt(2*M_PI);
//     // for (int r=0; r<R; r++) {
//     for(auto & r : which_valid_sim_prev_iter){
//       for (int mv=0; mv<M_l; mv++) {
//         f(l,r) += exp(-0.5*pow((prev_sim[r]-prevalence_map_l_valid[mv])/sd, 2))/norm_const[mv];
//       }
//     }
//   }
// }
// // // // // // finish old version

