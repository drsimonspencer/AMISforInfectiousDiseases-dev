#include <RcppArmadillo.h>
#include "functions.h"

//' @title Empirical estimator for the likelihood using Uniform kernel
//' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map.
//' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
//' @param delta Bandwidth value.
//' @param which_valid_sim_prev_iter Vector showing which simulated values are valid.
//' @param which_valid_prev_map_t List showing which samples are valid for each location at a time point. 
//' @param boundaries Vector of length two.
//' @return A matrix with L rows containing the empirical estimates for the likelihood.
//' @noRd
// [[Rcpp::export]]
NumericMatrix f_estimator_uniform(NumericMatrix prevalence_map, 
                                  NumericVector prev_sim, 
                                  double delta, 
                                  NumericVector which_valid_sim_prev_iter,
                                  List& which_valid_prev_map_t,
                                  NumericVector boundaries, 
                                  bool logar){
  double left_boundary = boundaries[0];
  double right_boundary = boundaries[1];
  int R = prev_sim.length();
  int L = prevalence_map.nrow();
  int M = prevalence_map.ncol();
  double half_delta = delta/2;
  NumericMatrix f(L, R);
  NumericVector prev_sim_lo = prev_sim - half_delta;
  NumericVector prev_sim_up = prev_sim + half_delta;
  NumericVector boundaryAdj(R, 1.0);
  for(auto & r : which_valid_sim_prev_iter){
    if(prev_sim_lo[r]<left_boundary){
      boundaryAdj[r] *= delta/(delta - (left_boundary - prev_sim_lo[r]));
    }
    if(prev_sim_up[r]>right_boundary){
      boundaryAdj[r] *= delta/(delta - (prev_sim_up[r] - right_boundary));
    }
  }
  NumericVector prevalence_map_l(M);
  NumericVector prevalence_map_l_valid;
  double prev_sim_up_r;
  double prev_sim_lo_r;
  double c; 
  int M_l;
  int count;
  for (int l=0; l<L; l++) {
    IntegerVector valid_samples_t_l = which_valid_prev_map_t[l];
    M_l = valid_samples_t_l.length();
    if(M_l>0L){
      prevalence_map_l = prevalence_map(l,_);
      prevalence_map_l_valid = prevalence_map_l[valid_samples_t_l];
      c = 1.0/(delta*(double)M_l);
      for(auto & r : which_valid_sim_prev_iter){
        prev_sim_up_r = prev_sim_up[r];
        prev_sim_lo_r = prev_sim_lo[r];
        count = 0L;
        for (int m=0; m<M_l; m++) {
          if((prev_sim_up_r >= prevalence_map_l_valid[m])&&(prev_sim_lo_r <= prevalence_map_l_valid[m])){
            count++;
          }
        }
        f(l,r) = c * boundaryAdj[r] * (double)count;
      }
    }
  }
  if(logar){
    for (int l=0; l<L; l++) {
      for (int r=0; r<R; r++) {
        f(l,r) = log(f(l,r));
      }
    }
  }
  return(f);
}

// // // // 
// // Vectorisation  as line below slows down code (in both Rcpp and RcppArmadillo)
// f(l,r) = c*sum((prev_sim_up[r] >= prevalence_map_l_valid) && (prev_sim_lo[r] <= prevalence_map_l_valid));

// // Equivalent function based on arma objects is a bit slower:
// 
// arma::mat f_estimator_uniform(arma::mat& prevalence_map, 
//                               arma::vec& prev_sim, 
//                               double delta, 
//                               arma::uvec& which_valid_sim_prev_iter,
//                               List& which_valid_prev_map_t,
//                               arma::vec& boundaries, 
//                               bool logar){
//   double left_boundary = boundaries[0];
//   double right_boundary = boundaries[1];
//   int R = prev_sim.n_elem;
//   int L = prevalence_map.n_rows;
//   double half_delta = delta/2;
//   arma::mat f = arma::zeros<arma::mat>(L, R);
//   arma::vec prev_sim_lo = prev_sim - half_delta;
//   arma::vec prev_sim_up = prev_sim + half_delta;
//   arma::rowvec prevalence_map_l_valid;
//   double c; 
//   int M_l;
//   int count;
//   arma::uvec loc = arma::zeros<arma::uvec>(1L);
//   for (int l=0; l<L; l++) {
//     loc = l;
//     arma::uvec valid_samples_t_l = which_valid_prev_map_t[l];
//     M_l = valid_samples_t_l.n_elem;
//     if(M_l>0L){
//       prevalence_map_l_valid = prevalence_map(loc, valid_samples_t_l);
//       c = 1.0/(delta*(double)M_l);
//       for(auto & r : which_valid_sim_prev_iter){
//         count = 0L;
//         for (int m=0; m<M_l; m++) {
//           if((prev_sim_up[r] >= prevalence_map_l_valid[m])&&(prev_sim_lo[r] <= prevalence_map_l_valid[m])){
//             count++;
//           }
//         }
//         f(l,r) = c*(double)count;
//         if(prev_sim_lo[r]<left_boundary){
//           f(l,r) *= delta/(delta - (left_boundary - prev_sim_lo[r]));
//         }
//         if(prev_sim_up[r]>right_boundary){
//           f(l,r) *= delta/(delta - (prev_sim_up[r] - right_boundary));
//         }
//       }
//     }
//   }
//   if(logar){
//     f = log(f);
//   }
//   return(f);
// }

// // // // // // old versions
// 
// arma::mat f_estimator_uniform_Rcpp(arma::mat& prevalence_map, 
//                                    arma::vec& prev_sim, 
//                                    double delta){
//   int R = prev_sim.n_elem;
//   int L = prevalence_map.n_rows;
//   int M = prevalence_map.n_cols;
//   double half_delta = delta/2;
//   arma::mat f = arma::zeros<arma::mat>(L, R);
//   arma::vec prev_sim_lo = prev_sim - half_delta;
//   arma::vec prev_sim_up = prev_sim + half_delta;
//   arma::rowvec prevalence_map_l = arma::zeros<arma::rowvec>(M);
//   double c = 1.0/(delta*M);
//   for (int l=0; l<L; l++) {
//     prevalence_map_l = prevalence_map.row(l);
//     for (int r=0; r<R; r++) {
//       f(l,r) = sum((prev_sim_up[r] >= prevalence_map_l) && (prev_sim_lo[r] <= prevalence_map_l));
//     }
//   }
//   f = f*c;
//   return(f);
// }


// 
// int M = prevalence_map.n_cols;
// double q = 0.0;
// arma::uvec valid = arma::zeros<arma::uvec>(M);
// 
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
//     c = 1.0/(delta*(double)M_l);
//     // for (int r=0; r<R; r++) {
//     for(auto & r : which_valid_sim_prev_iter){
//       f(l,r) = c*sum((prev_sim_up[r] >= prevalence_map_l_valid) && (prev_sim_lo[r] <= prevalence_map_l_valid));
//       if(prev_sim_lo[r]<left_boundary){
//         delta_adj = delta - (left_boundary - prev_sim_lo[r]);
//         f(l,r) *= delta/delta_adj;
//       }
//       if(prev_sim_up[r]>right_boundary){
//         delta_adj = delta - (prev_sim_up[r] - right_boundary);
//         f(l,r) *= delta/delta_adj;
//       }
//     }
//   }
// }
// // // // // // finish old version

