#include <RcppArmadillo.h>
#include "functions.h"

//' @title Check which prevalence samples are valid for each location at each time point.
//' @param prevalence_map List where each element is prevalence_map data for a time point.
//' @param boundaries Vector of length two.
//' @param sd Parameter sigma of the Gaussian kernel.
//' @return An (n_tims x n_locs x m) array.
//' @noRd
// [[Rcpp::export]]
arma::cube calc_log_norm_const_gaussian(const List& prevalence_map, 
                                        NumericVector boundaries, 
                                        double sd){
  int TT = prevalence_map.size();
  List prev_t_list = prevalence_map(0);
  NumericMatrix prev_t = prev_t_list["data"];
  int L = prev_t.nrow();
  int M = prev_t.ncol();
  arma::cube norm_const_gaussian = arma::zeros<arma::cube>(TT,L,M);
  LogicalVector m_idx(M);
  double p;
  for (int t=0; t<TT; t++) {
    List prev_t_list = prevalence_map(t);
    NumericMatrix prev_t = prev_t_list["data"];
    for (int l=0; l<L; l++) {
      for (int m=0; m<M; m++) {
        p = prev_t(l,m);
        m_idx[m] = !NumericVector::is_na(p) && (p>=boundaries[0]) && (p<=boundaries[1]);
        if(m_idx[m]){
          norm_const_gaussian(t, l, m) = R::pnorm(boundaries[1], p, sd, 1, 0) - 
            R::pnorm(boundaries[0], p, sd, 1, 0);
        }
      }
      double M_l = sum(m_idx);
      if(M_l>0){
        for (int m=0; m<M; m++) {
          norm_const_gaussian(t, l, m) *= sd*M_l*sqrt(2*M_PI);
        }
      }
    }
  }
  return(log(norm_const_gaussian));
}
