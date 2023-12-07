#include <RcppArmadillo.h>
#include "functions.h"

//' @title Check which prevalence samples are valid for each location at each time point.
//' @param prevalence_map List where each element is prevalence_map data for a time point.
//' @param boundaries Vector of length two.
//' @return A list where each element corresponds to a time point. For each time point, there will 
//' be a list of n_locs elements, where each one shows which samples (out of M) are valid. 
//' @noRd
// [[Rcpp::export]]
List get_which_valid_prev_map(const List& prevalence_map,
                              NumericVector boundaries){
  int TT = prevalence_map.size();
  List prev_t_list = prevalence_map(0);
  NumericMatrix prev_t = prev_t_list["data"];
  int L = prev_t.nrow();
  int M = prev_t.ncol();
  LogicalVector m_idx(M);
  Rcpp::IntegerVector v = Rcpp::seq(0, M-1);
  Rcpp::List prev_map_valid_samples(TT);
  double p;
  for (int t=0; t<TT; t++) {
    List prev_t_list = prevalence_map(t);
    NumericMatrix prev_t = prev_t_list["data"];
    Rcpp::List valid_samples_t(L);
    for (int l=0; l<L; l++) {
      for (int m=0; m<M; m++) {
        p = prev_t(l,m);
        m_idx[m] = (!NumericVector::is_na(p)) &&
          (p>=boundaries[0]) && (p<=boundaries[1]) &&
          (!Rcpp::traits::is_infinite<REALSXP>(p));
      }
      IntegerVector valid_samples_t_l = v[m_idx];
      valid_samples_t[l] = valid_samples_t_l;
    }
    prev_map_valid_samples[t] = valid_samples_t;
  }
  return(prev_map_valid_samples);
}

// // // Code where -Inf and Inf values are valid
// List get_which_valid_prev_map(const List& prevalence_map, 
//                               NumericVector boundaries){
//   int TT = prevalence_map.size();
//   List prev_t_list = prevalence_map(0);
//   NumericMatrix prev_t = prev_t_list["data"];
//   int L = prev_t.nrow();
//   int M = prev_t.ncol();
//   LogicalVector m_idx(M);
//   Rcpp::IntegerVector v = Rcpp::seq(0, M-1);
//   Rcpp::List prev_map_valid_samples(TT);
//   double left_b = boundaries[0];
//   double right_b = boundaries[1];
//   bool inf_left_b = true;
//   if(!Rcpp::traits::is_infinite<REALSXP>(left_b)){
//     inf_left_b = false;
//   }
//   bool inf_right_b = true;
//   if(!Rcpp::traits::is_infinite<REALSXP>(right_b)){
//     inf_right_b = false;
//   }
//   double p;
//   for (int t=0; t<TT; t++) {
//     List prev_t_list = prevalence_map(t);
//     NumericMatrix prev_t = prev_t_list["data"];
//     Rcpp::List valid_samples_t(L);
//     for (int l=0; l<L; l++) {
//       for (int m=0; m<M; m++) {
//         p = prev_t(l,m);
//         m_idx[m] = (!NumericVector::is_na(p)) && (inf_left_b || p>=left_b) && (inf_right_b || p<=right_b);
//       }
//       IntegerVector valid_samples_t_l = v[m_idx];
//       valid_samples_t[l] = valid_samples_t_l;
//     }
//     prev_map_valid_samples[t] = valid_samples_t;
//   }
//   return(prev_map_valid_samples);
// }



//' @title Check, at each time, which locations have valid data
//' @param which_valid_prev_map List obtained by get_which_valid_prev_map
//' @param n_tims Number of time points
//' @param n_locs Number of locations
//' @noRd
// [[Rcpp::export]]
List get_which_valid_locs_prev_map(List& which_valid_prev_map, 
                                   int n_tims, int n_locs){
  List which_valid_locs_prev_map(n_tims);
  for (int t=0; t<n_tims; t++) {
    List which_valid_prev_map_t = which_valid_prev_map[t];
    arma::uvec is_valid = arma::zeros<arma::uvec>(n_locs);
    for (int l=0; l<n_locs; l++) {
      arma::uvec valid_samples_t_l = which_valid_prev_map_t[l];
      if(valid_samples_t_l.n_elem>0){
        is_valid[l] = 1L;
      }
    }
    arma::uvec which_valid_locs = arma::find(is_valid==1L);
    which_valid_locs_prev_map[t] = which_valid_locs;
  }
  return(which_valid_locs_prev_map);
}


//' @title Determine first time each location appears in the data
//' @param which_valid_locs_prev_map List obtained by get_which_valid_locs_prev_map
//' @param n_tims Number of time points
//' @param n_locs Number of locations
//' @noRd
// [[Rcpp::export]]
arma::ivec get_locations_first_t(List& which_valid_locs_prev_map, 
                                 int n_tims, int n_locs){
 int t = 0L;
 arma::ivec locations_first_t = arma::zeros<arma::ivec>(n_locs);
 locations_first_t.fill(-1L);
 for (unsigned int l=0; l<(unsigned int)n_locs; l++) {
   t = 0L;
   while((locations_first_t[l]==-1L) && (t<n_tims)){
     arma::uvec which_valid_prev_map_t = which_valid_locs_prev_map[t];
     if(any(which_valid_prev_map_t==l)){
       locations_first_t[l] = t;
     }
     t++;
   }
 }
 return(locations_first_t);
}

//' @title Determine, at which time point, which locations are updated using induced prior
//' @param locations_first_t Vector obtained by locations_first_t
//' @param n_tims Number of time points
//' @noRd
// [[Rcpp::export]]
List get_locs_with_g(arma::ivec& locations_first_t, int n_tims){
  List locs_with_g(n_tims);
  for (int t=0; t<n_tims; t++) {
    arma::uvec idx = arma::find(locations_first_t == t);
    if(idx.n_elem>0L){
      arma::uvec locs_empirical_t_arma = idx;
      locs_with_g[t] = wrap(locs_empirical_t_arma); 
    }
  }
return(locs_with_g);
}

//' @title Determine, at which time point, which locations are updated without using induced prior
//' @param locations_first_t Vector obtained by locations_first_t
//' @param n_tims Number of time points
//' @noRd
// [[Rcpp::export]]
List get_locs_without_g(arma::ivec& locations_first_t, int n_tims){
 List locs_without_g(n_tims);
 for (int t=0; t<n_tims; t++) {
   arma::uvec idx = arma::find((locations_first_t < t) && (locations_first_t != -1L));
   if(idx.n_elem>0L){
     arma::uvec locs_t_arma = idx;
     locs_without_g[t] = wrap(locs_t_arma);
   }
 }
 return(locs_without_g);
}
