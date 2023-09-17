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
arma::mat compute_weight_matrix_empirical_uniform(const arma::mat& likelihoods, 
                                                  const arma::vec& prev_sim, 
                                                  List amis_params,
                                                  const arma::mat& weight_matrix,
                                                  arma::uvec& bool_valid_sim_prev,
                                                  arma::uvec& which_valid_sim_prev,
                                                  arma::uvec& which_invalid_sim_prev,
                                                  arma::uvec& locs){

  int n_sims = likelihoods.n_rows;
  bool logar = amis_params["log"];
  double delta = amis_params["delta"];
  double half_delta = delta/2;
  arma::vec prev_sim_lo = prev_sim - half_delta;
  arma::vec prev_sim_up = prev_sim + half_delta;
  arma::vec boundaries = amis_params["boundaries"];
  double left_boundary = boundaries[0];
  double right_boundary = boundaries[1];
  
  arma::mat new_weights = weight_matrix;
  
  arma::vec kern = arma::zeros<arma::vec>(n_sims);
  arma::vec norm_const = arma::zeros<arma::vec>(n_sims);
  norm_const.fill(arma::datum::inf);
  
  if(which_valid_sim_prev.n_elem>0){
    
    for(auto & r : which_valid_sim_prev){
        norm_const[r] = delta;
        if(prev_sim_lo[r]<left_boundary){
          norm_const[r] = delta - (left_boundary - prev_sim_lo[r]);
        }
        if(prev_sim_up[r]>right_boundary){
          norm_const[r] = delta - (prev_sim_up[r] - right_boundary);
        }
    }
    arma::vec g_weight = arma::zeros<arma::vec>(n_sims);
    double sum_g_l = 0.0;
    double M = 0.0;
    double se = 0.0;
    arma::uvec i_row = arma::zeros<arma::uvec>(1L);
    arma::uvec l_col = arma::zeros<arma::uvec>(1L);
    for(auto & i : which_valid_sim_prev){
      i_row = i;
      kern.fill(0.0);
      // arma::uvec wh = arma::find(abs(prev_sim-prev_sim[i]) <= half_delta);
      arma::uvec wh = arma::find(bool_valid_sim_prev && (abs(prev_sim-prev_sim[i])<=half_delta));
      kern(wh) = 1.0/norm_const(wh);
      
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
// --------------------------------------------------------------
// Rcout << "se: " << se << std::endl;
// Rcout << "sum_g_l: " << sum_g_l << std::endl;
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

//     for(auto & l : locs){
//       l_col = l;
//       if(logar){
//         g_weight = log(kern(wh)) + weight_matrix.submat(wh,l_col);
//         M = max(g_weight);
//         if(!(M == -arma::datum::inf)){
//           se = sum(exp(g_weight-M));
//           new_weights(i,l) = weight_matrix(i,l) + likelihoods(i,l) - M - log(se);
//         }else{
//           new_weights(i,l) = -arma::datum::inf;
//         }
//       }else{
//         sum_g_l = dot(kern(wh), weight_matrix.submat(wh,l_col));
//         if(sum_g_l>0){
//           new_weights(i,l) = weight_matrix(i,l)*likelihoods(i,l)/sum_g_l;
//         }else{
//           new_weights(i,l) = 0.0;
//         }
//       }
//     }
//




// // ############## non-vectorised version   

// for(int r=0; r<n_sims; r++){
//   if(fabs(prev_sim[r]-prev_sim[i]) <= half_delta){
//     kern[r] = 1.0/norm_const[r];
//   }
// }

// if(logar){
//   g_weight = log(kern) + weight_matrix.col(l);
//   M = max(g_weight);
//   if(!(M == -arma::datum::inf)){
//     se = 0.0;
//     for(int i=0; i<n_sims; i++){
//       se += exp(g_weight[i]-M);
//     }
//     new_weights(i,l) = weight_matrix(i,l) + likelihoods(i,l) - M - log(se);
//   }else{
//     new_weights(i,l) = -arma::datum::inf;
//   }
// }else{
//   g_weight = kern % weight_matrix.col(l);
//   sum_g_l = sum(g_weight);
//   if(sum_g_l>0){
//     new_weights(i,l) = weight_matrix(i,l)*likelihoods(i,l)/sum_g_l;
//   }else{
//     new_weights(i,l) = 0.0;
//   }
// }
// // ##############

// No need to use below code again because which_valid_sim_prev is used
// arma::vec kern = arma::zeros<arma::vec>(n_sims);
// if((prev_sim[i] >= left_boundary) && (prev_sim[i] <= right_boundary)){
//   kern = c;
// }else{
//   kern.fill(0.0);
// }


// arma::vec norm_const = arma::zeros<arma::vec>(n_sims);
// norm_const.fill(arma::datum::inf);
// for(int r=0; r<n_sims; r++){
//   if((prev_sim[r] >= left_boundary) && (prev_sim[r] <= right_boundary)){
//     norm_const[r] = delta;
//     if(prev_sim_lo[r]<left_boundary){
//       norm_const[r] = delta - (left_boundary - prev_sim_lo[r]);
//     }
//     if(prev_sim_up[r]>right_boundary){
//       norm_const[r] = delta - (prev_sim_up[r] - right_boundary);
//     }
//   }
// }

// if((prev_sim[i] >= left_boundary) && (prev_sim[i] <= right_boundary)){
//   for(int r=0; r<n_sims; r++){
//     kern[r] = c[r];
//   }
// }else{
//   kern.fill(0.0);
// }


// // // // // // old
// arma::mat compute_weight_matrix_empirical_Rcpp(const arma::mat& likelihoods, 
//                                                const arma::vec& prev_sim, 
//                                                List amis_params,
//                                                const arma::mat& weight_matrix,
//                                                arma::uvec& which_valid_sim_prev,
//                                                arma::uvec& which_invalid_sim_prev,
//                                                arma::uvec& locs){
//   
//   int n_sims = likelihoods.n_rows;
//   // int n_locs = likelihoods.n_cols;
//   
//   double delta = amis_params["delta"];
//   bool logar = amis_params["log"];
//   
//   arma::uvec i_row = arma::zeros<arma::uvec>(1L);
//   arma::mat new_weights = weight_matrix;
//   
//   for(int i=0; i<n_sims; i++){
//     
//     int n_locs_updated = locs.n_elem;
//     
//     i_row = i;
//     arma::uvec wh = arma::find(abs(prev_sim-prev_sim[i]) <= delta/2);
//     //  wh.n_elem is always >0 because of prev_sim[i] itself!!
//     arma::uvec which_nonfinite;
//     arma::mat g_terms = weight_matrix.submat(wh, locs);
//     arma::uvec which_finite;
//     arma::uvec non_zero_locs;
//     
//     if(logar){
//       arma::vec M = arma::zeros<arma::vec>(n_locs_updated);
//       arma::vec log_col_sums = arma::zeros<arma::vec>(n_locs_updated);
//       for(int j=0; j<n_locs_updated; j++){
//         M[j] = max(g_terms.col(j));
//         log_col_sums[j] = log(sum(exp(g_terms.col(j) - M[j])));
//       }
//       which_nonfinite = arma::find_nonfinite(M); // is it possible to have +Inf in the prior? If so, then this needs to change
//       which_finite = arma::find_finite(M);
//       non_zero_locs = locs(which_finite);
//       new_weights(i_row, non_zero_locs) = weight_matrix(i_row, non_zero_locs) + 
//         likelihoods(i_row, non_zero_locs) - M(which_finite).t() - log_col_sums(which_finite).t() + log(delta);
//     }else{
//       arma::vec col_sums = arma::zeros<arma::vec>(n_locs_updated);
//       for(int j=0; j<n_locs_updated; j++){
//         col_sums[j] = sum(g_terms.col(j));
//       }
//       which_nonfinite = arma::find(col_sums==0);
//       which_finite = arma::find(col_sums>0);
//       non_zero_locs = locs(which_finite);
//       new_weights(i_row, non_zero_locs) = delta * weight_matrix(i_row, non_zero_locs) % likelihoods(i_row, non_zero_locs) / (col_sums(which_finite).t());
//     }
// 
//     if(which_nonfinite.n_elem>0){
//       arma::uvec zero_locs = locs(which_nonfinite);
//       arma::vec value_for_zero_locs = arma::zeros<arma::vec>(zero_locs.n_elem);
//       if(logar){
//         value_for_zero_locs.fill(log(0.0));
//       }
//       new_weights(i_row, zero_locs) = value_for_zero_locs;
//     }
//     
//   }
//   
//   return(new_weights);
// }
// 
// // not_na_lik.fill(0L);
// // for (int j=0; j<n_locs; j++) {
// //   if (!R_IsNA(likelihoods(i,j))) {
// //     not_na_lik[j] = 1L;
// //   }
// // }
// // arma::uvec locs = arma::find(not_na_lik==1L);
// 
// // arma::uvec not_na_lik = arma::zeros<arma::uvec>(n_locs);
// // for (int l=0; l<n_locs; l++) {
// //   arma::uvec valid_samples_t_l = valid_prev_map_t[l];
// //   if(valid_samples_t_l.n_elem>0){
// //     not_na_lik[l] = 1L;
// //   }
// // }
// // arma::uvec locs = arma::find(not_na_lik==1L);
