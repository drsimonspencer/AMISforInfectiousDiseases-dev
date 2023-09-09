#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat compute_weight_matrix_empirical_gauss(const arma::mat& likelihoods, 
                                                const arma::vec& prev_sim, 
                                                List amis_params,
                                                const arma::mat& weight_matrix,
                                                arma::uvec& sim_within_boundaries,
                                                arma::uvec& sim_outside_boundaries,
                                                arma::uvec& locs);

arma::mat compute_weight_matrix_empirical_histogram(const arma::mat& likelihoods, 
                                                    const arma::vec& prev_sim, 
                                                    List amis_params,
                                                    const arma::mat& weight_matrix,
                                                    arma::uvec& is_within_boundaries,
                                                    arma::uvec& sim_within_boundaries,
                                                    arma::uvec& sim_outside_boundaries,
                                                    arma::uvec& locs);
                                                  
arma::mat compute_weight_matrix_empirical_uniform(const arma::mat& likelihoods, 
                                                  const arma::vec& prev_sim, 
                                                  List amis_params,
                                                  const arma::mat& weight_matrix,
                                                  arma::uvec& is_within_boundaries,
                                                  arma::uvec& sim_within_boundaries,
                                                  arma::uvec& sim_outside_boundaries,
                                                  arma::uvec& locs);
                                                  
arma::mat f_estimator_uniform(arma::mat& prevalence_map, 
                              arma::vec& prev_sim, 
                              double delta, 
                              arma::uvec& sim_within_boundaries,
                              List& which_valid_prev_map_t,
                              double left_boundary, 
                              double right_boundary);

arma::mat f_estimator_Gaussian(arma::mat& prevalence_map, 
                               arma::vec& prev_sim, 
                               double sd, 
                               arma::uvec& sim_within_boundaries,
                               List& which_valid_prev_map_t,
                               arma::mat& log_norm_const_gaussian_t,
                               double left_boundary, 
                               double right_boundary);

arma::mat f_estimator_histogram(arma::mat& prevalence_map, 
                                arma::vec& prev_sim, 
                                arma::vec& breaks, 
                                List& which_valid_prev_map_t);
                                
arma::mat compute_weight_matrix_nonRN_Rcpp(const arma::mat& likelihoods, 
                                      List amis_params,
                                      const arma::mat& weight_matrix,
                                      arma::uvec& sim_within_boundaries,
                                      arma::uvec& sim_outside_boundaries,
                                      arma::uvec& locs);
                                
List get_which_valid_prev_map(const List& prevalence_map, 
                              NumericVector boundaries);
                              
List get_which_valid_locs_prev_map(List& which_valid_prev_map, 
                                   int n_tims, int n_locs);
                                   
arma::ivec get_locations_first_t(List& which_valid_locs_prev_map, 
                                 int n_tims, int n_locs);
                                 
List get_locs_empirical(arma::ivec& locations_first_t, int n_tims);
                                 
List get_locs_bayesian(arma::ivec& locations_first_t, int n_tims);

arma::cube calc_log_norm_const_gaussian(const List& prevalence_map, 
                                        NumericVector boundaries, 
                                        double sd);
                          
#endif
