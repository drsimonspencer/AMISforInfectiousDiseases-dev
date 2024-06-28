#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat compute_weight_matrix_empirical_gauss(const arma::mat& likelihoods, 
                                                const arma::vec& prev_sim, 
                                                List amis_params,
                                                const arma::mat& weight_matrix,
                                                arma::uvec& which_valid_sim_prev,
                                                arma::uvec& which_invalid_sim_prev,
                                                arma::uvec& locs);

arma::mat compute_weight_matrix_empirical_histogram(const arma::mat& likelihoods, 
                                                    const arma::vec& prev_sim, 
                                                    List amis_params,
                                                    const arma::mat& weight_matrix,
                                                    arma::uvec& bool_valid_sim_prev,
                                                    arma::uvec& which_valid_sim_prev,
                                                    arma::uvec& which_invalid_sim_prev,
                                                    arma::uvec& locs);
                                                  
arma::mat compute_weight_matrix_empirical_uniform(const arma::mat& likelihoods, 
                                                  const arma::vec& prev_sim, 
                                                  List amis_params,
                                                  const arma::mat& weight_matrix,
                                                  arma::uvec& bool_valid_sim_prev,
                                                  arma::uvec& which_valid_sim_prev,
                                                  arma::uvec& which_invalid_sim_prev,
                                                  arma::uvec& locs);

arma::mat compute_weight_matrix_without_g(const arma::mat& likelihoods, 
                                           List amis_params,
                                           const arma::mat& weight_matrix,
                                           arma::uvec& which_valid_sim_prev,
                                           arma::uvec& which_invalid_sim_prev,
                                           arma::uvec& locs);
                                      
NumericMatrix f_estimator_uniform(NumericMatrix prevalence_map, 
                                  NumericVector prev_sim, 
                                  double delta, 
                                  NumericVector which_valid_sim_prev_iter,
                                  List& which_valid_prev_map_t,
                                  NumericVector boundaries, 
                                  bool logar);
                              
arma::mat f_estimator_Gaussian(arma::mat& prevalence_map, 
                               arma::vec& prev_sim, 
                               double sd, 
                               arma::uvec& which_valid_sim_prev_iter,
                               List& which_valid_prev_map_t,
                               arma::mat& log_norm_const_gaussian_t, 
                               bool logar);
                                    
arma::mat f_estimator_histogram(arma::mat& prevalence_map, 
                                arma::vec& prev_sim, 
                                arma::vec& breaks, 
                                List& which_valid_prev_map_t, 
                                bool logar);

arma::mat f_user_defined(Rcpp::Function likelihood_fun, 
                         NumericMatrix prevalence_map,
                         List& which_valid_prev_map_t,
                         NumericVector prev_sim, 
                         NumericVector which_valid_sim_prev_iter,
                         bool logar);
                                    
List get_which_valid_prev_map(const List& prevalence_map, 
                              NumericVector boundaries);
                              
List get_which_valid_locs_prev_map(List& which_valid_prev_map, 
                                   int n_tims, int n_locs);
                                   
arma::ivec get_locations_first_t(List& which_valid_locs_prev_map, 
                                 int n_tims, int n_locs);
                                 
List get_locs_with_g(arma::ivec& locations_first_t, int n_tims);
                                 
List get_locs_without_g(arma::ivec& locations_first_t, int n_tims);

arma::cube calc_log_norm_const_gaussian(const List& prevalence_map, 
                                        NumericVector boundaries, 
                                        double sd);
                          
#endif
