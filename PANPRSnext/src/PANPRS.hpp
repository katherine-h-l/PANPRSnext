#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Main CPP function
//' @param summary_betas matrix
//' @param ld_J vector
//' @param num_iter_vec vector
//' @param index_matrix matrix
//' @param index_J vector
//' @param ld_vec vector
//' @param upper_val double
//' @param SD_vec matrix
//' @param tuning_matrix matrix
//' @param beta_matrix matrix
//' @param lambda0_vec vector
//' @param z_matrix matrix
//' @param all_tuning_matrix matrix
//' @param lambda_vec_func vector
//' @param func_lambda matrix
//' @param Ifunc_SNP vector
//' @param dims vector
//' @param params vector
// [[Rcpp::export]]
Rcpp::List gsfPEN_cpp(
  arma::Mat<double> summary_betas,
  arma::Col<int> ld_J,
  arma::Col<int> num_iter_vec,
  arma::Mat<int> index_matrix,
  arma::Col<int> index_J,
  arma::Col<double> ld_vec,
  double upper_val,
  arma::Mat<double> SD_vec,
  arma::Mat<double> tuning_matrix,
  arma::Mat<double> beta_matrix,
  arma::Col<double> lambda0_vec,
  arma::Mat<double> z_matrix,
  arma::Mat<double> all_tuning_matrix,
  arma::Col<double> lambda_vec_func,
  arma::Mat<int> func_lambda,
  arma::Col<int> Ifunc_SNP,
  arma::Col<int> dims,
  arma::Col<int> params
);
