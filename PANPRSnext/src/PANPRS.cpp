#include <RcppArmadillo.h>

#include "PANPRS.hpp"

int test_me(int i) {
  return i * i;
}

Rcpp::List gsfPEN_cpp(
  arma::mat summary_betas, // matrix
  arma::vec ld_J, // not matrix
  arma::vec num_iter_vec, // not matrix
  arma::mat index_matrix, // matrix
  arma::vec index_J, // not matrix
  arma::vec ld_vec, // not matrix
  arma::vec upper_val, // not matrix
  arma::mat SD_vec, // matrix
  arma::mat tuning_matrix,// matrix
  arma::mat beta_matrix,// matrix
  arma::vec lambda0_vec, // not matrix
  arma::mat z_matrix, // matrix
  arma::mat all_tuning_matrix, // matrix
  arma::vec lambda_vec_func, // not matrix
  arma::mat func_lambda, // matrix
  arma::vec Ifunc_SNP, // not matrix
  arma::vec dims, // not matrix
  arma::vec params // not matrix
) {
  arma::mat Matrix1;
  arma::vec Vector1;
  arma::mat Matrix2;

  dims(0) = 99;

  dims.raw_print();

  return Rcpp::List::create(
    Rcpp::Named("BetaMatrix") = Matrix1,
    Rcpp::Named("IterationVector") = Vector1,
    Rcpp::Named("AllTuningMatrix") = Matrix2
  );
}
