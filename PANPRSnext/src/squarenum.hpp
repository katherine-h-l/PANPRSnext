#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Squares number
//' @param i number to square
// [[Rcpp::export]]
int test_me (int i);

//' Gets eigen values
//' @param M Matrix to get eigen values for
// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M);
