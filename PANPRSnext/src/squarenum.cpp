#include <RcppArmadillo.h>

#include "squarenum.hpp"

int test_me(int i) {
  return i * i;
}

arma::vec getEigenValues(arma::mat M) {
  return arma::eig_sym(M);
}
