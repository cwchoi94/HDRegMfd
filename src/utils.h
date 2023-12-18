#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


vec SEXP_to_vec(SEXP X);
mat SEXP_to_mat(SEXP X);

double L2_norm(arma::mat Y, arma::vec Ymu, Function inner);
int get_min_idx(arma::vec x, double threshold=1e-10, int type=0);

#endif

