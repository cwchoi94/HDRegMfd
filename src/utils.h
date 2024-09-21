#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


vec SEXP_to_vec(SEXP X);
mat SEXP_to_mat(SEXP X);

double L2_norm(arma::mat Y, arma::vec Ymu, String Yspace);
double L2_norm_real(arma::mat Y);
int get_min_idx(arma::vec x, double threshold=1e-10, int type=0);

List Make_reduce_dim_matrix(List Xorg, int Xdim_max);
List Make_dimension_indices(vec Xdims_cumul);


#endif

