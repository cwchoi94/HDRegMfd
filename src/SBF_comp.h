#ifndef SBF_COMP_H
#define SBF_COMP_H

#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;
using namespace std;
using namespace arma;


mat Reduced_X_mat(mat X, mat index_mat, int Xdim_max);
List Reduced_X_list(List X_list, mat index_mat, int Xdim_max);

List SBF_preprocessing(mat X, mat LogY, vec bandwidths, vec grids, vec weights, int degree = 0, String Kdenom_method = "numeric");
List SBF_preprocessing_reduce_dim(List SBF_comp, double Xdim_max, mat index_mat);

#endif

