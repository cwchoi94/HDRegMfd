#ifndef SBF_UTILS_H
#define SBF_UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


int multi_4d_ind_to_single(int j1, int j2, int k1, int k2, int p, int g);
IntegerVector single_ind_to_multi_4d(int j3, int p, int g);
IntegerVector multi_2d_ind_to_single_range(int j1, int j2, int p, int g);
IntegerVector multi_3d_ind_to_single_range(int j1, int j2, int k1, int p, int g);

double numerical_integral_1d(vec weights, vec func);
mat numerical_integral_3d_to_2d(vec weights, cube func);
cube numerical_integral_3d(vec weights, cube func1, cube func2);
mat numerical_integral_2d(mat func1, mat func2);




#endif

