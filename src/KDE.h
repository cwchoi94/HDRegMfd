#ifndef KDE_H
#define KDE_H

#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;
using namespace std;
using namespace arma;



List normalized_Kernel(mat X, vec bandwidths, vec grids, vec weights, int degree, String Kdenom_method = "numeric");

List KDE_(mat X, vec bandwidths, vec grids, vec weights, int degree = 0, String Kdenom_method = "numeric", bool is_proj = false);

#endif

