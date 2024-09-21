#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;


vec inner(mat u, mat v, vec p, String space);


#endif

