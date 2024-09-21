#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;


arma::vec inner(arma::mat u, arma::mat v, arma::vec p, String space);


#endif

