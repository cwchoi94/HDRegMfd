#ifndef PENALTY_H
#define PENALTY_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


double Penalty_ftn(arma::vec u, arma:: vec lambda_vec, String penalty, double gamma);
arma::vec Penalty_diff_ftn(arma::vec u, arma:: vec lambda_vec, String penalty, double gamma);

#endif

