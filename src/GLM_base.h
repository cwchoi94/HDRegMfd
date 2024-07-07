#ifndef GLM_BASE_H
#define GLM_BASE_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


mat Link(mat u, String link = "binomial");
mat Inv_Link(mat u, String link = "binomial");
mat Psi(mat u, String link = "binomial");
mat Psi_1d(mat u, String link = "binomial");
mat Psi_2d(mat u, String link = "binomial");

#endif

