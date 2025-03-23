#ifndef QM_BASE_H
#define QM_BASE_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


mat sQRloss(mat u, double tau, double h, String kernel = "Gaussian");
mat sQRloss_diff(mat u, double tau, double h, String kernel = "Gaussian");

#endif

