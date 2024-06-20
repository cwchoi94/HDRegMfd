#ifndef PCA_LIST_H
#define PCA_LIST_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List PCA_list(List Xall, double alpha=0.95);
List predict_PCA_list(List pca, List Xnew);

#endif

