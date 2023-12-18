#ifndef LR_H
#define LR_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List Make_reduce_dim_matrix(List Xorg, int Xdim_max);
List LR_each(List Xorg, mat LogY, vec Ymu, Function inner, double lambda, int Xdim_max, double R=100,
             double phi = 1, String penalty="LASSO", double gamma=0, double eta=1e-3, int max_iter=500, double threshold=1e-10);

double get_loss_LR(List X, mat LogY, List Xnew_, mat LogYnew, vec Ymu, Function inner,
                   double lambda, int Xdim_max, double R=100, double phi=1, String penalty="LASSO",double gamma=0);

#endif

