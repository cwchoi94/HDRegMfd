#ifndef LM_H
#define LM_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List LM_each(List Xorg, mat LogY, vec Ymu, Function inner, double lambda, int Xdim_max, double R=100,
             String penalty="LASSO", double phi = 1, double gamma=0, double eta=1e-3, int max_iter=500, double threshold=1e-10);

double get_loss_LM(List X, mat LogY, List Xnew_, mat LogYnew, vec Ymu, Function inner,
                   double lambda, int Xdim_max, double R=100, String penalty="LASSO", double phi=1, double gamma=0);

#endif

