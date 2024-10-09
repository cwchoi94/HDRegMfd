#ifndef AM_H
#define AM_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List AM_each(List kernel_list, vec Ymu, String Yspace, double lambda, double R=100, String penalty="LASSO", double gamma=0, 
             double phi = 1, double eta=1e-3, int max_iter=500, double threshold=1e-10);


double get_loss_CV_AM(List SBF_comp, mat Xnew, mat LogYnew, vec Ymu, String Yspace, double lambda, double R = 100, String cv_type = "AIC", String penalty = "LASSO", double gamma = 0);
double get_loss_CV_AM2(List SBF_comp, mat Xnew, mat LogYnew, vec Ymu, String Yspace, double lambda, double R = 100, String cv_type = "AIC", String penalty = "LASSO", double gamma = 0);


#endif

