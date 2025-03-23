#ifndef QM_H
#define QM_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List QM_each(List Xorg, mat LogY, double lambda, int Xdim_max, double tau = 0.5, double h = -1.0, String kernel = "Gaussian", String penalty = "LASSO", double gamma = 0,
             double phi0 = 1e-4, double c_phi = 1.1, int max_iter = 500, double threshold = 1e-10);

double get_loss_QM(List X, mat Y, List Xnew_, mat LogYnew, double lambda, int Xdim_max, 
                   double tau = 0.5, double h = -1.0, String kernel = "Gaussian", String penalty = "LASSO", double gamma = 0);

double get_loss_CV_QM(List X_, mat Y, double lambda, int Xdim_max, String cv_type = "AIC", 
                      double tau = 0.5, double h = -1.0, String kernel = "Gaussian", String penalty = "LASSO", double gamma = 0, double cv_const = 2.0);

#endif

