#ifndef GLM_H
#define GLM_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;


List GLM_each(List Xorg, mat LogY, double lambda, int Xdim_max, double R=100, String penalty="LASSO", String link="binomial", 
              double gamma = 0, double phi = 1, double eta = 1e-3, int max_iter = 500, double threshold = 1e-10);

double get_loss_GLM(List X, mat Y, List Xnew_, mat LogYnew, double lambda, int Xdim_max, double R=100, 
                    String penalty="LASSO", String link="binomial", double gamma = 0, double phi = 1);

double get_loss_CV_GLM(List X_, mat Y, double lambda, int Xdim_max, double R = 100, String cv_type = "AIC", 
                       String penalty = "LASSO", String link = "binomial", double gamma = 0, double phi = 1, double cv_const = 2.0);

#endif

