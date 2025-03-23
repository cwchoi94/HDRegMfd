// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "Penalty.h"
#include "PenaltySol.h"
#include "utils.h"
#include "PCA_list.h"
#include "QM_base.h"
#include "QM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



// [[Rcpp::export]]
List QM_CV(List X, arma::mat Y, arma::vec lambda_list, arma::vec Xdim_max_list, String cv_type = "AIC", double tau = 0.5, double h = -1.0,
           String kernel = "Gaussian", String penalty = "LASSO", double gamma = 0, double cv_const = 2.0, int max_cv_iter = 20, double threshold = 1e-10) {

    int r1 = lambda_list.size();
    int r2 = Xdim_max_list.size();

    double opt_lambda = lambda_list(r1 - 1);
    double opt_Xdim_max = Xdim_max_list(r2 - 1);
    double opt_lambda_old = opt_lambda;
    double opt_Xdim_max_old = opt_Xdim_max;

    // coordinate-wise cross validation
    List loss_list(max_cv_iter);
    mat parameter_list(2 * max_cv_iter, 2, fill::zeros);
    int iter = 0;
    int opt_idx = 0;
    vec parameters(2);
    vec loss1(r1);
    vec loss2(r2);
    double loss_min = 10000;
    while (iter < max_cv_iter) {
        List loss_iter(2);

        // lambda update
        // compute loss
        if (r1 > 1 || iter==0) {
            for (int i = 0; i < r1; i++) {
                double lambda = lambda_list(i);
                loss1(i) = get_loss_CV_QM(X, Y, lambda, opt_Xdim_max, cv_type, tau, h, kernel, penalty, gamma, cv_const);
            }
        }
        else {
            loss1(0) = loss_min;
        }
        
        // find opt_idx with allowing threshold error (for fast computation)
        opt_idx = get_min_idx(loss1, threshold, 0);
        opt_lambda = lambda_list(opt_idx);
        loss_min = loss1(opt_idx);
        loss_iter[0] = loss1;

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max };
        parameter_list.row(2 * iter + 0) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_lambda == opt_lambda_old)) || (r2 == 1)) {
            loss_list(iter) = loss_iter;
            break;
        }


        // Xdim_max update
        // compute loss
        if (r2 > 1) {
            for (int i = 0; i < r2; i++) {
                double Xdim_max = Xdim_max_list(i);
                loss2(i) = get_loss_CV_QM(X, Y, opt_lambda, Xdim_max, cv_type, tau, h, kernel, penalty, gamma, cv_const);
            }
        }
        else {
            loss2(0) = loss_min;
        }

        // find opt_idx with allowing threshold error (for fast computation)
        opt_idx = get_min_idx(loss2, threshold, 1);
        opt_Xdim_max = Xdim_max_list(opt_idx);
        loss_min = loss2(opt_idx);
        loss_iter[1] = loss2;        
        loss_list[iter] = loss_iter;

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max };
        parameter_list.row(2 * iter + 1) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_Xdim_max == opt_Xdim_max_old)) || (r1 == 1)) {
            loss_list(iter) = loss_iter;
            break;
        }

        opt_lambda_old = opt_lambda;
        opt_Xdim_max_old = opt_Xdim_max;

        iter = iter + 1;
    }


    List result = List::create(Named("opt.lambda") = opt_lambda, Named("opt.Xdim.max") = opt_Xdim_max, 
        Named("loss.list") = loss_list, Named("parameter.list") = parameter_list);

    return(result);
}