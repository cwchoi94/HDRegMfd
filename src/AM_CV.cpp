// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "utils.h"
#include "SBF_utils.h"
#include "SBF_comp.h"
#include "AM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;





// [[Rcpp::export]]
List AM_CV_average(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace,
                   arma::vec lambda_list, arma::vec Xdim_max_list, arma::vec R_list, arma::mat index_mat, String cv_type = "AIC",
                   String penalty = "LASSO", double gamma = 0, int max_cv_iter = 20, double threshold = 1e-10) {

    int r1 = lambda_list.size();
    int r2 = Xdim_max_list.size();
    int r3 = R_list.size();

    double opt_lambda = lambda_list(r1 - 1);
    double opt_Xdim_max = Xdim_max_list(r2 - 1);
    double opt_R = R_list(r3 - 1);
    double opt_lambda_old = opt_lambda;
    double opt_Xdim_max_old = opt_Xdim_max;
    double opt_R_old = opt_R;

    // coordinate-wise cross validation
    List loss_list(max_cv_iter);
    mat parameter_list(3 * max_cv_iter, 3, fill::zeros);
    int iter = 0;
    int opt_idx = 0;
    vec parameters(3);
    vec loss1(r1);
    vec loss2(r2);
    vec loss3(r3);
    double loss_min = 10000;
    while (iter < max_cv_iter) {
        List loss_iter(3);

        // reduce sizes of components
        mat Xnew_r = Reduced_X_mat(Xnew, index_mat, opt_Xdim_max);
        List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);

        // lambda update
        // compute loss
        if (r1 > 1 || iter==0) {
            for (int i = 0; i < r1; i++) {
                double lambda = lambda_list(i);
                loss1(i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, lambda, opt_R, cv_type, penalty, gamma);
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
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 0) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_lambda == opt_lambda_old)) || ((r2 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }



        // Xdim_max update
        // compute loss
        if (r2 > 1) {
            for (int i = 0; i < r2; i++) {
                double Xdim_max = Xdim_max_list(i);

                // reduce sizes of components
                mat Xnew_r = Reduced_X_mat(Xnew, index_mat, Xdim_max);
                List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, Xdim_max, index_mat);

                loss2(i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, opt_R, cv_type, penalty, gamma);
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

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 1) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_Xdim_max == opt_Xdim_max_old)) || ((r1 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        // reduce sizes of components
        Xnew_r = Reduced_X_mat(Xnew, index_mat, opt_Xdim_max);
        SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);



        // R update
        // compute loss
        if (r3 > 1) {
            for (int i = 0; i < r3; i++) {
                double R = R_list(i);
                loss3(i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, R, cv_type, penalty, gamma);
            }
        }
        else {
            loss3(0) = loss_min;
        }
        
        // find opt_idx with allowing threshold error (for fast computation)
        opt_idx = get_min_idx(loss3, threshold, 0);
        opt_R = R_list(opt_idx);
        loss_min = loss3(opt_idx);
        loss_iter[2] = loss3;
        loss_list[iter] = loss_iter;

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 2) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_R == opt_R_old)) || ((r1 == 1) && (r2 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        opt_lambda_old = opt_lambda;
        opt_Xdim_max_old = opt_Xdim_max;
        opt_R_old = opt_R;

        iter = iter + 1;
    }


    List result = List::create(Named("opt.lambda") = opt_lambda, Named("opt.Xdim.max") = opt_Xdim_max, Named("opt.R") = opt_R,
        Named("loss.list") = loss_list, Named("parameter.list") = parameter_list);

    return(result);
}








// [[Rcpp::export]]
List AM_CV_integral(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace,
                    arma::vec lambda_list, arma::vec Xdim_max_list, arma::vec R_list, arma::mat index_mat, String cv_type = "AIC",
                    String penalty = "LASSO", double gamma = 0, int max_cv_iter = 20, double threshold = 1e-10) {

    int r1 = lambda_list.size();
    int r2 = Xdim_max_list.size();
    int r3 = R_list.size();

    double opt_lambda = lambda_list(r1 - 1);
    double opt_Xdim_max = Xdim_max_list(r2 - 1);
    double opt_R = R_list(r3 - 1);
    double opt_lambda_old = opt_lambda;
    double opt_Xdim_max_old = opt_Xdim_max;
    double opt_R_old = opt_R;

    // coordinate-wise cross validation
    List loss_list(max_cv_iter);
    mat parameter_list(3 * max_cv_iter, 3, fill::zeros);
    int iter = 0;
    int opt_idx = 0;
    vec parameters(3);
    vec loss1(r1);
    vec loss2(r2);
    vec loss3(r3);
    double loss_min = 10000;
    while (iter < max_cv_iter) {
        List loss_iter(3);

        // reduce sizes of components
        mat Xnew_r = Reduced_X_mat(Xnew, index_mat, opt_Xdim_max);
        List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);


        // lambda update
        // compute loss
        if (r1 > 1 || iter == 0) {
            for (int i = 0; i < r1; i++) {
                double lambda = lambda_list(i);
                loss1(i) = get_loss_CV_AM_integral(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, lambda, opt_R, cv_type, penalty, gamma);
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
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 0) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_lambda == opt_lambda_old)) || ((r2 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }


        // Xdim_max update
        // compute loss
        if (r2 > 1) {
            for (int i = 0; i < r2; i++) {
                double Xdim_max = Xdim_max_list(i);

                // reduce sizes of components
                mat Xnew_r = Reduced_X_mat(Xnew, index_mat, Xdim_max);
                List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, Xdim_max, index_mat);

                loss2(i) = get_loss_CV_AM_integral(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, opt_R, cv_type, penalty, gamma);
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

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 1) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_Xdim_max == opt_Xdim_max_old)) || ((r1 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        // reduce sizes of components
        Xnew_r = Reduced_X_mat(Xnew, index_mat, opt_Xdim_max);
        SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);



        // R update
        // compute loss
        if (r3 > 1) {
            for (int i = 0; i < r3; i++) {
                double R = R_list(i);
                loss3(i) = get_loss_CV_AM_integral(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, R, cv_type, penalty, gamma);
            }
        }
        else {
            loss3(0) = loss_min;
        }

        // find opt_idx with allowing threshold error (for fast computation)
        opt_idx = get_min_idx(loss3, threshold, 0);
        opt_R = R_list(opt_idx);
        loss_min = loss3(opt_idx);
        loss_iter[2] = loss3;
        loss_list[iter] = loss_iter;

        // parameter update
        parameters = { opt_lambda, opt_Xdim_max, opt_R };
        parameter_list.row(3 * iter + 2) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_R == opt_R_old)) || ((r1 == 1) && (r2 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        opt_lambda_old = opt_lambda;
        opt_Xdim_max_old = opt_Xdim_max;
        opt_R_old = opt_R;

        iter = iter + 1;
    }


    List result = List::create(Named("opt.lambda") = opt_lambda, Named("opt.Xdim.max") = opt_Xdim_max, Named("opt.R") = opt_R,
        Named("loss.list") = loss_list, Named("parameter.list") = parameter_list);

    return(result);
}