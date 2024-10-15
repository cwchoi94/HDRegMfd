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
List AM_CBS_kfold(List X_list, List LogY_list, List Xnew_list, List LogYnew_list, List Ymu_list, String Yspace, int kfold, arma::mat bandwidths_mat, arma::vec grids, arma::vec weights,
                  arma::vec lambda_list, arma::vec Xdim_max_list, arma::vec R_list, arma::mat index_mat, String cv_type = "AIC",
                  String penalty = "LASSO", double gamma = 0, int degree = 0, String Kdenom_method = "numeric", int max_cv_iter = 20, double threshold = 1e-10) {

    int r1 = lambda_list.size();
    int r2 = Xdim_max_list.size();
    int r3 = R_list.size();
    int r0 = bandwidths_mat.n_rows;
    int p = bandwidths_mat.n_cols;

    double opt_lambda = lambda_list(r1 - 1);
    double opt_Xdim_max = Xdim_max_list(r2 - 1);
    double opt_R = R_list(r3 - 1);
    double opt_lambda_old = opt_lambda;
    double opt_Xdim_max_old = opt_Xdim_max;
    double opt_R_old = opt_R;
    vec opt_bandwidths = bandwidths_mat.row(r0 - 1).t();
    vec opt_bandwidths_old = opt_bandwidths;


    // coordinate-wise cross validation
    List loss_list(max_cv_iter);
    mat parameter_list(4 * max_cv_iter, p + 3, fill::zeros);
    int iter = 0;
    int opt_idx = 0;
    vec parameters(p+3);
    vec tmp_parameters(3);
    mat loss0(r0, p);
    vec loss1(r1);
    vec loss2(r2);
    vec loss3(r3);
    mat loss0_j_mat(kfold, r0);
    mat loss1_mat(kfold, r1);
    mat loss2_mat(kfold, r2);
    mat loss3_mat(kfold, r3);
    double loss_min = 10000;
    while (iter < max_cv_iter) {
        List loss_iter(4);

        // compute X_r
        List Xnew_r_list = Reduced_X_list(Xnew_list, index_mat, opt_Xdim_max);

        // bandwidths update
        if (r0 > 1 || iter == 0) {
            for (int j = 0; j < p; j++) {
                // compute loss
                vec tmp_bandwidths = opt_bandwidths;
                for (int idx = 0; idx < kfold; idx++) {
                    mat X = X_list[idx];
                    mat LogY = LogY_list[idx];
                    mat Xnew_r = Xnew_r_list[idx];
                    mat LogYnew = LogYnew_list[idx];
                    vec Ymu = Ymu_list[idx];

                    for (int i = 0; i < r0; i++) {
                        // compute SBF_comp and SBF_comp_r
                        tmp_bandwidths(j) = bandwidths_mat(i, j);

                        List SBF_comp = SBF_preprocessing(X, LogY, tmp_bandwidths, grids, weights, degree, Kdenom_method);
                        List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);

                        loss0_j_mat(idx, i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, opt_R, cv_type, penalty, gamma);
                    }
                }

                // find opt_idx with allowing threshold error (for fast computation)
                vec loss0_j = arma::conv_to<vec>::from(mean(loss0_j_mat, 0));
                opt_idx = get_min_idx(loss0_j, threshold, 0);
                opt_bandwidths(j) = bandwidths_mat(opt_idx, j);
                loss_min = loss0_j(opt_idx);
                loss0.col(j) = loss0_j;
            }
        }
        else {
            for (int i = 0; i < r0; i++) {
                for (int j = 0; j < p; j++) {
                    loss0(i, j) = loss_min;
                }
            }
        }
        
        loss_iter[0] = loss0;

        // parameter update
        tmp_parameters = { opt_lambda,opt_Xdim_max,opt_R };
        parameters = arma::join_cols(opt_bandwidths, tmp_parameters);
        parameter_list.row(4 * iter + 0) = trans(parameters);

        // check convergence
        if (((iter >= 1) && arma::all(opt_bandwidths == opt_bandwidths_old))) {
            loss_list(iter) = loss_iter;
            break;
        }


        // compute SBF_comp and SBF_comp_r
        List SBF_comp_list(kfold);
        List SBF_comp_r_list(kfold);
        for (int idx = 0; idx < kfold; idx++) {
            mat X = X_list[idx];
            mat LogY = LogY_list[idx];

            List SBF_comp = SBF_preprocessing(X, LogY, opt_bandwidths, grids, weights, degree, Kdenom_method);
            List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);

            SBF_comp_list[idx] = SBF_comp;
            SBF_comp_r_list[idx] = SBF_comp_r;
        }


        // lambda update
        // compute loss
        if (r1 > 1) {
            for (int idx = 0; idx < kfold; idx++) {
                List SBF_comp_r = SBF_comp_r_list[idx];
                mat Xnew_r = Xnew_r_list[idx];
                mat LogYnew = LogYnew_list[idx];
                vec Ymu = Ymu_list[idx];
                for (int i = 0; i < r1; i++) {
                    double lambda = lambda_list(i);
                    loss1_mat(idx, i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, lambda, opt_R, cv_type, penalty, gamma);
                }
            }            
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss1_mat(idx, 0) = loss_min;
            }
        }
        
        // find opt_idx with allowing threshold error (for fast computation)
        loss1 = arma::conv_to<vec>::from(mean(loss1_mat, 0));
        opt_idx = get_min_idx(loss1, threshold, 0);
        opt_lambda = lambda_list(opt_idx);
        loss_min = loss1(opt_idx);
        loss_iter[1] = loss1;

        // parameter update
        tmp_parameters = { opt_lambda,opt_Xdim_max,opt_R };
        parameters = arma::join_cols(opt_bandwidths, tmp_parameters);
        parameter_list.row(4 * iter + 1) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_lambda == opt_lambda_old)) || ((r0 == 1) && (r2 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }



        // Xdim_max update
        // compute loss
        if (r2 > 1) {
            for (int idx = 0; idx < kfold; idx++) {
                List SBF_comp = SBF_comp_list[idx];
                mat Xnew = Xnew_r_list[idx];
                mat LogYnew = LogYnew_list[idx];
                vec Ymu = Ymu_list[idx];
                for (int i = 0; i < r2; i++) {
                    double Xdim_max = Xdim_max_list(i);

                    // reduce sizes of components
                    mat Xnew_r = Reduced_X_mat(Xnew, index_mat, Xdim_max);
                    List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, Xdim_max, index_mat);

                    loss2_mat(idx, i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, opt_R, cv_type, penalty, gamma);
                }
            }
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss2_mat(idx, 0) = loss_min;
            }
        }

        // find opt_idx with allowing threshold error (for fast computation)
        loss2 = arma::conv_to<vec>::from(mean(loss2_mat, 0));
        opt_idx = get_min_idx(loss2, threshold, 1);
        opt_Xdim_max = Xdim_max_list(opt_idx);
        loss_min = loss2(opt_idx);
        loss_iter[2] = loss2;        

        // parameter update
        tmp_parameters = { opt_lambda,opt_Xdim_max,opt_R };
        parameters = arma::join_cols(opt_bandwidths, tmp_parameters);
        parameter_list.row(4 * iter + 2) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_Xdim_max == opt_Xdim_max_old)) || ((r0 == 1) && (r1 == 1) && (r3 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        // reduce sizes of components
        Xnew_r_list = Reduced_X_list(Xnew_list, index_mat, opt_Xdim_max);
        for (int idx = 0; idx < kfold; idx++) {
            List SBF_comp = SBF_comp_list[idx];
            List SBF_comp_r = SBF_preprocessing_reduce_dim(SBF_comp, opt_Xdim_max, index_mat);

            SBF_comp_r_list[idx] = SBF_comp_r;
        }



        // R update
        // compute loss
        if (r3 > 1) {
            for (int idx = 0; idx < kfold; idx++) {
                List SBF_comp_r = SBF_comp_r_list[idx];
                mat Xnew_r = Xnew_r_list[idx];
                mat LogYnew = LogYnew_list[idx];
                vec Ymu = Ymu_list[idx];
                for (int i = 0; i < r3; i++) {
                    double R = R_list(i);
                    loss3_mat(idx, i) = get_loss_CV_AM_average(SBF_comp_r, Xnew_r, LogYnew, Ymu, Yspace, opt_lambda, R, cv_type, penalty, gamma);
                }
            }
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss3_mat(idx, 0) = loss_min;
            }
        }
        
        // find opt_idx with allowing threshold error (for fast computation)
        loss3 = arma::conv_to<vec>::from(mean(loss3_mat, 0));
        opt_idx = get_min_idx(loss3, threshold, 0);
        opt_R = R_list(opt_idx);
        loss_min = loss3(opt_idx);
        loss_iter[3] = loss3;
        loss_list[iter] = loss_iter;

        // parameter update
        tmp_parameters = { opt_lambda,opt_Xdim_max,opt_R };
        parameters = arma::join_cols(opt_bandwidths, tmp_parameters);
        parameter_list.row(4 * iter + 3) = trans(parameters);

        // check convergence
        if (((iter >= 1) && (opt_R == opt_R_old)) || ((r0 == 1) && (r1 == 1) && (r2 == 1))) {
            loss_list(iter) = loss_iter;
            break;
        }

        opt_bandwidths_old = opt_bandwidths;
        opt_lambda_old = opt_lambda;
        opt_Xdim_max_old = opt_Xdim_max;
        opt_R_old = opt_R;

        iter = iter + 1;
    }


    List result = List::create(Named("opt.bandwidths") = opt_bandwidths, Named("opt.lambda") = opt_lambda, Named("opt.Xdim.max") = opt_Xdim_max, Named("opt.R") = opt_R,
        Named("loss.list") = loss_list, Named("parameter.list") = parameter_list);

    return(result);
}



