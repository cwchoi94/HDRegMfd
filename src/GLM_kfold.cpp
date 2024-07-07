// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "PenaltySol.h"
#include "utils.h"
#include "PCA_list.h"
#include "GLM_base.h"
#include "GLM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;




// [[Rcpp::export]]
List GLM_Kfold(List X_list, List Y_list, List Xnew_list, List Ynew_list, int kfold, arma::vec lambda_list, arma::vec Xdim_max_list, arma::vec R_list, 
              String penalty="LASSO", String link="binomial", double phi = 1, double gamma = 0, int max_cv_iter = 20, double threshold = 1e-10) {
  
    int r1 = lambda_list.size();
    int r2 = Xdim_max_list.size();
    int r3 = R_list.size();
        
    // CBS cross validation  
    double opt_lambda = lambda_list(r1-1);
    double opt_Xdim_max = Xdim_max_list(r2-1);
    double opt_R = R_list(r3-1);
    double opt_lambda_old = opt_lambda;
    double opt_Xdim_max_old = opt_Xdim_max;
    double opt_R_old = opt_R;
  
    List loss_list(max_cv_iter);
    mat parameter_list(3*max_cv_iter,3,fill::zeros);
    int iter = 0;
    int opt_idx = 0;
    vec parameters(3);
    vec loss1(r1);
    vec loss2(r2);
    vec loss3(r3);
    mat loss1_mat(kfold,r1);
    mat loss2_mat(kfold,r2);
    mat loss3_mat(kfold,r3);
    vec loss_min(kfold);
    while (iter<max_cv_iter){
        List loss_iter(3);
    
        // lambda update
        // compute loss
        if (r1 > 1 || iter == 0) {
            for (int idx = 0; idx < kfold; idx++) {
                List X = X_list[idx];
                mat Y = Y_list[idx];
                List Xnew = Xnew_list[idx];
                mat Ynew = Ynew_list[idx];
                for (int i = 0; i < r1; i++) {
                    double lambda = lambda_list(i);
                    loss1_mat(idx, i) = get_loss_GLM(X, Y, Xnew, Ynew, lambda, opt_Xdim_max, opt_R, penalty, link, phi, gamma);
                }
            }
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss1_mat(idx, 0) = loss_min(idx);
            }            
        }
        
    
        // find opt_idx with allowing threshold error (for fast computation)
        loss1 = arma::conv_to<vec>::from(mean(loss1_mat,0));
        opt_idx = get_min_idx(loss1, threshold, 0);
        opt_lambda = lambda_list(opt_idx);
        loss_min = loss1_mat.col(opt_idx);
        loss_iter[0] = loss1;
    
        // parameter update
        parameters = {opt_lambda, opt_Xdim_max, opt_R};
        parameter_list.row(3*iter+0) = trans(parameters);
    
        // check convergence
        if ((iter>=1) & (opt_lambda==opt_lambda_old)){
            loss_list(iter) = loss_iter;
            break;
        }
    
        // Xdim_max update
        // compute loss
        if (r2 > 1) {
            for (int idx = 0; idx < kfold; idx++) {
                List X = X_list[idx];
                mat Y = Y_list[idx];
                List Xnew = Xnew_list[idx];
                mat Ynew = Ynew_list[idx];
                for (int i = 0; i < r2; i++) {
                    double Xdim_max = Xdim_max_list(i);
                    loss2_mat(idx, i) = get_loss_GLM(X, Y, Xnew, Ynew, opt_lambda, Xdim_max, opt_R, penalty, link, phi, gamma);
                }
            }
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss2_mat(idx, 0) = loss_min(idx);
            }
        }
        
    
        // find opt_idx with allowing threshold error (for fast computation)
        loss2 = arma::conv_to<vec>::from(mean(loss2_mat,0));
        opt_idx = get_min_idx(loss2, threshold, 1);
        opt_Xdim_max = Xdim_max_list(opt_idx);
        loss_min = loss2_mat.col(opt_idx);
        loss_iter[1] = loss2;
    
        // parameter update
        parameters = {opt_lambda, opt_Xdim_max, opt_R};
        parameter_list.row(3*iter+1) = trans(parameters);
    
        // check convergence
        if ((iter>=1) & (opt_Xdim_max==opt_Xdim_max_old)){
            loss_list(iter) = loss_iter;
            break;
        }
    
        // R update
        // compute loss
        if (r3 > 1) {
            for (int idx = 0; idx < kfold; idx++) {
                List X = X_list[idx];
                mat Y = Y_list[idx];
                List Xnew = Xnew_list[idx];
                mat Ynew = Ynew_list[idx];
                for (int i = 0; i < r3; i++) {
                    double R = R_list(i);
                    loss3_mat(idx, i) = get_loss_GLM(X, Y, Xnew, Ynew, opt_lambda, opt_Xdim_max, R, penalty, link, phi, gamma);
                }
            }
        }
        else {
            for (int idx = 0; idx < kfold; idx++) {
                loss3_mat(idx, 0) = loss_min(idx);
            }
        }
        
    
        // find opt_idx with allowing threshold error (for fast computation)
        loss3 = arma::conv_to<vec>::from(mean(loss3_mat,0));
        opt_idx = get_min_idx(loss3, threshold, 0);
        opt_R = R_list(opt_idx);
        loss_min = loss3_mat.col(opt_idx);
        loss_iter[2] = loss3;
        loss_list[iter] = loss_iter;
    
        // parameter update
        parameters = {opt_lambda, opt_Xdim_max, opt_R};
        parameter_list.row(3*iter+2) = trans(parameters);
    
        // check convergence
        if ((iter>=1) & (opt_R==opt_R_old)){
            loss_list(iter) = loss_iter;
            break;
        }
    
        opt_lambda_old = opt_lambda;
        opt_Xdim_max_old = opt_Xdim_max;
        opt_R_old = opt_R;
    
        iter = iter+1;
    }

    List result = List::create(Named("opt.lambda") = opt_lambda, Named("opt.Xdim.max") = opt_Xdim_max, Named("opt.R") = opt_R,
                               Named("loss.list") = loss_list, Named("parameter.list") = parameter_list);
    
    return(result);
}


