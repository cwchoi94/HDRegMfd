// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "Penalty.h"
#include "PenaltySol.h"
#include "utils.h"
#include "QM_base.h"
#include "QM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



// [[Rcpp::export]]
List QM_each(List Xorg, arma::mat Yorg, double lambda, int Xdim_max, double tau, double h, String kernel, String penalty, double gamma, double phi0, double c_phi, int max_iter, double threshold){

    int p = Xorg.size();
    int n = Yorg.n_rows;
    int m = Yorg.n_cols;
    
    // Join each Xj in Xorg
    List Xlist = Make_reduce_dim_matrix(Xorg,Xdim_max);
    mat X = Xlist["X"];
    vec Xdims = Xlist["Xdims"];
    int P = X.n_cols;

    // Define basic parameters
    vec Xdims_cumul_tmp = cumsum(Xdims);
    int nXdims_tmp = Xdims_cumul_tmp.size();
    vec Xdims_cumul(nXdims_tmp+1); // add 0 at first in Xdims_cumul
    for (int i=1; i<=nXdims_tmp; i++){
        Xdims_cumul(i) = Xdims_cumul_tmp(i-1);
    }
    vec Xdims_sqrt = sqrt(Xdims);

    // Compute a basic h if h<0 (default)
    if (h < 0.0) {
        h = max(0.05, tau * (1.0 - tau) * pow(log(sum(Xdims)) / n, 1.0 / 4.0));
    }

    // compute indices
    // ind1: indices involving beta_j
    // ind2: indices involving beta_k for k!=j
    List ind_list = Make_dimension_indices(Xdims_cumul);
    List ind_list1 = ind_list["ind1"];
    List ind_list2 = ind_list["ind2"];
        
    // Preprocessing and compute XX and XY
    vec W(P);
    for (int i=0; i<P; i++){
        W(i) = mean(square(X.col(i)));
    }
    uvec W_nonzero = find(W != 0);    
    int n_W_nonzero = W_nonzero.size();
    for (int i=0; i<n_W_nonzero; i++){
        int j = W_nonzero(i);
        X.col(j) = X.col(j) / sqrt(W(j));
    }
  
    mat XX = trans(X) * X / n;
    mat XY = trans(X) * Yorg / n;
    vec Ymean = mean(Yorg, 0);

    
    // I-LAMM algorithm
    mat beta(P,m,fill::zeros);
    mat beta_old = beta;
    vec beta0(m,fill::zeros);
    vec beta0_old = beta0;
   
    vec beta_norm(p,fill::zeros);
    vec lambda_Xdim_sqrt = lambda * Xdims_sqrt;
    vec lambda_vec = lambda_Xdim_sqrt;

    double phi = phi0;
    int iter = 0;
    int iter_inner = 0;
    vec iter_inner_vec(max_iter);
    while (iter<max_iter){
        iter = iter+1;
        beta_old = beta;
        beta0_old = beta0;
        lambda_vec = Penalty_diff_ftn(beta_norm, lambda_Xdim_sqrt, penalty, gamma);
        

        // sove locally adaptive LASSO problem
        mat beta_inner = beta;
        vec beta0_inner = beta0;
        vec beta_norm_inner = beta_norm;

        iter_inner = 0;
        while (iter_inner < max_iter) {
            iter_inner = iter_inner + 1;
            phi = max(phi0, phi / c_phi);
            phi = min(phi, 1.0 / phi0);

            mat beta_inner_old = beta_inner;
            mat beta0_inner_old = beta0_inner;
            vec beta_norm_inner_old = beta_norm_inner;

            // compute loss_inner_old and loss_1d_inner_old
            mat res_inner_old = Yorg - repelem(beta0_inner_old, n, 1) - X * beta_inner_old;
            mat loss_mat_inner_old = sQRloss(res_inner_old, tau, h, kernel); // (n,m) mat
            double penalty_term_inner_old = Penalty_ftn(beta_norm_inner_old, lambda_vec, "LASSO", 0);
            double loss_inner_old = mean(mean(loss_mat_inner_old, 0)) + penalty_term_inner_old;

            mat loss_1d_inner_old = - 1.0 *  sQRloss_diff(res_inner_old, tau, h, kernel); // (n,m) mat, loss(Y - Xbeta) -> it needs to multiply -1
            mat X_loss_1d_inner_old = trans(X) * loss_1d_inner_old / n; // (P,m) mat
            

            // update beta and beta0 when the corresponding loss decreases
            int iter_inner_inner = 0;
            while (iter_inner_inner < max_iter) {
                iter_inner_inner = iter_inner_inner + 1;

                // beta0 update
                mat w0 = beta0_inner_old - mean(loss_1d_inner_old, 0) / phi;
                double normw0 = L2_norm_real(w0);
                beta0_inner = LASSO_sol(w0, 1, normw0, 0, 0, 0);


                // betaj update
                for (int j = 0; j < p; j++) {
                    // define indices
                    uvec ind1 = ind_list1[j];
                    uvec ind2 = ind_list2[j];

                    mat wj = beta_inner_old.rows(ind1) - X_loss_1d_inner_old.rows(ind1) / phi;
                    double normwj = L2_norm_real(wj);
                    mat betaj_inner = LASSO_sol(wj, 1, normwj, lambda_vec[j] / phi, 0, 0);
                    beta_inner.rows(ind1) = betaj_inner;

                    beta_norm_inner(j) = L2_norm_real(betaj_inner);
                }

                // compute loss_inner
                mat res_inner = Yorg - repelem(beta0_inner, n, 1) - X * beta_inner;
                mat loss_mat_inner = sQRloss(res_inner, tau, h, kernel); // (n,m) mat
                double penalty_term_inner = Penalty_ftn(beta_norm_inner, lambda_vec, "LASSO", 0);
                double loss_inner = mean(mean(loss_mat_inner, 0)) + penalty_term_inner;
                

                if (loss_inner > loss_inner_old) {
                    phi = c_phi * phi;
                }
                else {
                    break;
                }                
            }

            double residual_inner = L2_norm_real(beta_inner - beta_inner_old);
            if (abs(residual_inner) < threshold) {
                break;
            }
        }

        // beta and beta0 update
        beta = beta_inner;
        beta0 = beta0_inner;
        beta_norm = beta_norm_inner;
        iter_inner_vec(iter) = iter_inner;
        

        double residual = L2_norm_real(beta-beta_old);
        if (abs(residual)<threshold){
            break;
        }
    }
  
    // X: normalized -> restore the actual beta
    for (int i=0; i<n_W_nonzero; i++){
        int j = W_nonzero(i);
        beta.row(j) = beta.row(j) / sqrt(W(j));
    }
  
    List result = List::create(Named("beta") = beta, Named("beta0") = beta0, Named("beta.norm") = beta_norm, 
                               Named("Xdims") = Xdims, Named("iter") = iter, Named("iter.inner") = iter_inner_vec, Named("lambda") = lambda, Named("Xdim.max") = Xdim_max,
                               Named("tau") = tau, Named("h") = h, Named("kernel") = kernel, Named("phi") = phi, Named("penalty") = penalty, Named("gamma") = gamma);
  
    return(result);
}





double get_loss_QM(List X, arma::mat Y, List Xnew_, arma::mat Ynew, double lambda, int Xdim_max, double tau, double h, String kernel, String penalty, double gamma) {

    // model training
    List model = QM_each(X, Y, lambda, Xdim_max, tau, h, kernel, penalty, gamma);
    mat beta = model["beta"];
    vec beta0 = model["beta0"];
    double h_ = model["h"];

    // Xnew data 
    int n2 = Ynew.n_rows;
    List Xnewlist = Make_reduce_dim_matrix(Xnew_, Xdim_max);
    mat Xnew = Xnewlist["X"];

    // compute loss
    mat res_new = Ynew - repelem(beta0, n2, 1) - Xnew * beta;
    mat loss_mat = sQRloss(res_new, tau, h_, kernel); // (n2,m) mat
    double loss = mean(mean(loss_mat, 0));
    loss = log(loss);

    return(loss);
}




double get_loss_CV_QM(List X_, arma::mat Y, double lambda, int Xdim_max, String cv_type, double tau, double h, String kernel, String penalty, double gamma, double cv_const) {

    // model training
    List model = QM_each(X_, Y, lambda, Xdim_max, tau, h, kernel, penalty, gamma);
    mat beta = model["beta"];
    vec beta0 = model["beta0"];
    double h_ = model["h"];

    // Xnew data 
    int n = Y.n_rows;
    List Xlist = Make_reduce_dim_matrix(X_, Xdim_max);
    mat X = Xlist["X"];

    // compute loss
    mat res = Y - repelem(beta0, n, 1) - X * beta;
    mat loss_mat = sQRloss(res, tau, h_, kernel); // (n,m) mat
    double loss = mean(mean(loss_mat, 0));
    loss = log(loss);

    // compute additional penalty term for AIC
    double aic = 0;
    double bic = 0;

    vec beta_norm = model["beta.norm"];
    vec Xdims = model["Xdims"];

    if (cv_type == "AIC" || cv_type == "ABIC") {
        aic = cv_const * sum(Xdims.elem(find(beta_norm != 0))) / n;
    }
    if (cv_type == "BIC" || cv_type == "ABIC") {
        bic = sum(Xdims.elem(find(beta_norm != 0))) * log(n) * cv_const / (2 * n);
    }

    loss = loss + aic + bic;


    return(loss);
}

