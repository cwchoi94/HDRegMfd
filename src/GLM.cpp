// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "PenaltySol.h"
#include "utils.h"
#include "GLM_base.h"
#include "GLM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



// [[Rcpp::export]]
List GLM_each(List Xorg, arma::mat Yorg, double lambda, int Xdim_max, double R, String penalty, String link, double gamma, double phi, double eta, int max_iter, double threshold){

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

    
    // ADMM algorithm
    mat beta(P,m,fill::zeros);
    mat beta_old = beta;
    vec beta0(m,fill::zeros);
    vec beta0_old = beta0;
    double A = R;
    double B = 0;
  
    vec beta_norm(p,fill::zeros);
    vec beta_norm_Xdim_sqrt = Xdims_sqrt % beta_norm; // sqrt(Kj)|beta_j|_2
  
    double residual = 1e5;
    int iter = 0;
    while (iter<max_iter){
        iter = iter+1;
        beta_old = beta;
        beta0_old = beta0;
        

        // beta0 update
        mat theta0 = X * beta + repelem(beta0, n, 1);
        mat psi_1d_theta0 = Psi_1d(theta0, link);
        double phi0 = max(Psi_2d(theta0, link).max(), phi);

        mat w0 = phi0 * beta0 + Ymean - mean(psi_1d_theta0, 0);
        double sigma0 = phi0;
        double normw0 = L2_norm_real(w0);

        beta0 = betaj_update(w0, sigma0, normw0, penalty, 0, 0, 0, beta0, 1);

        // beta update
        for (int j=0; j<p; j++){
            // define indices
            uvec ind1 = ind_list1[j];
            uvec ind2 = ind_list2[j];
            
            // compute components
            mat thetaj = X * beta + repelem(beta0, n, 1);
            mat psi_1d_thetaj = Psi_1d(thetaj, link);
            double phij = max(Psi_2d(thetaj, link).max(),phi);

            mat wj = phij * beta.rows(ind1) + XY.rows(ind1) - trans(X.cols(ind1)) * psi_1d_thetaj / n;
            double sigmaj = phij + eta * Xdims(j);
            double normwj = L2_norm_real(wj);
            double lambdaj = Xdims_sqrt(j) * lambda;
            double nuj = B * Xdims_sqrt(j) + eta*(sum(beta_norm_Xdim_sqrt)-beta_norm_Xdim_sqrt(j)+A-R);

            // betaj update
            mat betaj = betaj_update(wj,sigmaj,normwj,penalty,lambdaj,nuj,gamma,beta.rows(ind1),beta_norm(j));
            beta.rows(ind1) = betaj;
            beta_norm(j) = L2_norm_real(betaj);
            beta_norm_Xdim_sqrt(j) = Xdims_sqrt(j) * beta_norm(j);
        }

        // A update
        if (eta != 0) {
            A = max(R - sum(beta_norm_Xdim_sqrt) - B / eta, 0.0);
        }
        else {
            A = 0.0;
        }
        
    
        // B update
        B = B + eta*(sum(beta_norm_Xdim_sqrt) + A - R);
    
        residual = L2_norm_real(X * (beta-beta_old)) / sqrt(n);
        if (abs(residual)<threshold){
            break;
        }
    }
  
    // X: normalized -> restore the actual beta
    for (int i=0; i<n_W_nonzero; i++){
        int j = W_nonzero(i);
        beta.row(j) = beta.row(j) / sqrt(W(j));
    }
  
    List result = List::create(Named("beta") = beta, Named("beta0") = beta0, Named("beta.norm") = beta_norm, Named("A") = A, Named("B") = B,
                               Named("Xdims") = Xdims, Named("iter") = iter, Named("lambda") = lambda, Named("Xdim.max") = Xdim_max, Named("R") = R, 
                               Named("phi") = phi, Named("penalty") = penalty, Named("gamma") = gamma);
  
    return(result);
}



double get_loss_GLM(List X, arma::mat Y, List Xnew_, arma::mat Ynew, double lambda, int Xdim_max, double R, String penalty, String link, double gamma, double phi){
         
    // model training
    List model = GLM_each(X, Y, lambda, Xdim_max, R, penalty, link, gamma, phi);
    mat beta = model["beta"];
    vec beta0 = model["beta0"];
  
    // Xnew data 
    int n2 = Ynew.n_rows;
    List Xnewlist = Make_reduce_dim_matrix(Xnew_,Xdim_max);
    mat Xnew = Xnewlist["X"];
  
    // compute theta
    mat theta = Xnew * beta + repelem(beta0, n2, 1);
    mat psi_theta = Psi(theta, link);

    // compute loss
    double loss = accu(psi_theta - Ynew % theta) / n2;

    return(loss);
}



double get_loss_CV_GLM(List X_, arma::mat Y, double lambda, int Xdim_max, double R, String cv_type, String penalty, String link, double gamma, double phi, double cv_const) {

    // model training
    List model = GLM_each(X_, Y, lambda, Xdim_max, R, penalty, link, gamma, phi);
    mat beta = model["beta"];
    vec beta0 = model["beta0"];

    // Xnew data 
    int n = Y.n_rows;
    List Xlist = Make_reduce_dim_matrix(X_, Xdim_max);
    mat X = Xlist["X"];

    // compute theta
    mat theta = X * beta + repelem(beta0, n, 1);
    mat psi_theta = Psi(theta, link);

    // compute loss
    double loss = accu(psi_theta - Y % theta) / n;

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
