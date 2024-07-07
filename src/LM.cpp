// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "PenaltySol.h"
#include "utils.h"
#include "LM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
List LM_each(List Xorg, arma::mat LogY, arma::vec Ymu, Function inner, double lambda, int Xdim_max, double R, String penalty, double phi, double gamma, double eta, int max_iter, double threshold){

    int p = Xorg.size();
    int n = LogY.n_rows;
    int m = LogY.n_cols;
    
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
    mat XY = trans(X) * LogY / n;
        
    
    // ADMM algorithm
    mat beta(P,m,fill::zeros);
    mat beta_old = beta;
    double A = R;
    double B = 0;
  
    vec beta_norm(p,fill::zeros);
    vec beta_norm_Xdim_sqrt = Xdims_sqrt % beta_norm; // sqrt(Kj)|beta_j|_2
  
    double residual = 1e5;
    int iter = 0;
    while (iter<max_iter){
        iter = iter+1;
        beta_old = beta;
    
        // beta update
        for (int j=0; j<p; j++){
            // define indices
            uvec ind1 = ind_list1[j];
            uvec ind2 = ind_list2[j];
            
            // compute components
            mat wj = phi * beta.rows(ind1) + XY.rows(ind1) - XX.rows(ind1) * beta;
            double sigmaj = phi + eta * Xdims(j);
            double normwj = L2_norm(wj, Ymu, inner);
            double lambdaj = Xdims_sqrt(j) * lambda;
            double nuj = B * Xdims_sqrt(j) + eta*(sum(beta_norm_Xdim_sqrt)-beta_norm_Xdim_sqrt(j)+A-R);

            // betaj update
            mat betaj = betaj_update(wj,sigmaj,normwj,penalty,lambdaj,nuj,gamma,beta.rows(ind1),beta_norm(j));
            beta.rows(ind1) = betaj;
            beta_norm(j) = L2_norm(betaj,Ymu,inner);
            beta_norm_Xdim_sqrt(j) = Xdims_sqrt(j) * beta_norm(j);
        }
    
        // A update
        A = max(R - sum(beta_norm_Xdim_sqrt) - B / eta, 0.0);
    
        // B update
        B = B + eta*(sum(beta_norm_Xdim_sqrt) + A - R);
    
        residual = L2_norm(X * (beta-beta_old),Ymu,inner) / sqrt(n);
        if (abs(residual)<threshold){
            break;
        }
    }
  
    // X: normalized -> restore the actual beta
    for (int i=0; i<n_W_nonzero; i++){
        int j = W_nonzero(i);
        beta.row(j) = beta.row(j) / sqrt(W(j));
    }
  
    List result = List::create(Named("beta") = beta, Named("A") = A, Named("B") = B, Named("Xdims") = Xdims, Named("iter") = iter,                                
                               Named("lambda") = lambda, Named("Xdim.max") = Xdim_max, Named("R") = R, Named("phi") = phi,
                               Named("penalty") = penalty, Named("gamma") = gamma);
  
    return(result);
}


double get_loss_LM(List X, arma::mat LogY, List Xnew_, arma::mat LogYnew, arma::vec Ymu, Function inner, double lambda, int Xdim_max, double R, String penalty, double phi, double gamma){
         
    // model training
    List model = LM_each(X,LogY,Ymu,inner,lambda,Xdim_max,R,penalty,phi,gamma);
    mat beta = model["beta"];
  
    // Xnew data 
    int n2 = LogYnew.n_rows;
    List Xnewlist = Make_reduce_dim_matrix(Xnew_,Xdim_max);
    mat Xnew = Xnewlist["X"];
  
    // compute Yhat
    mat LogYhat = Xnew * beta;

    // compute rmse
    double rmse = L2_norm(LogYhat-LogYnew,Ymu,inner)/sqrt(n2);

    return(rmse);
}




