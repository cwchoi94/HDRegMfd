// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "PenaltySol.h"
#include "utils.h"
#include "Linear.h"

using namespace Rcpp;
using namespace std;
using namespace arma;


List Make_reduce_dim_matrix(List Xorg, int Xdim_max){
    // Join each Xj in Xorg
    // Also make vector of dimensions
    int p = Xorg.size();
    vec Xdims(p);
    mat Xj = Xorg[0];
    int pj = Xj.n_cols;
    if (pj>Xdim_max){
        Xj = Xj.cols(0,Xdim_max-1);
    }
    mat X = Xj;
    Xdims(0) = Xj.n_cols;
    for (int j=1; j<p; j++){
        mat Xj = Xorg[j];
        int pj = Xj.n_cols;
        if (pj>Xdim_max){
            Xj = Xj.cols(0,Xdim_max-1);
        }
        Xdims(j) = Xj.n_cols;
        X = join_rows(X,Xj);
    }

    List result = List::create(Named("X") = X, Named("Xdims") = Xdims);

    return(result);
}


List Make_dimension_indices(arma::vec Xdims_cumul){
    int p = Xdims_cumul.size()-1;
    int P = Xdims_cumul(p);

    List ind_list1(p);
    List ind_list2(p);
    for (int j=0; j<p; j++){
        // generate indices
        // ind1 = a:b, ind2 = -(a:b)
        int a = Xdims_cumul[j];
        int b = Xdims_cumul[j+1]-1;
        std::vector<int> ind1_vec;
        std::vector<int> ind2_vec;
        for (int i=0; i<P; i++){
            if ((a<=i) & (i<=b)){
                ind1_vec.push_back(i);
            }
            else{
                ind2_vec.push_back(i);
            }
        }
        uvec ind1 = arma::conv_to<uvec>::from(ind1_vec);
        uvec ind2 = arma::conv_to<uvec>::from(ind2_vec);
        
        ind_list1[j] = ind1;
        ind_list2[j] = ind2;
    }

    List ind_list = List::create(Named("ind1") = ind_list1, Named("ind2") = ind_list2);
    
    return (ind_list);
}


arma::mat Linear_betaj_update(arma::mat wj, double sigmaj, double normwj, String penalty, double lambdaj, double kappaj, double gamma, mat betaj, double betaj_norm){
    mat betaj_new;
    bool is_sol2 = (kappaj+lambdaj<0) & (normwj ==0);

    if (penalty=="LASSO"){
        if (is_sol2){
            betaj_new = LASSO_sol2(betaj,sigmaj,betaj_norm,lambdaj,kappaj,gamma);
        }
        else{
            betaj_new = LASSO_sol(wj,sigmaj,normwj,lambdaj,kappaj,gamma);
        }
    }
    else if (penalty=="SCAD"){
        if (is_sol2){
            betaj_new = SCAD_sol2(betaj,sigmaj,betaj_norm,lambdaj,kappaj,gamma);
        }
        else{
            betaj_new = SCAD_sol(wj,sigmaj,normwj,lambdaj,kappaj,gamma);
        }
    }
    else if (penalty=="MCP"){
        if (is_sol2){
            betaj_new = MCP_sol2(betaj,sigmaj,betaj_norm,lambdaj,kappaj,gamma);
        }
        else{
            betaj_new = MCP_sol(wj,sigmaj,normwj,lambdaj,kappaj,gamma);
        }
    }

    return(betaj_new);
}


// [[Rcpp::export]]
List LM_each(List Xorg, arma::mat LogY, arma::vec Ymu, Function inner, double lambda, int Xdim_max, double R, double phi, String penalty, double gamma, double eta, int max_iter, double threshold){

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
            double kappaj = B * Xdims_sqrt(j) + eta*(sum(beta_norm_Xdim_sqrt)-beta_norm_Xdim_sqrt(j)+A-R);

            // betaj update
            mat betaj = Linear_betaj_update(wj,sigmaj,normwj,penalty,lambdaj,kappaj,gamma,beta.rows(ind1),beta_norm(j));
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


double get_loss_LM(List X, arma::mat LogY, List Xnew_, arma::mat LogYnew, arma::vec Ymu, Function inner, double lambda, int Xdim_max, double R, double phi, String penalty, double gamma){
         
    // model training
    List model = LM_each(X,LogY,Ymu,inner,lambda,Xdim_max,R,phi,penalty,gamma);
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




