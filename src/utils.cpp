// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include "utils.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



arma::vec SEXP_to_vec(SEXP X){
    NumericVector X_ = X;
    vec X__ = as<arma::vec>(wrap(X_));
    return (X__);
}


arma::mat SEXP_to_mat(SEXP X){
    NumericMatrix X_ = X;
    mat X__ = arma::mat(X_.begin(),X_.nrow(),X_.ncol(),false);
    return (X__);
}


double L2_norm(arma::mat Y, arma::vec Ymu, Function inner){    

    int n = Y.n_rows;
    double z = 0;

    for (int i = 0; i < n; i++) {
        SEXP x = inner(Y.row(i),Y.row(i),Ymu);
        z = z + pow(*REAL(x), 1);
    }
    z = sqrt(z);
    
    return (z);
}


int get_min_idx(arma::vec x, double threshold, int type){
    int r = x.size();
    int opt_idx = 0;
    double xmin = min(x);

    if (type==0){
        // descending
        for (int i=r-1; i>=0; i--){
            if (x(i)<xmin+threshold){
                opt_idx = i;
                break;
            }
        }
    }
    else{
        // ascending
        for (int i=0; i<r; i++){
            if (x(i)<xmin+threshold){
                opt_idx = i;
                break;
            }
        }
    }

    return (opt_idx);
}



