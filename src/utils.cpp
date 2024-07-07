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


double L2_norm_real(arma::mat Y){
    double z = sqrt(accu(square(Y)));
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




List Make_reduce_dim_matrix(List Xorg, int Xdim_max) {
    // Join each Xj in Xorg
    // Also make vector of dimensions
    int p = Xorg.size();
    vec Xdims(p);
    mat Xj = Xorg[0];
    int pj = Xj.n_cols;
    if (pj > Xdim_max) {
        Xj = Xj.cols(0, Xdim_max - 1);
    }
    mat X = Xj;
    Xdims(0) = Xj.n_cols;
    for (int j = 1; j < p; j++) {
        mat Xj = Xorg[j];
        int pj = Xj.n_cols;
        if (pj > Xdim_max) {
            Xj = Xj.cols(0, Xdim_max - 1);
        }
        Xdims(j) = Xj.n_cols;
        X = join_rows(X, Xj);
    }

    List result = List::create(Named("X") = X, Named("Xdims") = Xdims);

    return(result);
}


List Make_dimension_indices(arma::vec Xdims_cumul) {
    int p = Xdims_cumul.size() - 1;
    int P = Xdims_cumul(p);

    List ind_list1(p);
    List ind_list2(p);
    for (int j = 0; j < p; j++) {
        // generate indices
        // ind1 = a:b, ind2 = -(a:b)
        int a = Xdims_cumul[j];
        int b = Xdims_cumul[j + 1] - 1;
        std::vector<int> ind1_vec;
        std::vector<int> ind2_vec;
        for (int i = 0; i < P; i++) {
            if ((a <= i) & (i <= b)) {
                ind1_vec.push_back(i);
            }
            else {
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

