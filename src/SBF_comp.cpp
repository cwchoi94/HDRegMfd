// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include "utils.h"
#include "SBF_utils.h"
#include "KDE.h"
#include "SBF_comp.h"

using namespace Rcpp;
using namespace std;
using namespace arma;




arma::mat Reduced_X_mat(arma::mat X, arma::mat index_mat, int Xdim_max) {
    // find indices of (j,k) with k <= Xdim_max
    uvec col_indices_uvec = arma::find(index_mat.col(2) <= Xdim_max);

    // compute the reduced X with the above row_indices
    mat X_reduced = X.cols(col_indices_uvec);

    return X_reduced;
}


List Reduced_X_list(List X_list, arma::mat index_mat, int Xdim_max) {

    int kfold = X_list.size();

    List X_list_reduced(kfold);
    for (int i = 0; i < kfold; i++) {
        mat X = X_list[i];
        mat X_reduced = Reduced_X_mat(X, index_mat, Xdim_max);
        X_list_reduced[i] = X_reduced;
    }

    return X_list_reduced;
}




 // [[Rcpp::export]]
 List SBF_preprocessing(arma::mat X, arma::mat LogY, arma::vec bandwidths, arma::vec grids, arma::vec weights, int degree, String Kdenom_method) {

     List kde = KDE_(X, bandwidths, grids, weights, degree, Kdenom_method, true);

     int n = X.n_rows;
     int p = X.n_cols;
     int g = kde["ngrid"];
     int r = degree + 1;
     int m = LogY.n_cols;

     cube proj = kde["proj"];
     cube kde_1d = kde["kde.1d"];
     cube kde_1d_inv = kde["kde.1d.inv"];
     List kvalues_all = kde["kvalues"];


     // make kde_1d and kde_1d_inv lists
     List kde_1d_list(p);
     List kde_1d_inv_list(p);
     for (int j = 0; j < p; j++) {
         int ind1 = g * j;
         int ind2 = g * (j + 1) - 1;
         cube kde_1d_j = kde_1d.rows(ind1, ind2);
         cube kde_1d_inv_j = kde_1d_inv.rows(ind1, ind2);
         kde_1d_list[j] = kde_1d_j;
         kde_1d_inv_list[j] = kde_1d_inv_j;
     }

     // compute a tildem list
     List tildem(p);
     for (int j = 0; j < p; j++) {
         cube kvalues_all_j = kvalues_all[j]; // (g,,n,r) cube
         cube kde_1d_inv_j = kde_1d_inv_list[j];   // (g,r,r) cube

         cube tildem_j(g, r, m);
         for (int k = 0; k < g; k++) {
             mat kde_1d_inv_jk = kde_1d_inv_j.row(k);
             mat kvalues_all_jk = kvalues_all_j.row(k);

             if (degree == 0) {
                 tildem_j.row(k) = kde_1d_inv_jk * kvalues_all_jk * LogY / n;
             }
             else if (degree > 0) {
                 tildem_j.row(k) = kde_1d_inv_jk * kvalues_all_jk.t() * LogY / n;
             }

         }

         tildem[j] = tildem_j;
     }

     return List::create(Named("tildem") = tildem, Named("kde.1d") = kde_1d_list, Named("proj") = proj, Named("bandwidths") = bandwidths,
         Named("grids") = grids, Named("weights") = weights);
 }



 // [[Rcpp::export]]
 List SBF_preprocessing_reduce_dim(List SBF_comp, double Xdim_max, arma::mat index_mat) {

     uvec col_indices_uvec = arma::find(index_mat.col(2) <= Xdim_max);

     List tildem = SBF_comp["tildem"];
     List kde_1d = SBF_comp["kde.1d"];
     cube proj = SBF_comp["proj"];
     vec bandwidths = SBF_comp["bandwidths"];
     vec grids = SBF_comp["grids"];
     vec weights = SBF_comp["weights"];

     int P = kde_1d.size();
     int p = col_indices_uvec.size();
     int g = weights.size();
     int r = proj.n_cols;

     vec bandwidths_reduced = bandwidths.elem(col_indices_uvec);
     List tildem_reduced(p);
     List kde_1d_reduced(p);
     cube proj_reduced(p * p * g * g, r, r);
     for (int i = 0; i < p; i++) {
         int j = col_indices_uvec(i);
         cube tildem_j = tildem[j];
         cube kde_1d_j = kde_1d[j];

         tildem_reduced[i] = tildem_j;
         kde_1d_reduced[i] = kde_1d_j;

         for (int i2 = 0; i2 < p; i2++) {
             int j2 = col_indices_uvec(i2);
             IntegerVector ind_range = multi_2d_ind_to_single_range(j, j2, P, g);
             IntegerVector ind_range_reduced = multi_2d_ind_to_single_range(i, i2, p, g);

             cube proj_jj2 = proj.rows(ind_range(0), ind_range(1));
             proj_reduced.rows(ind_range_reduced(0), ind_range_reduced(1)) = proj_jj2;
         }
     }

     return List::create(Named("tildem") = tildem_reduced, Named("kde.1d") = kde_1d_reduced, Named("proj") = proj_reduced, Named("bandwidths") = bandwidths_reduced,
         Named("grids") = grids, Named("weights") = weights);
 }