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




arma::cube sqrt_mat_cube(arma::cube x, double tol = 1e-8) {
    // x: (p,r,r) matrix
    // 
    // return matrix-wise inverse, i.e., y[i] = inv(x[i]) for 1<=i<=p
    // y: (p,r,r) matrix

    int p = x.n_rows;
    int r = x.n_cols;
    mat tol_mat = tol * arma::eye(r, r);

    cube y(p, r, r);
    for (int i = 0; i < p; i++) {
        mat xi = x.row(i);
        y.row(i) = arma::sqrtmat_sympd(xi + tol_mat);
    }

    return y;
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
     List kde_1d_sq_list(p);
     List kde_1d_inv_list(p);
     for (int j = 0; j < p; j++) {
         int ind1 = g * j;
         int ind2 = g * (j + 1) - 1;

         cube kde_1d_j = kde_1d.rows(ind1, ind2);
         cube kde_1d_sq_j = sqrt_mat_cube(kde_1d_j);
         cube kde_1d_inv_j = kde_1d_inv.rows(ind1, ind2);

         kde_1d_list[j] = kde_1d_j;
         kde_1d_sq_list[j] = kde_1d_sq_j;
         kde_1d_inv_list[j] = kde_1d_inv_j;
     }

     // change the dimension of proj 
     // (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2, r1, r2)  ->  (p1*p2, g1*r1, g2*r2) cube
     // also multiply weights so that we can compute the integral quickly
     cube proj_new(p * p, g * r, g * r);
     for (int j1 = 0; j1 < p; j1++) {
         for (int j2 = 0; j2 < p; j2++) {

             mat proj_jj2_new(g * r, g * r);
             for (int k1 = 0; k1 < g; k1++) {
                 for (int k2 = 0; k2 < g; k2++) {
                     int single_ind = multi_4d_ind_to_single(j1, j2, k1, k2, p, g);
                     mat projg_jj2_kk2 = proj.row(single_ind); // (r1,r2) mat
                     proj_jj2_new.submat(r * k1, r * k2, r * (k1 + 1) - 1, r * (k2 + 1) - 1) = weights[k2] * projg_jj2_kk2;
                 }
             }

             int ind_new = p * j1 + j2;
             proj_new.row(ind_new) = proj_jj2_new;
         }
     }


     // compute a tildem list
     // p list - (g,r,m) -> (p,g*r,m) cube
     cube tildem(p, g * r, m);
     for (int j = 0; j < p; j++) {
         cube kvalues_all_j = kvalues_all[j]; // (g,n,r) cube
         cube kde_1d_inv_j = kde_1d_inv_list[j];   // (g,r,r) cube

         mat tildem_j(g * r, m);
         for (int k = 0; k < g; k++) {
             mat kde_1d_inv_jk = kde_1d_inv_j.row(k);
             mat kvalues_all_jk = kvalues_all_j.row(k);

             if (degree == 0) {
                 tildem_j.rows(r * k, r * (k + 1) - 1) = kde_1d_inv_jk * kvalues_all_jk * LogY / n;
             }
             else if (degree > 0) {
                 tildem_j.rows(r * k, r * (k + 1) - 1) = kde_1d_inv_jk * kvalues_all_jk.t() * LogY / n;
             }

         }

         tildem.row(j) = tildem_j;
     }

     return List::create(Named("kde.1d") = kde_1d_list, Named("kde.1d.sq") = kde_1d_sq_list, 
         Named("tildem") = tildem, Named("proj") = proj_new, Named("bandwidths") = bandwidths,
         Named("grids") = grids, Named("weights") = weights, Named("r") = r);
 }



 // [[Rcpp::export]]
 List SBF_preprocessing_reduce_dim(List SBF_comp, double Xdim_max, arma::mat index_mat) {

     uvec col_indices_uvec = arma::find(index_mat.col(2) <= Xdim_max);

     List kde_1d = SBF_comp["kde.1d"];
     List kde_1d_sq = SBF_comp["kde.1d.sq"];
     cube tildem = SBF_comp["tildem"];
     cube proj = SBF_comp["proj"];
     vec bandwidths = SBF_comp["bandwidths"];
     vec grids = SBF_comp["grids"];
     vec weights = SBF_comp["weights"];
     int r = SBF_comp["r"];

     int P = kde_1d.size();
     int p = col_indices_uvec.size();
     int g = weights.size();
     int m = tildem.n_slices;


     vec bandwidths_reduced = bandwidths.elem(col_indices_uvec);
     List kde_1d_reduced(p);
     List kde_1d_sq_reduced(p);
     cube tildem_reduced(p, g * r, m);
     cube proj_reduced(p * p, g * r, g * r);
     for (int i = 0; i < p; i++) {
         int j = col_indices_uvec(i);

         // kde_1d_j
         cube kde_1d_j = kde_1d[j];
         cube kde_1d_sq_j = kde_1d_sq[j];

         kde_1d_reduced[i] = kde_1d_j;
         kde_1d_sq_reduced[i] = kde_1d_j;

         // tildem_j
         mat tildem_j = tildem.row(j);
         tildem_reduced.row(i) = tildem_j;

         // proj_jj2
         for (int i2 = 0; i2 < p; i2++) {
             int j2 = col_indices_uvec(i2);
             mat proj_jj2 = proj.row(P * j + j2);
             proj_reduced.row(p * i + i2) = proj_jj2;
         }
     }

     return List::create(Named("kde.1d") = kde_1d_reduced, Named("kde.1d.sq") = kde_1d_sq_reduced, 
         Named("tildem") = tildem_reduced, Named("proj") = proj_reduced, Named("bandwidths") = bandwidths_reduced,
         Named("grids") = grids, Named("weights") = weights, Named("r") = r);
 }