// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include "utils.h"
#include "SBF_utils.h"
#include "KDE.h"

using namespace Rcpp;
using namespace std;
using namespace arma;





arma::rowvec power_vec(double u, int r) {
    rowvec z(r,fill::ones);
    for (int i = 1; i < r; i++) {
        z(i) = pow(u, i);
    }

    return z;
}




arma::cube inv_cube(arma::cube x, double tol = 1e-8) {
    // x: (p,r,r) matrix
    // 
    // return matrix-wise inverse, i.e., y[i] = inv(x[i]) for 1<=i<=p
    // y: (p,r,r) matrix

    int p = x.n_rows;
    int r = x.n_cols;
    mat tol_mat = tol * arma::eye(r, r);

    cube y(p,r,r);
    for (int i = 0; i < p; i++) {
        mat xi = x.row(i);
        y.row(i) = inv_sympd(xi + tol_mat);
    }

    return y;
}


arma::cube transpose_cube(arma::cube x) {
    // x: (p*q,r,r') matrix (p=q and r=r')
    // 
    // return:
    // t(x): (q*p,r',r) matrix

    int p = sqrt(x.n_rows);
    int r = x.n_cols;

    cube y(p * p, r, r);
    for (int j = 0; j < p; j++) {
        for (int k = 0; k < p; k++) {
            mat x_each = x.row(p * j + k);
            y.row(p * k + j) = x_each.t();
        }
    }

    return y;
}




arma::mat compute_rot_mat(arma::vec x, arma::vec y) {

    double c = 0;
    int p = x.n_elem;
    mat R(p, p);

    if (norm(x, 2) != 0) {
        c = norm(y, 2) / norm(x, 2);

        if (p > 1) {
            // Normalize x and y
            x /= norm(x, 2);
            y /= norm(y, 2);

            // Compute cos(theta) and sin(theta)
            double cos_theta = dot(x, y);
            double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

            // Calculate u as a unit vector orthogonal to x in the direction of y
            vec u = y - cos_theta * x;
            u /= norm(u, 2);

            // Define the rotation matrix
            R = arma::eye(p, p) - x * x.t() - u * u.t();
            mat rotation_sub = { {cos_theta, -sin_theta}, {sin_theta, cos_theta} };
            R += join_horiz(x, u) * rotation_sub * join_horiz(x, u).t();
        }
        else {
            R = arma::eye(p, p);
        }
    }

    return c * R;
}




// kernel functions

double K0(double x) {
    // epanechnikov kernel with support [-1,1]
    if (std::abs(x) < 1) {
        return 0.75 * (1 - x * x);
    }
    else {
        return 0.0;
    }
}


double intK0(double x) {
    // integrate the epanechnikov kernel from -infty to x
    if (x < -1) {
        return 0.0;
    }
    else if (std::abs(x) < 1) {
        return 0.25 * (3 * x - x * x * x + 2);
    }
    else if (x > 1) {
        return 1.0;
    }
    return 0.0;
}


double Kh(double v, double h) {
    return K0(v / h) / h;
}


double Kh_denom_exact(double v, double h) {
    return intK0((1 - v) / h) - intK0((0 - v) / h);
}




////// KDE code


 // [[Rcpp::export]]
 List normalized_Kernel(arma::mat X, arma::vec bandwidths, arma::vec grids, arma::vec weights, int degree, String Kdenom_method) {
    // compute the normalized kernel K_h(x,X)
    // 
    // Input:
    // - X: (n,p) matrix
    // - grids, weights: g vector
    // - bandwidths: p vector
    //  
    // Return:
    // - kvalues_all: p list - (g,n,degree+1) cube with (k,i,_) entries of power_vec((x_k-X_ij)/h_j) * K_h(x_k,X_ij)

    int n = X.n_rows;
    int p = X.n_cols;
    int g = grids.size();
    int r = degree + 1;

    List vec_components_all(p);
    List kvalues_all(p);
    for (int j = 0; j < p; j++) {
        vec Xj = X.col(j);
        double h = bandwidths[j];

        cube vec_components(g, n, r);
        cube kvalues(g, n, r);
        double Kdenom = 0;
        for (int i = 0; i < n; i++) {
            mat vec_components_each(g, r);
            mat kvalues_each(g, r);
            for (int k = 0; k < g; k++) {
                double u = (Xj[i] - grids[k]) / h;
                rowvec u_vec = power_vec(u, r);

                vec_components_each.row(k) = u_vec;
                kvalues_each.row(k) = K0(u) * u_vec / h;
            }

            if (Kdenom_method == "numeric") {
                Kdenom = numerical_integral_1d(weights, kvalues_each.col(0));
            }
            else if (Kdenom_method == "exact") {
                Kdenom = Kh_denom_exact(Xj[i], h);
            }

            if (Kdenom < 1e-8) {
                kvalues_each *= 0; // set to zero
            }
            else {
                kvalues_each /= Kdenom;
            }
            kvalues.col(i) = kvalues_each;
            vec_components.col(i) = vec_components_each;
        }

        vec_components_all[j] = vec_components;
        kvalues_all[j] = kvalues;
    }

    return List::create(Named("kvalues") = kvalues_all, Named("vec_components") = vec_components_all);
}




 // [[Rcpp::export]]
 List KDE_(arma::mat X, arma::vec bandwidths, arma::vec grids, arma::vec weights, int degree, String Kdenom_method, bool is_proj) {
     
     int n = X.n_rows;
     int p = X.n_cols;
     int g = grids.size();
     int r = degree + 1;

     vec unit_vec(r);
     unit_vec(0) = 1;

     // compute normalized kernel values and vectorized components: p list - (g,n,r) cube
     List object = normalized_Kernel(X, bandwidths, grids, weights, degree, Kdenom_method);
     List kvalues_all = object["kvalues"];
     List vec_components_all = object["vec_components"];

     // kde_1d: (p,g,r,r) = (p*g,r,r) cube
     cube kde_1d(p * g, r, r);
     for (int j = 0; j < p; j++) {
         cube kvalues_all_j = kvalues_all[j];
         cube vec_components_all_j = vec_components_all[j];

         // kde_1d_j: (g,r,r) cube
         cube kde_1d_j(g, r, r);
         vec kde_1d_j_00(g);
         for (int k = 0; k < g; k++) {
             mat kvalues_all_jk = kvalues_all_j.row(k);
             mat vec_components_all_jk = vec_components_all_j.row(k);

             int ind = g * j + k;

             // don't have to divide by n
             if (degree == 0) {
                 kde_1d.row(ind) = kvalues_all_jk * vec_components_all_jk.t();
             }
             else if (degree > 0) {
                 kde_1d.row(ind) = kvalues_all_jk.t() * vec_components_all_jk;
             }

             kde_1d_j_00(k) = kde_1d(ind, 0, 0);
         }

         double denom = numerical_integral_1d(weights, kde_1d_j_00);
         if (denom != 0) {
             for (int k = 0; k < g; k++) {
                 int ind = g * j + k;
                 kde_1d.row(ind) = kde_1d.row(ind) / denom;
             }
         }
     }
     
     cube kde_1d_inv = inv_cube(kde_1d);
     


     // kde_2d, proj: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2, r1, r2) cubes, with the identification of indices using functions in utils_SBF 
     // In fact, p1=p2, g1=g2, r1=r2. We only denote to distingush them.

     cube kde_2d(p * p * g * g, r, r, fill::zeros);
     cube proj(p * p * g * g, r, r, fill::zeros);
     for (int j1 = 0; j1 < p; j1++) {
         cube kvalues_all_j1 = kvalues_all[j1];

         for (int j2 = 0; j2 < p; j2++) {
             cube kvalues_all_j2 = kvalues_all[j2];

             if (j2 < j1) {
                 IntegerVector ind_range = multi_2d_ind_to_single_range(j1, j2, p, g);
                 IntegerVector ind_range_inv = multi_2d_ind_to_single_range(j2, j1, p, g);
                 cube kde_2d_tmp = kde_2d.rows(ind_range_inv(0), ind_range_inv(1));
                 kde_2d.rows(ind_range(0), ind_range(1)) = transpose_cube(kde_2d_tmp);
             }
             if (j1 == j2) {
                 continue;
             }
             else if (j1 < j2) {
                 for (int k1 = 0; k1 < g; k1++) {
                     for (int k2 = 0; k2 < g; k2++) {
                         mat kvalues_all_jk1 = kvalues_all_j1.row(k1); // (n,r1) matrix
                         mat kvalues_all_jk2 = kvalues_all_j2.row(k2); // (n,r2) matrix
                         
                         // don't have to divide by n
                         int single_ind = multi_4d_ind_to_single(j1, j2, k1, k2, p, g);
                         if (degree == 0) {
                             kde_2d.row(single_ind) = kvalues_all_jk1 * kvalues_all_jk2.t() / n;
                         }
                         else if (degree > 0) {
                             kde_2d.row(single_ind) = kvalues_all_jk1.t() * kvalues_all_jk2 / n;
                         }
                     }
                 }
             }


             // adjust constants so that  numerical_integral_3d_to_2d to satisfy: integral of kde_1d_j_inv * kde_2djj2 * u_0 = u_0
             for (int k1 = 0; k1 < g; k1++) {
                 mat kde_1d_jk = kde_1d.row(j1 * g + k1); // (r1,r1) matrix
                 mat kde_1d_jk_inv = kde_1d_inv.row(j1 * g + k1); // (r1,r1) matrix
                 IntegerVector ind_range = multi_3d_ind_to_single_range(j1, j2, k1, p, g);
                 cube kde_2d_jj2_k = kde_2d.rows(ind_range(0), ind_range(1)); // (g2,r1,r2) cube

                 // compute a rotation matrx for kde_2d
                 mat denom_mat_2d = numerical_integral_3d_to_2d(weights, kde_2d_jj2_k); // (r1,r2) matrix
                 denom_mat_2d = kde_1d_jk_inv * denom_mat_2d; // (r1,r2) matrix
                 mat rot_mat_2d = compute_rot_mat(denom_mat_2d.col(0), unit_vec);

                 for (int k2 = 0; k2 < g; k2++) {
                     int single_ind = multi_4d_ind_to_single(j1, j2, k1, k2, p, g);
                     mat kde_2d_jj2_kk2 = kde_2d.row(single_ind);
                     kde_2d_jj2_kk2 = rot_mat_2d * kde_2d_jj2_kk2;

                     kde_2d.row(single_ind) = kde_2d_jj2_kk2;
                     proj.row(single_ind) = kde_1d_jk_inv * kde_2d_jj2_kk2;
                 }

                 if (is_proj) {
                     cube proj_jj2_k = proj.rows(ind_range(0), ind_range(1));

                     // compute a rotation matrx for kde_2d
                     mat denom_mat_proj = numerical_integral_3d_to_2d(weights, proj_jj2_k); // (r1,r2) matrix
                     mat rot_mat_proj = compute_rot_mat(denom_mat_proj.col(0), unit_vec);

                     for (int k2 = 0; k2 < g; k2++) {
                         int single_ind = multi_4d_ind_to_single(j1, j2, k1, k2, p, g);
                         mat proj_jj2_kk2 = proj.row(single_ind);

                         proj.row(single_ind) = rot_mat_proj * proj_jj2_kk2;
                     }

                 }

             }
         }
     }

     if (is_proj) {
         return List::create(Named("kvalues") = kvalues_all, Named("kde.1d") = kde_1d, Named("kde.1d.inv") = kde_1d_inv, Named("proj") = proj,
             Named("bandwidths") = bandwidths, Named("grids") = grids, Named("weights") = weights, Named("ngrid") = g);
     }
     else {
         return List::create(Named("kvalues") = kvalues_all, Named("kde.1d") = kde_1d, Named("kde.1d.inv") = kde_1d_inv, Named("kde.2d") = kde_2d,
             Named("bandwidths") = bandwidths, Named("grids") = grids, Named("weights") = weights, Named("ngrid") = g);
     }
 }



