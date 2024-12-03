// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "PenaltySol.h"
#include "Geometry.h"
#include "utils.h"
#include "SBF_utils.h"
#include "SBF_comp.h"
#include "AM.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



arma::cube clone_cube(const arma::cube& original) {
    arma::cube clone = original;
    return clone;
}


double L2_mat_inner_SBF(arma::mat Yj1, arma::mat Yj2, arma::cube kde_1d_sq_j, arma::vec weights, arma::vec Ymu, String Yspace) {
    // Yj1, Yj2: (g*r,m) mats
    // kde_1d_sq_j: (g,r,r) mat

    int g = kde_1d_sq_j.n_rows;
    int r = kde_1d_sq_j.n_slices;
    int m = Yj1.n_cols;

    double z = 0;
    for (int k = 0; k < g; k++) {
        mat Y_jk1 = Yj1.rows(r * k, r * (k + 1) - 1);
        mat Y_jk2 = Yj2.rows(r * k, r * (k + 1) - 1);

        mat kde_1d_sq_jk = kde_1d_sq_j.row(k); // (r,r) matrix      
        mat Y_jk1_ = kde_1d_sq_jk * Y_jk1; // (r,m) matrix
        mat Y_jk2_ = kde_1d_sq_jk * Y_jk2; // (r,m) matrix

        vec z_k = inner(Y_jk1_, Y_jk2_, Ymu, Yspace);

        z += weights[k] * sum(z_k);
    }

    return z;
}


double L2_mat_norm_SBF(arma::mat Yj, arma::cube kde_1d_sq_j, arma::vec weights, arma::vec Ymu, String Yspace) {
    double z = L2_mat_inner_SBF(Yj, Yj, kde_1d_sq_j, weights, Ymu, Yspace);
    return sqrt(z);
}


arma::mat hatmj_update(arma::mat tmp_hatmj, double sigma, double norm_hatmj, String penalty, double lambda, double nu, double gamma, arma::mat mhat_j, double norm_mhat_j) {
    mat mhat_j_new = betaj_update(tmp_hatmj, sigma, norm_hatmj, penalty, lambda, nu, gamma, mhat_j, norm_mhat_j);
    return(mhat_j_new);
}



// [[Rcpp::export]]
List AM_each(List SBF_comp, arma::vec Ymu, String Yspace, double lambda, double R, String penalty, double gamma, double phi, double eta, int max_iter, double threshold) {

    // p = the number of (j,k) with k <= Xdim_max
    // tildem: (p, g*r, m) cube
    // kde_1d_sq: p list - (g,r,r) cube
    // proj: (p1,p2,g1,g2,r1,r2) = (p1*p2, g1*r1, g2*r2) cubes, with the identification of indices using functions in SBF_comp

    cube tildem = SBF_comp["tildem"];
    List kde_1d_sq = SBF_comp["kde.1d.sq"];
    cube proj = SBF_comp["proj"];
    vec bandwidths = SBF_comp["bandwidths"];
    vec grids = SBF_comp["grids"];
    vec weights = SBF_comp["weights"];

    int p = kde_1d_sq.size();
    int g = weights.size();
    int r = SBF_comp["r"];
    int m = tildem.n_slices;

    // initial values
    cube mhat(p, g * r, m, fill::zeros); // (p, g*r, m)
    cube mhat_old(p, g * r, m, fill::zeros); // (p, g*r, m)
    double A = R;
    double B = 0;

    vec mhat_norm(p, fill::zeros);
    double sigma = phi + eta;

    // ADMM algorithm
    double residual;
    int iter = 0;
    mat tmp_hatm_jj2(g * r, m, fill::zeros);
    mat mhat_j(g * r, m, fill::zeros);
    mat mhat_j2(g * r, m, fill::zeros);
    mat mhat_j_old(g * r, m, fill::zeros);
    while (iter < max_iter) {
        iter = iter + 1;
        residual = 0;
        mhat_old = clone_cube(mhat);

        // mhat update
        for (int j = 0; j < p; j++) {
            mat tmp_hatmj(g * r, m, fill::zeros);
            cube kde_1d_sq_j = kde_1d_sq[j];

            for (int j2 = 0; j2 < p; j2++) {
                mat proj_jj2 = proj.row(p * j + j2);

                if (j2 == j) {
                    if (m == 1) {
                        mat tmp_tmp_hatm_jj2 = tildem.row(j);
                        tmp_hatm_jj2 = tmp_tmp_hatm_jj2.t();
                    }
                    else {
                        tmp_hatm_jj2 = tildem.row(j);
                    }

                    tmp_hatmj += tmp_hatm_jj2;
                }
                else if (mhat_norm(j2) != 0) {
                    if (m == 1) {
                        mat tmp_mhat_j2 = mhat.row(j2);
                        mhat_j2 = tmp_mhat_j2.t();
                    }
                    else {
                        mhat_j2 = mhat.row(j2);
                    }
                    tmp_hatm_jj2 = numerical_integral_2d(proj_jj2, mhat_j2);
                    tmp_hatmj -= tmp_hatm_jj2;
                }
            }


            // compute components
            double norm_hatmj = L2_mat_norm_SBF(tmp_hatmj, kde_1d_sq_j, weights, Ymu, Yspace);
            double nuj = B + eta * (sum(mhat_norm) - mhat_norm(j) + A - R);


            // mhat_j update
            if (m == 1) {
                mat tmp_mhat_j = mhat.row(j);
                mhat_j = tmp_mhat_j.t();
            }
            else {
                mhat_j = mhat.row(j);
            }
            mhat_j = hatmj_update(tmp_hatmj, sigma, norm_hatmj, penalty, lambda, nuj, gamma, mhat_j, mhat_norm(j));
            mhat.row(j) = mhat_j;
            mhat_norm(j) = L2_mat_norm_SBF(mhat_j, kde_1d_sq_j, weights, Ymu, Yspace);


            // compute the L2 norm between mhat_j and mhat_j_old
            if (m == 1) {
                mat tmp_mhat_j_old = mhat_old.row(j);
                mhat_j_old = tmp_mhat_j_old.t();
            }
            else {
                mhat_j_old = mhat_old.row(j);
            }
            double residual_j = L2_mat_norm_SBF(mhat_j - mhat_j_old, kde_1d_sq_j, weights, Ymu, Yspace);
            residual += residual_j;
        }

        // A update
        if (eta != 0) {
            A = max(R - sum(mhat_norm) - B / eta, 0.0);
        }
        else {
            A = 0.0;
        }

        // B update
        B = B + eta * (sum(mhat_norm) + A - R);

        if (abs(residual) < threshold) {
            break;
        }
    }


    // change the dimension of mhat
    // (p, g*r, m) cube -> p list - (g,r,m) cube
    List mhat_final(p);
    for (int j = 0; j < p; j++) {
        if (m == 1) {
            mat tmp_mhat_j = mhat.row(j);
            mhat_j = tmp_mhat_j.t();
        }
        else {
            mhat_j = mhat.row(j);
        }

        cube mhat_final_j(g, r, m, fill::zeros);
        for (int k = 0; k < g; k++) {
            mat mhat_jk = mhat_j.rows(r * k, r * (k + 1) - 1);
            mhat_final_j.row(k) = mhat_jk;
        }

        mhat_final[j] = mhat_final_j;
    }


    List result = List::create(Named("mhat") = mhat_final, Named("mhat.org") = mhat, 
        Named("mhat.norm") = mhat_norm, Named("bandwidths") = bandwidths,
        Named("A") = A, Named("B") = B, Named("iter") = iter, Named("grids") = grids, 
        Named("lambda") = lambda, Named("R") = R, Named("phi") = phi,
        Named("penalty") = penalty, Named("gamma") = gamma);

    return(result);
}


//////////////////////////////////////////////////////////
///////////////// prediction


arma::mat approx_rcpp(arma::vec x, arma::mat y, arma::vec xnew) {
    int g = x.size();
    int n2 = xnew.size();
    int m = y.n_cols;

    mat yout(n2, m);
    for (int i = 0; i < n2; i++) {
        // Find the interval in which xnew[i] lies
        for (int j = 0; j < g - 1; j++) {
            if (xnew[i] >= x[j] && xnew[i] <= x[j + 1]) {
                // Linear interpolation formula
                rowvec slope = (y.row(j + 1) - y.row(j)) / (x[j + 1] - x[j]);
                yout.row(i) = y.row(j) + slope * (xnew[i] - x[j]);
                break;
            }
        }
    }
    return yout;
}



arma::mat predict_AM(List object, arma::mat Xnew) {
    // Xnew: a new transformed score matrix

    List mhat = object["mhat"];
    vec mhat_norm = object["mhat.norm"];
    vec grids = object["grids"];
    int p = mhat.size();

    cube mhat0 = mhat[0];
    int n2 = Xnew.n_rows;
    int m = mhat0.n_slices;

    mat mhatnew(n2, m);
    for (int j = 0; j < p; j++) {
        cube mhat_j = mhat[j];
        vec Xnew_j = Xnew.col(j);
        if (mhat_norm[j] != 0) {
            mhatnew += approx_rcpp(grids, mhat_j.col(0), Xnew_j);
        }
    }

    return mhatnew;
}




//////////////////////////////////////////////////////////
///////////////// cross validation 


// [[Rcpp::export]]
double get_loss_CV_AM_average(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace, double lambda, double R, String cv_type, String penalty, double gamma, double cv_const) {

    // model training
    List model = AM_each(SBF_comp, Ymu, Yspace, lambda, R, penalty, gamma);

    // compute Yhat
    mat LogYhat = predict_AM(model, Xnew);

    // compute mse
    int n = LogYnew.n_rows;
    double loss = L2_norm(LogYhat - LogYnew, Ymu, Yspace) / sqrt(n);
    loss = pow(loss, 2) / 2;
    loss = log(loss);

    // compute additional penalty term for AIC and BIC
    double aic = 0;
    double bic = 0;

    vec mhat_norm = model["mhat.norm"];
    vec bandwidths = model["bandwidths"];
    vec nh = n * bandwidths;

    if (cv_type == "AIC" || cv_type == "ABIC") {
        vec nh_tmp = cv_const / nh;
        aic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //aic = 2.0 * sum(mhat_norm != 0) / n;
    }
    if (cv_type == "BIC" || cv_type == "ABIC") {
        vec nh_tmp = cv_const * log(nh) / (2 * nh);
        bic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //bic = sum(mhat_norm != 0) * log(n) / n;
    }

    loss = loss + aic + bic;

    return(loss);
}



// [[Rcpp::export]]
double get_loss_CV_AM_integral(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace, double lambda, double R, String cv_type, String penalty, double gamma, double cv_const) {

    // tildem: (p, g*r, m) cube
    // kde_1d: p list - (g,r,r) cube
    // proj: (p1,p2,g1,g2,r1,r2) = (p1*p2, g1*r1, g2*r2) cubes, with the identification of indices using functions in SBF_comp
    cube tildem = SBF_comp["tildem"];
    List kde_1d_sq = SBF_comp["kde.1d.sq"];
    cube proj = SBF_comp["proj"];
    vec weights = SBF_comp["weights"];

    int p = tildem.n_rows;
    int g = weights.size();
    int r = SBF_comp["r"];
    int n = LogYnew.n_rows;
    int m = LogYnew.n_cols;

    // model training
    List model = AM_each(SBF_comp, Ymu, Yspace, lambda, R, penalty, gamma);
    
    cube mhat = model["mhat.org"];
    
    mat mhat_j(g * r, m, fill::zeros);
    mat mhat_j2(g * r, m, fill::zeros);
    mat tildem_j(g * r, m, fill::zeros);
    double loss = 0;
    for (int j = 0; j < p; j++) {
        if (m == 1) {
            mat tmp_mhat_j = mhat.row(j);
            mhat_j = tmp_mhat_j.t();
        }
        else {
            mhat_j = mhat.row(j); // (g*r, m) mat
        }
        cube kde_1d_sq_j = kde_1d_sq[j]; // (g,r,r) cube

        // the square norms of hatm_j
        loss += L2_mat_inner_SBF(mhat_j, mhat_j, kde_1d_sq_j, weights, Ymu, Yspace) / 2;
        
        // the inner products of hatm_j and hatm_j2 for j!=j2
        for (int j2 = j + 1; j2 < p; j2++) {
            if (m == 1) {
                mat tmp_mhat_j2 = mhat.row(j2);
                mhat_j2 = tmp_mhat_j2.t();
            }
            else {
                mhat_j2 = mhat.row(j2); // (g*r, m) mat
            }

            mat proj_jj2 = proj.row(p * j + j2);
            mat tmp_mhat_jj2 = numerical_integral_2d(proj_jj2, mhat_j2);
            
            loss += L2_mat_inner_SBF(mhat_j, tmp_mhat_jj2, kde_1d_sq_j, weights, Ymu, Yspace);
        }
        
        // the inner products of hatm_j and tildem_j
        if (m == 1) {
            mat tmp_tildem_j = tildem.row(j);
            tildem_j = tmp_tildem_j.t();
        }
        else {
            tildem_j = tildem.row(j); // (g*r, m) mat
        }
        loss -= L2_mat_inner_SBF(mhat_j, tildem_j, kde_1d_sq_j, weights, Ymu, Yspace);
    }

    // the average square norm of LogY
    double Ysigma = L2_norm(LogYnew, Ymu, Yspace) / sqrt(n);
    loss += pow(Ysigma, 2) / 2;
    loss = log(loss);

    // compute additional penalty term for AIC and BIC
    double aic = 0;
    double bic = 0;

    vec mhat_norm = model["mhat.norm"];
    vec bandwidths = model["bandwidths"];
    vec nh = n * bandwidths;

    if (cv_type == "AIC" || cv_type == "ABIC") {
        vec nh_tmp = cv_const / nh;
        aic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //aic = 2.0 * sum(mhat_norm != 0) / n;
    }
    if (cv_type == "BIC" || cv_type == "ABIC") {
        vec nh_tmp = cv_const * log(nh) / (2 * nh);
        bic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //bic = sum(mhat_norm != 0) * log(n) / n;
    }

    loss = loss + aic + bic;

    return(loss);
}

