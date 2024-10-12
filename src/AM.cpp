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



// [[Rcpp::export]]
double L2_mat_inner_SBF(arma::cube Yj1, arma::cube Yj2, arma::cube kde_1d_j, arma::vec weights, arma::vec Ymu, String Yspace) {
    // Yj1, Yj2: (g,r,m) cubes

    int g = Yj1.n_rows;
    int m = Yj1.n_slices;

    double z = 0;
    for (int k = 0; k < g; k++) {
        mat Y_jk1 = Yj1.row(k);
        mat Y_jk2 = Yj2.row(k);
        if (m == 1) {
            Y_jk1 = Y_jk1.t(); // (r,m) matrix
            Y_jk2 = Y_jk2.t(); // (r,m) matrix
        }
        mat kde_1d_jk = kde_1d_j.row(k); // (r,r) matrix        
        mat Y_jk1_ = kde_1d_jk * Y_jk1; // (r,m) matrix

        vec z_k = inner(Y_jk1_, Y_jk2, Ymu, Yspace);

        z += weights[k] * sum(z_k);
    }

    return z;
}


// [[Rcpp::export]]
double L2_mat_norm_SBF(arma::cube Yj, arma::cube kde_1d_j, arma::vec weights, arma::vec Ymu, String Yspace) {
    
    double z = L2_mat_inner_SBF(Yj, Yj, kde_1d_j, weights, Ymu, Yspace);
    return sqrt(z);
}



// [[Rcpp::export]]
List AM_each(List SBF_comp, arma::vec Ymu, String Yspace, double lambda, double R, String penalty, double gamma, double phi, double eta, int max_iter, double threshold) {

    // p = the number of (j,k) with k <= Xdim_max
    // tildem: p list - (g,r,m) cube
    // kde_1d: p list - (g,r,r) cube
    // proj: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2, r1, r2) cubes, with the identification of indices using functions in utils_SBF 

    List tildem = SBF_comp["tildem"];
    List kde_1d = SBF_comp["kde.1d"];
    cube proj = SBF_comp["proj"];
    vec bandwidths = SBF_comp["bandwidths"];
    vec grids = SBF_comp["grids"];
    vec weights = SBF_comp["weights"];

    int p = tildem.size();
    int g = weights.size();

    // initial values
    cube tildem0 = tildem[0];
    List mhat(p);
    List mhat_old(p);
    for (int j = 0; j < p; j++) {
        cube mhat_j(size(tildem0), fill::zeros);
        mhat[j] = mhat_j;
    }
    double A = R;
    double B = 0;

    vec mhat_norm(p, fill::zeros);
    double sigma = phi + eta;

    // ADMM algorithm
    double residual;
    int iter = 0;
    while (iter < max_iter) {
        iter = iter + 1;
        residual = 0;
        mhat_old = Rcpp::clone(mhat);

        // mhat update
        for (int j = 0; j < p; j++) {
            cube tmp_hatmj(size(tildem0), fill::zeros);
            cube kde_1d_j = kde_1d[j];

            for (int j2 = 0; j2 < p; j2++) {
                IntegerVector ind_range = multi_2d_ind_to_single_range(j, j2, p, g);
                cube proj_jj2 = proj.rows(ind_range(0), ind_range(1));

                if (j2 == j) {
                    cube tmp_hatm_jj2 = tildem[j];
                    tmp_hatmj += tmp_hatm_jj2;
                }
                else if (mhat_norm(j2) != 0) {
                    cube mhat_j2 = mhat[j2];
                    cube tmp_hatm_jj2 = numerical_integral_3d(weights, proj_jj2, mhat_j2);
                    tmp_hatmj -= tmp_hatm_jj2;
                }
            }

            // compute components
            double norm_hatmj = L2_mat_norm_SBF(tmp_hatmj, kde_1d_j, weights, Ymu, Yspace);
            double nuj = B + eta * (sum(mhat_norm) - mhat_norm(j) + A - R);

            // mhat_j update
            cube mhat_j = mhat[j];
            mhat_j = hatmj_update(tmp_hatmj, sigma, norm_hatmj, penalty, lambda, nuj, gamma, mhat_j, mhat_norm(j));
            mhat[j] = mhat_j;
            mhat_norm(j) = L2_mat_norm_SBF(mhat_j, kde_1d_j, weights, Ymu, Yspace);

            // compute the L2 norm between mhat_j and mhat_j_old
            cube mhat_j_old = mhat_old[j];
            double residual_j = L2_mat_norm_SBF(mhat_j - mhat_j_old, kde_1d_j, weights, Ymu, Yspace);
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

    List result = List::create(Named("mhat") = mhat, Named("mhat.norm") = mhat_norm, Named("bandwidths") = bandwidths, 
        Named("A") = A, Named("B") = B, Named("iter") = iter, Named("grids") = grids, 
        Named("lambda") = lambda, Named("R") = R, Named("phi") = phi,
        Named("penalty") = penalty, Named("gamma") = gamma);

    return(result);
}



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



// [[Rcpp::export]]
double get_loss_CV_AM_average(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace, double lambda, double R, String cv_type, String penalty, double gamma) {

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
        vec nh_tmp = 1.0 / nh;
        aic = 2.0 * sum(nh_tmp.elem(find(mhat_norm != 0)));
        //aic = 2.0 * sum(mhat_norm != 0) / n;
    }
    if (cv_type == "BIC" || cv_type == "ABIC") {
        vec nh_tmp = log(nh) / nh;
        bic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //bic = sum(mhat_norm != 0) * log(n) / n;
    }

    loss = loss + aic + bic;

    return(loss);
}



// [[Rcpp::export]]
double get_loss_CV_AM_integral(List SBF_comp, arma::mat Xnew, arma::mat LogYnew, arma::vec Ymu, String Yspace, double lambda, double R, String cv_type, String penalty, double gamma) {

    List tildem = SBF_comp["tildem"];
    List kde_1d = SBF_comp["kde.1d"];
    cube proj = SBF_comp["proj"];
    vec weights = SBF_comp["weights"];

    int p = tildem.size();
    int g = weights.size();
    int n = LogYnew.n_rows;

    // model training
    List model = AM_each(SBF_comp, Ymu, Yspace, lambda, R, penalty, gamma);
    
    List mhat = model["mhat"];    
    

    double loss = 0;
    for (int j = 0; j < p; j++) {
        cube mhat_j = mhat[j];
        cube kde_1d_j = kde_1d[j];

        // the square norms of hatm_j
        loss += L2_mat_inner_SBF(mhat_j, mhat_j, kde_1d_j, weights, Ymu, Yspace) / 2;

        // the inner products of hatm_j and hatm_j2 for j!=j2
        for (int j2 = j + 1; j2 < p; j2++) {
            cube mhat_j2 = mhat[j2];

            IntegerVector ind_range = multi_2d_ind_to_single_range(j, j2, p, g);
            cube proj_jj2 = proj.rows(ind_range(0), ind_range(1));
            cube tmp_mhat_jj2 = numerical_integral_3d(weights, proj_jj2, mhat_j2);

            loss += L2_mat_inner_SBF(mhat_j, tmp_mhat_jj2, kde_1d_j, weights, Ymu, Yspace);
        }

        // the inner products of hatm_j and tildem_j
        cube tildem_j = tildem[j];
        loss -= L2_mat_inner_SBF(mhat_j, tildem_j, kde_1d_j, weights, Ymu, Yspace);
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
        vec nh_tmp = 1.0 / nh;
        aic = 2.0 * sum(nh_tmp.elem(find(mhat_norm != 0)));
        //aic = 2.0 * sum(mhat_norm != 0) / n;
    }
    if (cv_type == "BIC" || cv_type == "ABIC") {
        vec nh_tmp = log(nh) / nh;
        bic = sum(nh_tmp.elem(find(mhat_norm != 0)));
        //bic = sum(mhat_norm != 0) * log(n) / n;
    }

    loss = loss + aic + bic;

    return(loss);
}

