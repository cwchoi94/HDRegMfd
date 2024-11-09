// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include "utils.h"
#include "SBF_utils.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
// multi index - single index pairing



// [[Rcpp::export]]
int multi_4d_ind_to_single(int j1, int j2, int k1, int k2, int p, int g) {
    // for 0 <= j1, j2 < p and 0 <= k1, k2 < g, consider the identification (j1,j2,k1,k2) = j3 for some 0 <= j3 < p^2g^2
    // 
    // compute the corresponding index j3

    int j3 = j1 * p * g * g + j2 * g * g + k1 * g + k2;
    return j3;
}


// [[Rcpp::export]]
IntegerVector single_ind_to_multi_4d(int j3, int p, int g) {
    // for 0 <= j1, j2 < p and 0 <= k1, k2 < g, consider the identification (j1,j2,k1,k2) = j3 for some 0 <= j3 < p^2g^2
    // 
    // compute the corresponding index j3

    int j1, j2, k1, k2;

    // Extract j1 from j3
    j1 = j3 / (p * g * g);
    j3 = j3 % (p * g * g);

    // Extract j2
    j2 = j3 / (g * g);
    j3 = j3 % (g * g);

    // Extract k1
    k1 = j3 / g;

    // Extract k2
    k2 = j3 % g;

    return IntegerVector::create(j1, j2, k1, k2);
}


// [[Rcpp::export]]
IntegerVector multi_2d_ind_to_single_range(int j1, int j2, int p, int g) {
    // for 0 <= j1, j2 < p and 0 <= k1, k2 < g, consider the identification (j1,j2,k1,k2) = j3 for some 0 <= j3 < p^2g^2
    // 
    // find the range of (j1,j2,0,0) ~ (j1,j2,g-1,g-1)

    int j3 = j1 * p * g * g + j2 * g * g;
    int j4 = j3 + g * g - 1;

    return IntegerVector::create(j3, j4);
}

// [[Rcpp::export]]
IntegerVector multi_3d_ind_to_single_range(int j1, int j2, int k1, int p, int g) {
    // for 0 <= j1, j2 < p and 0 <= k1, k2 < g, consider the identification (j1,j2,k1,k2) = j3 for some 0 <= j3 < p^2g^2
    // 
    // find the range of (j1,j2,k1,0) ~ (j1,j2,k1,g-1)

    int j3 = j1 * p * g * g + j2 * g * g + k1 * g;
    int j4 = j3 + g - 1;

    return IntegerVector::create(j3, j4);
}




///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
// numerical integral



// [[Rcpp::export]]
double numerical_integral_1d(arma::vec weights, arma::vec func) {
    // weights, func: g vectors
    // compute 'func' on 'grids'

    return arma::dot(weights, func);
}



// [[Rcpp::export]]
arma::mat numerical_integral_3d_to_2d(arma::vec weights, arma::cube func) {
    // weights: g vector
    // func: (g,r,r') matrix
    // 
    // return:
    // - compute the integral: (r,r') matrix

    int g = weights.size();
    mat tmp = func.row(0);

    mat integral(size(tmp));
    for (int i = 0; i < g; ++i) {
        mat func_each = func.row(i);
        integral += weights[i] * func_each;
    }

    return integral;
}



// [[Rcpp::export]]
arma::cube numerical_integral_3d(arma::vec weights, arma::cube func1, arma::cube func2) {
    // weights: g2 vector
    // func1: (g1,g2,r1,r2) = (g1*g2,r1,r2) cube
    // func2: (g2,r2,r3) cube
    // 
    // return:
    // - compute the integral: (g1,r1,r3) cube

    int g2 = weights.size();
    int g1 = func1.n_rows / g2;
    int r1 = func1.n_cols;
    int r3 = func2.n_slices;

    cube integral(g1, r1, r3);
    for (int i1 = 0; i1 < g1; ++i1) {
        for (int i2 = 0; i2 < g2; ++i2) {
            int ind = g2 * i1 + i2;

            mat func1_each = func1.row(ind);
            mat func2_each = func2.row(i2);
            if (r3 == 1) {
                func2_each = func2_each.t();
            }            
            integral.row(i1) += weights[i2] * (func1_each * func2_each);
        }
    }

    return integral;
}




// [[Rcpp::export]]
arma::mat numerical_integral_2d(arma::mat func1, arma::mat func2) {
    // func1: (g1*r1, g2*r2) mat
    // func2: (g2*r2, p3) mat
    // 
    // return:
    // - compute the integral: (g1*r1, p3) mat

    mat integral = func1 * func2;

    return integral;
}









