// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include "GLM_base.h"


using namespace Rcpp;
using namespace std;
using namespace arma;


// Link function G
// [[Rcpp::export]]
arma::mat Link(arma::mat u, String link) {
    mat z;
    if (link == "binomial") {
        z = log(u / (1 - u));
    }
    else if (link == "poisson") {
        z = log(u);
    }
    else if (link == "exponential") {
        z = -1/u;
    }
    return(z);
}


// Inverse link function G^{-1}
// [[Rcpp::export]]
arma::mat Inv_Link(arma::mat u, String link) {
    mat z;
    if (link == "binomial") {
        z = 1 / (1 + exp(-u));
    }
    else if (link == "poisson") {
        z = exp(u);
    }
    else if (link == "exponential") {
        z = -1/u;
    }
    return(z);
}


// Log partition function Psi
// [[Rcpp::export]]
arma::mat Psi(arma::mat u, String link) {
    mat z;
    if (link == "binomial") {
        z = log(1 + exp(u));
    }
    else if (link == "poisson") {
        z = exp(u);
    }
    else if (link == "exponential") {
        z = -log(-u);
    }
    return(z);
}


// First derivative of log partition function Psi'=G^{-1}
// [[Rcpp::export]]
arma::mat Psi_1d(arma::mat u, String link) {
    mat z = Inv_Link(u, link);
    return(z);
}


// Second derivative of log partition function Psi''=G^{-1}'
// [[Rcpp::export]]
arma::mat Psi_2d(arma::mat u, String link) {
    mat z;
    if (link == "binomial") {
        z = exp(u);
        z = z / square(1 + z);
    }
    else if (link == "poisson") {
        z = exp(u);
    }
    else if (link == "exponential") {
        z = 1 / square(u);
    }
    return(z);
}



