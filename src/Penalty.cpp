// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include "Penalty.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



double SCAD(double u, double lambda, double gamma) {
    u = std::abs(u);

    if (u <= lambda) {
        return lambda * u;
    }
    else if (u <= gamma * lambda) {
        return (2.0 * gamma * lambda * u - u * u - lambda * lambda) / (2.0 * (gamma - 1));
    }
    else {
        return lambda * lambda * (gamma + 1) / 2.0;
    }
}

double SCAD_diff(double u, double lambda, double gamma) {
    u = std::abs(u);

    if (u <= lambda) {
        return lambda;
    }
    else if (u <= gamma * lambda) {
        return (gamma * lambda - u) / (gamma - 1);
    }
    else {
        return 0.0;
    }
}


double MCP(double u, double lambda, double gamma) {
    u = std::abs(u);

    if (u <= gamma *lambda) {
        return lambda * u - u * u / (2.0 * gamma);
    }
    else {
        return gamma * lambda * lambda / 2.0;
    }
}

double MCP_diff(double u, double lambda, double gamma) {
    u = std::abs(u);

    if (u <= gamma * lambda) {
        return lambda - u / gamma;
    }
    else {
        return 0.0;
    }
}



////////////////////////////////////////
////////////////////////////////////////




// sum of P_{lambda_j}(u_j)
// [[Rcpp::export]]
double Penalty_ftn(arma::vec u, arma::vec lambda_vec, String penalty, double gamma) {

    int p = u.n_elem;
    vec penalty_vec(p);

    if (penalty == "LASSO") {
        penalty_vec = lambda_vec % u;
    }
    else if (penalty == "SCAD") {
        for (int j = 0; j < p; j++) {
            penalty_vec[j] = SCAD(u[j], lambda_vec[j], gamma);
        }
    }
    else if (penalty == "MCP") {
        for (int j = 0; j < p; j++) {
            penalty_vec[j] = MCP(u[j], lambda_vec[j], gamma);
        }
    }

    return sum(penalty_vec);
}


// a vector of P_{lambda_j}'(u_j)
// [[Rcpp::export]]
arma::vec Penalty_diff_ftn(arma::vec u, arma::vec lambda_vec, String penalty, double gamma) {

    int p = u.n_elem;
    vec penalty_diff_vec(p);

    if (penalty == "LASSO") {
        penalty_diff_vec = lambda_vec;
    }
    else if (penalty == "SCAD") {
        for (int j = 0; j < p; j++) {
            penalty_diff_vec[j] = SCAD_diff(u[j], lambda_vec[j], gamma);
        }
    }
    else if (penalty == "MCP") {
        for (int j = 0; j < p; j++) {
            penalty_diff_vec[j] = MCP_diff(u[j], lambda_vec[j], gamma);
        }
    }

    return penalty_diff_vec;
}



