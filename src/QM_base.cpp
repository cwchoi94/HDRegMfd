// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
//#include "QM_base.h"


using namespace Rcpp;
using namespace std;
using namespace arma;



// Gaussian loss
arma::mat loss_Gaussian(arma::mat u) {
    return sqrt(2.0 / arma::datum::pi) * exp(-pow(u, 2) / 2.0) + u % (1.0 - 2.0 * arma::normcdf(-u));
}

arma::mat loss_diff_Gaussian(arma::mat u) {
    return -sqrt(2.0 / arma::datum::pi) * u % exp(-pow(u, 2) / 2.0) + (1 - 2.0 * arma::normcdf(-u)) + 2.0 * u % arma::normpdf(-u);
}


// Uniform loss
arma::mat loss_Uniform(arma::mat u) {
    mat cond = arma::conv_to<arma::mat>::from(arma::abs(u) < 1);

    mat loss1 = (pow(u, 2) + 1.0) / 2.0;
    mat loss2 = arma::abs(u);

    return loss1 % cond + loss2 % (1 - cond);
}

arma::mat loss_diff_Uniform(arma::mat u) {
    mat cond1 = arma::conv_to<arma::mat>::from(arma::abs(u) < 1);
    mat cond2 = arma::conv_to<arma::mat>::from(u <= -1);
    mat cond3 = arma::conv_to<arma::mat>::from(u >= 1);

    return u % cond1 - 1.0 * cond2 + 1.0 * cond3;
}


// Epanechnikov kernel
arma::mat loss_Epanechnikov(arma::mat u) {
    mat cond = arma::conv_to<arma::mat>::from(arma::abs(u) < 1);

    mat loss1 = 3.0 / 4.0 * pow(u, 2) - pow(u, 4) / 8.0 + 3.0 / 8.0;
    mat loss2 = arma::abs(u);

    return loss1 % cond + loss2 % (1 - cond);
}

arma::mat loss_diff_Epanechnikov(arma::mat u) {
    mat cond1 = arma::conv_to<arma::mat>::from(arma::abs(u) < 1);
    mat cond2 = arma::conv_to<arma::mat>::from(u <= -1);
    mat cond3 = arma::conv_to<arma::mat>::from(u >= 1);

    mat loss = 3.0 / 2.0 * u - pow(u, 3) / 4.0;

    return loss % cond1 - 1 * cond2 + 1 * cond3;
}


// Triangular kernel
arma::mat loss_Triangular(arma::mat u) {
    mat u_abs = arma::abs(u);
    mat cond = arma::conv_to<arma::mat>::from(arma::abs(u) < 1);

    mat loss1 = pow(u_abs, 2) - pow(u_abs, 3) / 3.0 + 1.0 / 3.0;
    mat loss2 = u_abs;

    return loss1 % cond + loss2 % (1 - cond);
}

arma::mat loss_diff_Triangular(arma::mat u) {
    mat cond1 = arma::conv_to<arma::mat>::from(u <= -1);
    mat cond2 = arma::conv_to<arma::mat>::from(u >= 1);
    mat cond3 = arma::conv_to<arma::mat>::from((-1 < u) && (u <= 0));
    mat cond4 = arma::conv_to<arma::mat>::from((0 < u) && (u < 1));

    mat loss3 = 2.0 * u + pow(u, 2);
    mat loss4 = 2.0 * u - pow(u, 2);

    return -1 * cond1 + 1 * cond2 + loss3 % cond3 + loss4 % cond4;
}


/////////////////////////////////////////////
/////////////////////////////////////////////


// smoothed quantile loss ell_{tau,h}(u)
// [[Rcpp::export]]
arma::mat sQRloss(arma::mat u, double tau, double h, String kernel) {
    mat loss;
    if (kernel == "Gaussian") {
        loss = loss_Gaussian(u / h);
    }
    else if (kernel == "Uniform") {
        loss = loss_Uniform(u / h);
    }
    else if (kernel == "Epanechnikov") {
        loss = loss_Epanechnikov(u / h);
    }
    else if (kernel == "Triangular") {
        loss = loss_Triangular(u / h);
    }


    return h / 2.0 * loss + (tau - 1.0 / 2.0) * u;
}


// first derivative of smoothed quantile loss ell_{tau,h}'(u)
// [[Rcpp::export]]
arma::mat sQRloss_diff(arma::mat u, double tau, double h, String kernel) {
    mat loss_diff;
    if (kernel == "Gaussian") {
        loss_diff = loss_diff_Gaussian(u / h);
    }
    else if (kernel == "Uniform") {
        loss_diff = loss_diff_Uniform(u / h);
    }
    else if (kernel == "Epanechnikov") {
        loss_diff = loss_diff_Epanechnikov(u / h);
    }
    else if (kernel == "Triangular") {
        loss_diff = loss_diff_Triangular(u / h);
    }

    return loss_diff / 2.0 + (tau - 1.0 / 2.0);
}

