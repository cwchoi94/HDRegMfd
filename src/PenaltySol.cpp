// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "PenaltySol.h"

using namespace Rcpp;
using namespace arma;


// betaj update for each iteration
arma::mat betaj_update(arma::mat wj, double sigmaj, double normwj, String penalty, double lambdaj, double nuj, double gamma, mat betaj, double betaj_norm) {
    mat betaj_new;
    bool is_sol2 = (nuj + lambdaj < 0) & (normwj == 0);

    if (penalty == "LASSO") {
        if (is_sol2) {
            betaj_new = LASSO_sol2(betaj, sigmaj, betaj_norm, lambdaj, nuj, gamma);
        }
        else {
            betaj_new = LASSO_sol(wj, sigmaj, normwj, lambdaj, nuj, gamma);
        }
    }
    else if (penalty == "SCAD") {
        if (is_sol2) {
            betaj_new = SCAD_sol2(betaj, sigmaj, betaj_norm, lambdaj, nuj, gamma);
        }
        else {
            betaj_new = SCAD_sol(wj, sigmaj, normwj, lambdaj, nuj, gamma);
        }
    }
    else if (penalty == "MCP") {
        if (is_sol2) {
            betaj_new = MCP_sol2(betaj, sigmaj, betaj_norm, lambdaj, nuj, gamma);
        }
        else {
            betaj_new = MCP_sol(wj, sigmaj, normwj, lambdaj, nuj, gamma);
        }
    }

    return(betaj_new);
}






// basic coordinatewise minimizer function
arma::mat S(arma::mat xy, double xx, double normxy, double lambda){
    double c =  0.0;  
    if (xx>0){
        c = max(NumericVector::create(0,1-lambda/normxy))/xx;
    }       
  
    mat z = c * xy;
    return(z);
}


// LASSO
// [[Rcpp::export]]
arma::mat LASSO_sol(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z = S(xy,xx,normxy,lambda+kappa);

    return(z);
}

// [[Rcpp::export]]
arma::mat LASSO_sol2(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z;

    if (normxy==0){
        z = xy;
    }       
    else{
        z = -(kappa+lambda)/(xx*normxy) * xy;
    }    
    
    return(z);
}


// SCAD
//[[Rcpp::export]]
arma::mat SCAD_sol(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z;

    if (normxy<=(1+xx)*lambda+kappa){
        z = S(xy,xx,normxy,lambda+kappa);
    }        
    else if (normxy<=gamma*xx*lambda+kappa){
        z = ((gamma-1)*xx)/((gamma-1)*xx-1) * S(xy,xx,normxy,lambda*gamma/(gamma-1)+kappa);
    }        
    else{
        z = S(xy,xx,normxy,kappa);
    }        

    return(z);
}

//[[Rcpp::export]]
arma::mat SCAD_sol2(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z;

    if (normxy==0){
        z = xy;
    }
    else if ((1+xx)*lambda+kappa>=0){
        z = -(lambda+kappa)/(xx*normxy) * xy;
    }
    else if (gamma*xx*lambda+kappa<0){
        z = -(gamma*lambda+(gamma-1)*kappa)/(((gamma-1)*xx-1)*normxy) * xy;
    }   
    else{
        z = -kappa/(xx*normxy)*xy;
    }        

    return(z);
}

// MCP
//[[Rcpp::export]]
arma::mat MCP_sol(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z;

    if (normxy<=gamma*xx*lambda+kappa){
        z = (gamma*xx)/(gamma*xx-1)*S(xy,xx,normxy,lambda+kappa);
    }        
    else{
        z = S(xy,xx,normxy,kappa);
    }
        

    return(z);
}

//[[Rcpp::export]]
arma::mat MCP_sol2(arma::mat xy, double xx, double normxy, double lambda, double kappa, double gamma){
    mat z;

    if (normxy==0){
        z = xy;
    }
    else if (gamma*xx*lambda+kappa>=0){
        z = -gamma*(lambda+kappa)/((gamma*xx-1)*normxy) * xy;
    }        
    else{
        z = -kappa/(xx*normxy)*xy;
    }       

    return(z);
}


