// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include <stdio.h>

using namespace Rcpp;
using namespace std;
using namespace arma;



// [[Rcpp::export]]
List PCA_list(List Xall){
        
    int p = Xall["p"];
    StringVector spaces = Xall["spaces"];

    // use PCA.manifold function defined in R
    Environment myEnv = Environment::global_env();
    // Environment myEnv = Environment::namespace_env("HilbertHLR");
    Function pca_ftn = myEnv["PCA.manifold"];
    
    // compute hpca for each Xj
    List pca(p);
    List pca_each;
    std::string space;
    for (int j=0; j<p; j++){
        space = spaces(j);
        pca_each = pca_ftn(Xall[j],space);
        pca[j] = pca_each;
    }

    //nuisance parameters will be written in R
    pca["p"] = p;
    pca["spaces"] = spaces;

    return(pca);
}

// [[Rcpp::export]]
List predict_PCA_list(List pca, List Xnew){

    int p = pca["p"];
    
    // use predict.PCA.manfold function defined in R
    Environment myEnv = Environment::global_env();
    // Environment myEnv = Environment::namespace_env("HilbertHLR");
    Function pca_predict = myEnv["predict.PCA.manifold"];

    List scores(p);
    List pca_each;
    SEXP score_each;
    for (int j=0; j<p; j++){
        score_each = pca_predict(pca[j],Xnew[j]);
        scores[j] = score_each;
    }
        
    return(scores);
}