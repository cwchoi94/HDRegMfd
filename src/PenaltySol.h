#ifndef PENALTYSOL_H
#define PENALTYSOL_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat LASSO_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat LASSO_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

mat SCAD_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat SCAD_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

mat MCP_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat MCP_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

#endif

