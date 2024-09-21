#ifndef PENALTYSOL_H
#define PENALTYSOL_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat betaj_update(mat wj, double sigmaj, double normwj, String penalty, double lambdaj, double kappaj, double gamma, mat betaj, double betaj_norm);
cube hatmj_update(cube tmp_hatmj, double sigma, double norm_hatmj, String penalty, double lambda, double nu, double gamma, cube mhat_j, double norm_mhat_j);

mat LASSO_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat LASSO_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

mat SCAD_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat SCAD_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

mat MCP_sol(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);
mat MCP_sol2(mat xy, double xx, double normxy, double lambda, double kappa=0, double gamma=0);

#endif

