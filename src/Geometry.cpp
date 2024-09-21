// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <math.h>
#include "Geometry.h"

using namespace Rcpp;
using namespace std;
using namespace arma;



/////// inner product for each manifold


vec inner_Euclid(mat u, mat v, vec p) {
	int n = u.n_rows;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z;
}


vec inner_simplex(mat u, mat v, vec p) {
	int n = u.n_rows;
	int m = u.n_cols;

	vec z(n);
	for (int i = 0; i < n; i++) {
		mat u_mat1 = repmat(u.row(i), m, 1);
		mat u_mat2 = repmat(u.row(i).t(), 1, m);
		mat v_mat1 = repmat(v.row(i), m, 1);
		mat v_mat2 = repmat(v.row(i).t(), 1, m);

		z[i] = accu((u_mat1-u_mat2) % (v_mat1 - v_mat2));
	}
	z /= 2 * m;

	return z;
}


vec inner_functional(mat u, mat v, vec p) {
	int n = u.n_rows;
	int m = u.n_cols;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z / m;
}


vec inner_BayesHilbert(mat u, mat v, vec p) {
	int n = u.n_rows;
	int m = u.n_cols;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z / m;
}


vec inner_sphere(mat u, mat v, vec p) {
	int n = u.n_rows;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z;
}


vec inner_SPD_LogEuclid(mat u, mat v, vec p) {
	int n = u.n_rows;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z;
}


vec inner_SPD_AffInv(mat u, mat v, vec p) {
	int n = u.n_rows;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z;
}


vec inner_Wasserstein(mat u, mat v, vec p) {
	int n = u.n_rows;
	int m = u.n_cols;

	vec z(n);
	for (int i = 0; i < n; i++) {
		z[i] = arma::dot(u.row(i), v.row(i));
	}

	return z / m;
}




////////// general inner product function


vec inner(mat u, mat v, vec p, String space) {
	if (space == "Euclid") {
		return inner_Euclid(u, v, p);
	}
	else if (space == "simplex") {
		return inner_simplex(u, v, p);
	}
	else if (space == "functional") {
		return inner_functional(u, v, p);
	}
	else if (space == "BayesHilbert") {
		return inner_BayesHilbert(u, v, p);
	}
	else if (space == "sphere") {
		return inner_sphere(u, v, p);
	}
	else if (space == "SPD.LogEuclid") {
		return inner_SPD_LogEuclid(u, v, p);
	}
	else if (space == "SPD.AffInv") {
		return inner_SPD_AffInv(u, v, p);
	}
	else {
		return inner_Wasserstein(u, v, p);
	}
}



