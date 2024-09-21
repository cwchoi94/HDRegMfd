library(manifold)
library(fdapace)
library(MASS)
library(expm)
library(pracma)
library(base)
library(plyr)
library(hms)
library(Rcpp)
library(RcppArmadillo)



# Check a penalty function.
Check.penalty = function(penalty){
  if (!(penalty %in% c('LASSO','SCAD','MCP'))){
    stop("The penalty must be one of 'LASSO','SCAD','MCP'")
  }
}

# Check a link function.
Check.link = function(link){
  if (!(link %in% c('binomial','poisson','exponential'))){
    stop("The link should be one of 'binomial','poisson' or 'exponential'. If you use an 'identity' or 'normal' link, please use the 'LM' function.")
  }
}

# Check a cv type.
Check.cv.type = function(cv.type){
  if (!(cv.type %in% c('AIC','BIC','ABIC'))){
    stop("The cv.type should be one of 'AIC','BIC' or 'ABIC'.")
  }
}



# if x is a p-vector, change x to (1,p) matrix
vec.to.mat = function(x){
  if (is.null(nrow(x))){
    x = t(as.matrix(x))
  }
  return(x)
}


# Make a matrix whose each row is x.
vec.duplicate = function(x,n=1){
  x = vec.to.mat(x)
  
  if (nrow(x)==1){
    x = matrix(x,n,ncol(x),byrow=TRUE)
  }
  return(x)
}


#' @title Vector Norm on the Tangent Space
#' 
#' @description
#' Computes the vector norm on the tangent space \eqn{T_p\mathcal{M}}.
#' This function is implemented for \eqn{\ell^2}-, \eqn{\ell^1}- and \eqn{\ell^\infty}-types of norms.
#' 
#' @param X a \eqn{p\times m} matrix of manifold-valued data.
#' @param p a base point of the tangent space \eqn{T_p\mathcal{M}}.
#' @param space the name of the underlying space \eqn{\mathcal{M}} of \eqn{X}.
#' @param type the type of vector norm. One of 'L2', 'L1', or 'Linf'.
#' 
#' @return the vector norm
#' @export
vector.norm = function(X,p,space='Euclid',type='L2'){
  if (!(type %in% c('L2','L1','Linf'))){
    stop("'type' should be one of 'L2', 'L1', or 'Linf'")
  }
  X = vec.to.mat(X)
  vec = norm.manifold(X,p,space)
  if (type=='L2'){
    vec = sqrt(sum(vec^2))
  } else if (type=='L1'){
    vec = sum(vec)
  } else if (type=='Linf'){
    vec = max(vec)
  }
  
  return(vec)
}



