library(manifold)
library(fdapace)
library(MASS)
library(expm)
library(pracma)
library(base)
library(hms)
library(Rcpp)
library(RcppArmadillo)




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


#' vector norm of matrices on the tangent space at \eqn{p}
#' 
#' @param X a \eqn{p\times m} matrix of manifold-valued data.
#' @param p a base point on the tangent space.
#' @param space an underlying space of \eqn{X}.
#' @param type a type of vector norm. One of 'L2', 'L1', or 'Linf'.
#' 
#' @return a vector norm
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





