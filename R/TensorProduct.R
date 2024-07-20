### utils for tensor product of Hilbert spaces


#' @title Make an element in tensor product spaces of \eqn{H_1} and \eqn{H_2}, \eqn{H_1\otimes H_2}
#' 
#' @description 
#' Make an element in tensor product spaces of \eqn{H_1} and \eqn{H_2}, \eqn{H_1\otimes H_2}
#' Let \eqn{d_1} and \eqn{d_2} be the dimensions of \eqn{H_1} and \eqn{H_2}.
#' 
#' @param element1 an \eqn{n}-by-\eqn{m_1} matrix. Each row is an element in \eqn{H_1}.
#' @param element2 an \eqn{n}-by-\eqn{m_2} matrix. Each row is an element in \eqn{H_2}.
#' @param space1,space2 the underlying spaces of \eqn{H_1} and \eqn{H_2}, respectively.
#' @param mu1,mu2 base elements in \eqn{H_1} and \eqn{H_2}, respectively. These will be used for computing norms.
#' 
#' @return a tensor product object, class "tensor"
#'    \describe{
#'       \item{n}{a number of elements \eqn{n}.}
#'       \item{mu1,mu2}{base elements in \eqn{H_1} and \eqn{H_2}, respectively.}
#'       \item{space1,space2}{the underlying space of \eqn{H_1} and \eqn{H_2}, respectively.}
#'       \item{element1}{an \eqn{n}-by-\eqn{m_1} matrix. Each row is an element in \eqn{H_1}.}
#'       \item{element2}{an \eqn{n}-by-\eqn{m_2} matrix. Each row is an element in \eqn{H_2}.}
#'       \item{dim}{the intrinsic input dimension of \eqn{H_1}.}
#' }
#' @export
make.tensor = function(element1,element2,space1,space2,mu1,mu2){
  element1 = vec.to.mat(element1)
  element2 = vec.to.mat(element2)
  
  z = list(mu1=mu1,mu2=mu2,element1=element1,element2=element2,
           space1=space1,space2=space2,dim=nrow(element1))
  class(z) = 'tensor'
  return(z)
}


#' @title Inner product of tensor product spaces
#' 
#' @description
#' Two elements must have the same underlying spaces \eqn{H_1} and \eqn{H_2}, and base elements \eqn{mu_1} and \eqn{mu_2}. 
#'
#' @param x,y elements in \eqn{H_1\otimes H_2} with \eqn{n_1} and \eqn{n_2}, created by \code{\link{make.tensor}}.
#' 
#' @return an inner product of x and y
#' @export
inner.tensor = function(x,y){
  space1 = x$space1
  space2 = x$space2
  mu1 = x$mu1
  mu2 = x$mu2
  
  dim1 = x$dim
  dim2 = y$dim
  
  # E,M: (dim1,dim2) matrices
  # inner product = tr(E %*% t(M))
  E = sapply(1:dim2,function(j){sapply(1:dim1,function(i){inner.manifold(x$element1[i,],y$element1[j,],mu1,space1)})})
  M = sapply(1:dim2,function(j){sapply(1:dim1,function(i){inner.manifold(x$element2[i,],y$element2[j,],mu2,space2)})})
  return(sum(diag(E %*% t(M))))
}

#' Norm of tensor product spaces
#' 
#' @param x an element in \eqn{H_1\otimes H_2}, created by \code{\link{make.TPspace}}.
#' 
#' @return a norm of x
#' @export
norm.tensor = function(x){
  return(sqrt(inner.tensor(x,x)))
}



#' @title Distance of tensor product spaces
#' 
#' @description
#' Two elements must have the same underlying spaces \eqn{H_1} and \eqn{H_2}, and base elements \eqn{mu_1} and \eqn{mu_2}.#' 
#'
#' @param x,y elements in \eqn{H_1\otimes H_2} with \eqn{n_1} and \eqn{n_2}, created by \code{\link{make.tensor}}.
#' 
#' @return a distance of x and y
#' @export
dist.tensor = function(x,y){
  distsq = inner.tensor(x,x) + inner.tensor(y,y) - 2*inner.tensor(x,y)
  return(sqrt(distsq))
}


#' Compute operator maps from \eqn{H_1} to \eqn{H_2}
#' 
#'
#' @param x a tensor product object, created by \code{\link{make.tensor}}.
#' @param y an \eqn{n}-by-\eqn{n_1} matrix. Each row is an element in \eqn{H_1}.
#' 
#' @return an \eqn{n}-by-\eqn{n_2} matrix. Each row is an element in \eqn{H_2}.
#' @export
operator.tensor = function(x,y){
  mu1 = x$mu1
  space1 = x$space1
  element1 = x$element1
  
  dim1 = x$dim
  dim2 = nrow(y)
  
  # E: (n2,p) matrix
  E = sapply(1:dim1,function(i){sapply(1:dim2,function(j){inner.manifold(element1[i,],y[j,],mu1,space1)})})
  return (E %*% x$element2)
}









