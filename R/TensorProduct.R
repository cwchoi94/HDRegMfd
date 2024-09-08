


#' @title Make an element in a Tensor Product Space 
#' 
#' @description 
#' Given two Hilbert spaces \eqn{\mathbb{H}_1} and \eqn{\mathbb{H}_2}, make an element in the tensor product space \eqn{H_1\otimes H_2}.
#' 
#' @param element1 an \eqn{n\times m_1} matrix, where each row is an element in \eqn{\mathbb{H}_1}.
#' @param element2 an \eqn{n\times m_2} matrix, where each row is an element in \eqn{\mathbb{H}_2}.
#' @param space1 the name of the underlying space \eqn{\mathbb{H}_1}.
#' @param space2 the name of the underlying space \eqn{\mathbb{H}_2}.
#' @param mu1 a base point of \eqn{\mathbb{H}_1}.
#' @param mu2 a base point of \eqn{\mathbb{H}_2}.
#' 
#' @return a "tensor" object with the following arguments:
#'    \describe{
#'       \item{n}{the number of elements.}
#'       \item{element1}{an \eqn{n\times m_1} matrix, where each row is an element in \eqn{\mathbb{H}_1}.}
#'       \item{element2}{an \eqn{n\times m_2} matrix, where each row is an element in \eqn{\mathbb{H}_2}.}
#'       \item{mu1}{a base point of \eqn{\mathbb{H}_1}.}
#'       \item{mu2}{a base point of \eqn{\mathbb{H}_2}.}
#'       \item{space1}{the name of the underlying space \eqn{\mathbb{H}_1}.}
#'       \item{space2}{the name of the underlying space \eqn{\mathbb{H}_2}.}
#'       \item{dim}{the intrinsic dimension of \eqn{\mathbb{H}_1}.}
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


#' @title Inner product on the Tensor Product Space
#' 
#' @description
#' Computes the inner product of two elements in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}.
#'
#' @param x,y elements in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}, created by \code{\link{make.tensor}}.
#' 
#' @return the inner product of x and y.
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

#' @title Norm on the Tensor Product Space
#' 
#' @description
#' Computes the inner product of an element in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}.
#' 
#' @param x an element in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}, created by \code{\link{make.tensor}}.
#' 
#' @return the norm of x.
#' @export
norm.tensor = function(x){
  return(sqrt(inner.tensor(x,x)))
}



#' @title Distance on the Tensor Product Space
#' 
#' @description
#' Computes the distance between two elements in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}.
#' 
#' @inheritParams inner.tensor
#' 
#' @return the distance between x and y.
#' @export
dist.tensor = function(x,y){
  distsq = inner.tensor(x,x) + inner.tensor(y,y) - 2*inner.tensor(x,y)
  return(sqrt(distsq))
}


#' @title Operator Map from \eqn{\mathbb{H}_1} to \eqn{\mathbb{H}_2}
#' 
#' @description
#' Given a tensor \eqn{x\in\mathbb{H}_1\otimes\mathbb{H}_2} and an element \eqn{y\in\mathbb{H}_1}, this function computes the element \eqn{x(y) \in \mathbb{H}_2}.
#' 
#' @param x an element in \eqn{\mathbb{H}_1\otimes \mathbb{H}_2}, created by \code{\link{make.tensor}}.
#' @param y an \eqn{n\times m_1} matrix where each row is an element in \eqn{\mathbb{H}_1}.
#' 
#' @return an \eqn{n\times m_2} matrix where each row is an element in \eqn{\mathbb{H}_2}.
#' @export
operator.tensor = function(x,y){
  mu1 = x$mu1
  space1 = x$space1
  element1 = x$element1
  
  dim1 = x$dim
  dim2 = nrow(y)
  
  # E: (m2,p) matrix
  E = sapply(1:dim1,function(i){sapply(1:dim2,function(j){inner.manifold(element1[i,],y[j,],mu1,space1)})})
  return (E %*% x$element2)
}









