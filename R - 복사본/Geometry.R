#' Basic functions for manifolds
#' 
#' Basic geometry structure:
#' inner product, norm in tangent spaces, distance
#' exponential and logarithmic maps
#' 
#' For PCA:
#' Frechet mean, PCA
#' 


#' @title Check if the given space is a manifold
#' 
#' @description
#' Check if the given space is a manifold.
#' The space name must be one of 'Euclid', 'simplex', 'sphere', 'SPD.LogEuclid', 'SPD.AffInv', 'functional', 'BayesHilbert', or 'Wasserstein'.
#' 
#' @param space the name of space, which must be one of 'Euclid', 'simplex', 'sphere', 'SPD.LogEuclid', 'SPD.AffInv', 'functional', 'BayesHilbert', or 'Wasserstein'.
#' 
#' @export
Check.manifold = function(space){
  if (!(space %in% c('Euclid','simplex','sphere','SPD.LogEuclid','SPD.AffInv','functional','BayesHilbert','Wasserstein'))){
    stop("The space must be one of 'Euclid', 'simplex', 'sphere', 'SPD.LogEuclid', 'SPD.AffInv', 'functional', 'BayesHilbert', or 'Wasserstein'.")
  }
}


#' @title Inner product on the tangent space 
#' 
#' @description
#' Compute the inner product on the tangent space \eqn{T_p\mathcal{M}} of the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param u,v \eqn{n\times m'} matrices, where each row is a tangent vector on \eqn{T_p\mathcal{M}}.
#' @param p an \eqn{m} vector representing the base point in \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n} vector of inner products for each corresponding row of \eqn{u} and \eqn{v}.
#' @export
inner.manifold = function(u,v,p,space='Euclid'){
  Check.manifold(space)
  inner = eval(parse(text=paste0('inner.',space)))
  z = inner(u,v,p)
  return(z)
}


#' @title Norm on the tangent space
#' 
#' @description
#' Compute the norm on the tangent space \eqn{T_p\mathcal{M}} of the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param u an \eqn{n\times m'} matrix, where each row is a tangent vector on \eqn{T_p\mathcal{M}}.
#' @param p an \eqn{m} vector representing the base point in \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n} vector of norms for each row of \eqn{u}.
#' @export
norm.manifold = function(u,p,space='Euclid'){
  Check.manifold(space)
  norm = eval(parse(text=paste0('norm.',space)))
  z = norm(u,p)
  return(z)
}



#' @title Geodesic distance on the manifold
#' 
#' @description
#' Compute the geodesic distance on the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param p,q \eqn{m} vectors or \eqn{n\times m} matrices, where each row is a point on \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n} vector of geodesic distances between corresponding rows of \eqn{p} and \eqn{q}.
#' @export
dist.manifold = function(p,q,space='Euclid'){
  Check.manifold(space)
  dist = eval(parse(text=paste0('dist.',space)))
  z = dist(p,q)
  return(z)
}


#' @title Riemannian exponential map on the manifold
#' 
#' @description
#' Compute the Riemnnain exponential map on the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param p an \eqn{m} vector representing the base point in \eqn{\mathcal{M}}.
#' @param u an \eqn{n\times m'} matrix, where each row is a tangent vector on \eqn{T_p\mathcal{M}}.
#' 
#' @return an \eqn{n\times m} matrix of Riemannian exponentials of each row of \eqn{u}.
#' @export
RieExp.manifold = function(p,u,space='Euclid'){
  Check.manifold(space)
  RieExp = eval(parse(text=paste0('RieExp.',space)))
  z = RieExp(p,u)
  return(z)
}


#' @title Riemannian logarithmic map on the manifold
#' 
#' @description
#' Compute the Riemnnain logarithmic map on the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param p,q \eqn{m} vectors or \eqn{n\times m} matrices, where each row is a point on \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n\times m} matrix of Riemannian logarithmics from  corresponding rows of \eqn{p} to those of \eqn{q}.
#' @export
RieLog.manifold = function(p,q,space='Euclid'){
  Check.manifold(space)
  RieLog = eval(parse(text=paste0('RieLog.',space)))
  z = RieLog(p,q)
  return(z)
}


#' @title Orthonormal basis on the tangent space of the manifold
#' 
#' @description
#' Compute the orthonormal basis of the tangent space \eqn{T_p\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param p an \eqn{m} vector representing the base point in \eqn{\mathcal{M}}.
#' @param dim the number of orthonormal basis, only used for infinite-dimensional \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n\times m} matrix where each row is an orthonormal basis of the tangent space \eqn{T_p\mathcal{M}}.
#' @export
basis.manifold = function(p,dim=50,space='Euclid'){
  Check.manifold(space)
  basis = eval(parse(text=paste0('basis.',space)))
  z = basis(p,dim)
  return(z)
}


#' @title The Frechet mean on the manifold
#' 
#' @description
#' Compute the Frechet mean of a random variable \eqn{X} taking values in the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams Check.manifold
#' 
#' @param X an \eqn{n\times m} matrix where each row \eqn{X_i} is a point in \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{m} vector of the Frechet mean \eqn{\mu} of \eqn{X}.
#' @export
FrechetMean.manifold = function(X,space='Euclid'){
  Check.manifold(space)
  FrechetMean = eval(parse(text=paste0('FrechetMean.',space)))
  z = FrechetMean(X)
  return(z)
}


#' @title Principal component analysis for manifold-valued data.
#' 
#' @description
#' Performs the Principal component analysis (spectral decomposition) for a random variable \eqn{X} taking values in the manifold \eqn{\mathcal{M}}.
#' 
#' @inheritParams FrechetMean.manifold
#' 
#' @param alpha a truncation parameter of the number of basis vectors, used only for infinite dimensional \eqn{\mathcal{M}}. Selects the first index where the cumulative variance is equal to or greater than \eqn{\alpha} of the total variance.
#' 
#' @return a 'PCA.manifold' object with the following components:
#' \describe{
#'       \item{space}{the name of the underlying space \eqn{\mathcal{M}}.}
#'       \item{values}{an \eqn{m} vector of eigenvalues.}
#'       \item{vectors}{a \eqn{K\times m} matrix where each row is a corresponding orthonormal basis.}
#'       \item{mu}{the Frechet mean of \eqn{X}.}
#'       \item{dim}{a number of eigenvectors \eqn{K}}
#'       \item{...}{additional arguments passed into specific methods.}
#' }
#' @export
PCA.manifold = function(X,space='Euclid',alpha=0.95){
  Check.manifold(space)
  PCA = eval(parse(text=paste0('PCA.',space)))
  z = PCA(X,alpha)
  z[['space']] = space
  class(z) = 'PCA.manifold'
  return(z)
}


#' @title Prediction of a score matrix of manifold-valued data.
#' 
#' @description
#' Compute a score matrix of the manifold-valued data using the principal component analysis (spectral decomposition) results made by \code{\link{PCA.manaifold}}.
#' 
#' @param object a \code{\link{PCA.manifold}} object.
#' @param Xnew an \eqn{n\times m} matrix where each row is a point in \eqn{\mathcal{M}}.
#' 
#' @return an \eqn{n\times m'} score matrix.
#' @export
predict.PCA.manifold = function(object,Xnew){
  n2 = nrow(Xnew)
  
  space = object[['space']]
  dim = object[['dim']]
  mu = object[['mu']]
  vectors = object[['vectors']]
  
  Xnew.tangent = RieLog.manifold(mu,Xnew,space)
  scores = sapply(1:dim,function(j){inner.manifold(Xnew.tangent,vectors[j,],mu,space)})
  return (scores)
}



