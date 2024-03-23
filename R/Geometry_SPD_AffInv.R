

#' Inner product on tangent space at p for SPD data with the Affine-invariant metric

#' @describeIn inner.manifold Method
#' @export
inner.SPD.AffInv = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  
  
  z = rowSums(u * v)
  return(z)
}

inner.each.SPD.AffInv = function(u,v,p=NULL){
  z = u %*% v
  return(z)
}


#' norm on tangent space at p for SPD data with the Affine-invariant metric

#' @describeIn norm.manifold Method
#' @export
norm.SPD.AffInv = function(u,p=NULL){
  z = sqrt(inner.SPD.AffInv(u,u,p))
  return(z)
}


#' Geodesic distance for SPD data with the Affine-invariant metric

#' @describeIn dist.manifold Method
#' @export
dist.SPD.AffInv = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = manifold::distance(manifold::createM('AffInv'),t(p),t(q))
  return(z)
}


#' Riemannian exponential map for SPD data with the Affine-invariant metric

#' @describeIn RieExp.manifold Method
#' @export
RieExp.SPD.AffInv = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  z = manifold::rieExp(manifold::createM('AffInv'),t(p),t(u))
  z = t(z)
  return(z)
}


#' Riemannian logarithmic map for SPD data with the Affine-invariant metric

#' @describeIn RieExp.manifold Method
#' @export
RieLog.SPD.AffInv = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = manifold::rieLog(manifold::createM('AffInv'),t(p),t(q))
  z = t(z)
  return(z)
}


#' Basis on the tangent space at the point for SPD data with the Affine-invariant metric

#' @describeIn Basis.manifold Method
#' @export
basis.SPD.AffInv = function(p,dim=50){
  m = length(p)
  dim = sqrt(m)
  p = matrix(p,dim,dim)
  
  ind = which(lower.tri(p, diag=TRUE), arr.ind=TRUE, useNames=FALSE)
  M = matrix(0,dim,dim)
  z = apply(ind, 1, function(ii) {
    M <- matrix(0, dim, dim)
    if (ii[1] == ii[2]) {
      M[ii[1], ii[2]] = 1
    } else {
      M[ii[1], ii[2]] = 1/sqrt(2)
      M[ii[2], ii[1]] = 1/sqrt(2)
    }
    M
  })
  z = t(z)
  
  return(z)
}


#' Frechet mean for SPD.AffInv data

#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.SPD.AffInv = function(X){
  X = vec.to.mat(X)
  
  mu = manifold::frechetMean(manifold::createM('AffInv'),t(X))
  mu = t(mu)[1,]
  return(mu)
}


#' Principal component analysis for SPD data with the Affine-invariant metric

#' @describeIn PCA.manifold Method
#' @export
PCA.SPD.AffInv = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  dim = sqrt(m)
  k = dim*(dim+1)/2
  
  mu = FrechetMean.SPD.AffInv(X)
  X.tangent = RieLog.SPD.AffInv(mu,X)
  
  Cov = t(X.tangent) %*% X.tangent / n
  C = eigen(Cov)
  
  values = C$values[1:k]
  vectors = t(C$vectors)[1:k,]
  
  result = list(values=values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' Prediction score matrix for SPD data with the Affine-invariant metric

#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.SPD.AffInv = function(object,Xnew){
  object[['space']] = 'SPD.AffInv'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}

