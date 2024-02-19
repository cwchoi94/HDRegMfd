



#' Inner product on tangent space at p for Euclid data

#' @describeIn inner.manifold Method
#' @export
inner.Euclid = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowSums(u * v)
  return(z)
}

inner.each.Euclid = function(u,v,p=NULL){
  z = u %*% v
  return(z)
}


#' norm on tangent space at p for Euclid data

#' @describeIn norm.manifold Method
#' @export
norm.Euclid = function(u,p=NULL){
  z = sqrt(inner.Euclid(u,u,p))
  return(z)
}


#' Geodesic distance for Euclid data

#' @describeIn dist.manifold Method
#' @export
dist.Euclid = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = sqrt(rowSums((p-q)^2))
  return(z)
}


#' Riemannian exponential map for Euclid data

#' @describeIn RieExp.manifold Method
#' @export
RieExp.Euclid = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  z = p+u
  return(z)
}


#' Riemannian logarithmic map for Euclid data

#' @describeIn RieExp.manifold Method
#' @export
RieLog.Euclid = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = q-p
  return(z)
}


#' Basis on the tangent space at the point for Euclid data

#' @describeIn basis.manifold Method
#' @export
basis.Euclid = function(p,dim=NULL){
  z = diag(length(p))
  return(z)
}


#' Frechet mean for Euclid data

#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.Euclid = function(X){
  X = vec.to.mat(X)
  
  mu = colMeans(X)
  return(mu)
}


#' Principal component analysis for Euclid data

#' @describeIn PCA.manifold Method
#' @export
PCA.Euclid = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  
  mu = FrechetMean.Euclid(X)
  X.tangent = RieLog.Euclid(mu,X)
  
  Cov = t(X.tangent) %*% X.tangent / n
  C = eigen(Cov)
  
  vectors = t(C$vectors)
  
  result = list(values=C$values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' Prediction score matrix for Euclid data

#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.Euclid = function(object,Xnew){
  object[['space']] = 'Euclid'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}



