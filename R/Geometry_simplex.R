
#' @title Centered log-ratio transform
#' 
#' @description 
#' Transform simplex or density data to R^m
#' clr(x)_i = log(x_i)-sum_j log(x_j)
#' 
#' @param x an (n,m) matrix of simplex-valued or Bayes-Hilbert space-valued data
#' 
#' @return A clr transformed (n,m) matrix
#' @export
clr = function(x){
  x = vec.to.mat(x)
  y = log(x)
  y = y - rowMeans(y)
  return (y)
}


#' @title Inverse centered log-transform for simplex data
#' 
#' @description 
#' Inverse transform R^m to simplex data
#' inv.clr(x)_i = exp(x_i)/(sum_j exp(x_i)) * normalizer
#' 
#' @param x A clr transformed (n,m) matrix
#' @param normalizer A normalizing constant, default=1 
#' 
#' @return A simplex-valued (n,m) matrix
#' @export
inv.clr.simplex = function(x,normalizer=1){
  x = vec.to.mat(x)
  x = x - rowMeans(x)
  y = exp(x)
  y = y/rowSums(y)*normalizer
  return(y)
}



#' Inner product on tangent space at p for simplex data
#' 
#' @describIn inner.manifold Method
#' @export
inner.simplex = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  n = nrow(u)
  z = sapply(1:n,function(i){inner.each.simplex(u[i,],v[i,])})
  return(z)
}

inner.each.simplex = function(u,v){
  m = length(u)
  a = outer(u,u,'-')
  b = outer(v,v,'-')
  z = sum(a*b)/(2*m)
  return(z)
}


#' norm on tangent space at p for simplex data
#' 
#' @describIn norm.manifold Method
#' @export
norm.simplex = function(u,p=NULL){
  z = sqrt(inner.simplex(u,u,p))
  return(z)
}


#' Geodesic distance for simplex data
#' 
#' @describIn dist.manifold Method
#' @export
dist.simplex = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  u = clr(p)
  v = clr(q)
  
  z = norm.simplex(u-v)
  return(z)
}


#' Riemannian exponential map for simplex data
#' 
#' @describIn RieExp.manifold Method
#' @export
RieExp.simplex = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  v = clr(p)
  z = inv.clr.simplex(v+u)
  return(z)
}


#' Riemannian logarithmic map for simplex data
#' 
#' @describIn RieExp.manifold Method
#' @export
RieLog.simplex = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  u = clr(p)
  v = clr(q)
  z = v-u
  return(z)
}


#' Basis on the tangent space at the point for simplex data
#' 
#' @describIn basis.manifold Method
#' @export
basis.simplex = function(p,dim=NULL){
  X = inv.clr.simplex(diag(length(p)))
  z = PCA.simplex(X)$vectors
  return(z)
}


#' Frechet mean for simplex data
#' 
#' @describIn FrechetMean.manifold Method
#' @export
FrechetMean.simplex = function(X){
  X = vec.to.mat(X)
  
  X.clr = clr(X)
  
  mu = colMeans(X.clr)
  mu = inv.clr.simplex(mu)[1,]
  return(mu)
}



# two below functions are useful for simplification
clr.to.real = function(x){
  y = apply(x,1,function(z){as.vector(outer(z,z,'-'))})
  return(t(y))
}

real.to.clr = function(y){
  m = sqrt(ncol(y))
  x = apply(y,1,function(z){rowMeans(matrix(z,m,m))})
  return(t(x))
}


#' Principal component analysis for simplex data
#' 
#' @describeIn PCA.manifold Method
#' @export
PCA.simplex = function(X){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  
  mu = FrechetMean.simplex(X)
  X.tangent = RieLog.simplex(mu,X)
  X.tangent = clr.to.real(X.tangent)
  
  Cov = t(X.tangent) %*% X.tangent / n
  C = eigen(Cov)
  
  values = C$values[1:(m-1)]/(2*m)
  vectors = sqrt(2*m)*real.to.clr(t(C$vectors[,1:(m-1)]))
  
  result = list(values=values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' Prediction score matrix for simplex data
#' 
#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.simplex = function(object,Xnew){
  object[['space']] = 'simplex'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


