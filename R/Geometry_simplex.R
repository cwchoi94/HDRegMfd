
#' @title Centered log-ratio transform 
#' 
#' @description 
#' Compute the centered log-ratio transform for compositional or density data.
#' Specifically, for each \eqn{x=(x_i)\in\mathbb{R}^m}, the transformation is give by
#' \deqn{\text{clr}(x)_i = \log(x_i) - \frac{1}{m}\sum_{j}\log(x_j)}.
#' 
#' @param x an \eqn{n\times m} matrix of simplex- or Bayes-Hilbert space-valued data.
#' 
#' @return an \eqn{n\times m} matrix of clr-transformed data.
#' @export
clr = function(x){
  x = vec.to.mat(x)
  y = log(x)
  y = y - rowMeans(y)
  return (y)
}


#' @title Inverse centered log-transform for simplex-valued data
#' 
#' @description 
#' Compute the inverse centered log-ratio transform to convert data from \eqn{\mathbb{R}^m} to compositional data, see \code{\link{clr}}.
#' Specifically, for each \eqn{x=(x_i)\in\mathbb{R}^m}, the transformation is give by
#' \deqn{\text{inv.clr}(x)_i = \frac{\exp(x_i)}{\sum_{j} \exp(x_j)} \times \text{normalizer}}.
#' 
#' @param x an \eqn{n\times m} matrix of the clr-transformed data.
#' @param normalizer a normalization constant, with a default value of 1.
#' 
#' @return a \eqn{n\times m} matrix of density data.
#' @export
inv.clr.simplex = function(x,normalizer=1){
  x = vec.to.mat(x)
  x = x - rowMeans(x)
  y = exp(x)
  y = y/rowSums(y)*normalizer
  return(y)
}


#' @describeIn inner.manifold Method
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

inner.each.simplex = function(u,v,p=NULL){
  m = length(u)
  a = outer(u,u,'-')
  b = outer(v,v,'-')
  z = sum(a*b)/(2*m)
  return(z)
}


#' @describeIn norm.manifold Method
#' @export
norm.simplex = function(u,p=NULL){
  z = sqrt(inner.simplex(u,u,p))
  return(z)
}


#' @describeIn dist.manifold Method
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


#' @describeIn RieExp.manifold Method
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


#' @describeIn RieLog.manifold Method
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


#' @describeIn basis.manifold Method
#' @export
basis.simplex = function(p,dim=NULL){
  X = inv.clr.simplex(diag(length(p)))
  z = PCA.simplex(X)$vectors
  return(z)
}


#' @describeIn FrechetMean.manifold Method
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


#' @describeIn PCA.manifold Method
#' @export
PCA.simplex = function(X,alpha=0.95){
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


#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.simplex = function(object,Xnew){
  object[['space']] = 'simplex'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


