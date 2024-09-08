
#' @title Inverse centered log-transform for density data
#' 
#' @description 
#' Compute the inverse centered log-ratio transform to convert data from \eqn{\mathbb{R}^m} to density data, see \code{\link{clr}}.
#' Specifically, for each \eqn{x=(x_i)\in\mathbb{R}^m}, the transformation is give by
#' \deqn{\text{inv.clr}(x)_i = \frac{\exp(x_i)}{\frac{1}{m}\sum_{j} \exp(x_j)} \times \text{normalizer}}.
#' 
#' @param x an \eqn{n\times m} matrix of the clr-transformed data.
#' @param normalizer a normalization constant, with a default value of 1.
#' 
#' @return a \eqn{n\times m} matrix of density data.
#' @export
inv.clr.BayesHilbert = function(x,normalizer=1){
  x = vec.to.mat(x)
  x = x - rowMeans(x)
  y = exp(x)
  y = y/rowMeans(y)*normalizer
  return(y)
}


#' @describeIn inner.manifold Method
#' @export
inner.BayesHilbert = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowMeans(u * v)
  return(z)
}

inner.each.BayesHilbert = function(u,v,p=NULL){
  m = length(u)
  z = u %*% v / m
  return(z)
}


#' @describeIn norm.manifold Method
#' @export
norm.BayesHilbert = function(u,p=NULL){
  z = sqrt(inner.BayesHilbert(u,u,p))
  return(z)
}


#' @describeIn dist.manifold Method
#' @export
dist.BayesHilbert = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = clr(p)
  q = clr(q)
  
  z = sqrt(rowMeans((p-q)^2))
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieExp.BayesHilbert = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  p = clr(p)
  z = inv.clr.BayesHilbert(p+u)
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieLog.BayesHilbert = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  p = clr(p)
  q = clr(q)
  
  z = q-p
  return(z)
}


#' @describeIn Basis.manifold Method
#' @export
basis.BayesHilbert = function(p,dim=50){
  ngrid = length(p)
  z = sapply(1:dim,function(j){basis.functional_(2*j,ngrid)})
  z = t(z)/apply(z,2,norm.BayesHilbert) # normalize
  return(z)
}


#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.BayesHilbert = function(X){
  X = vec.to.mat(X)
  
  X = clr(X)
  mu = colMeans(X)
  mu = inv.clr.BayesHilbert(mu)[1,]
  return(mu)
}


#' @describeIn PCA.manifold Method
#' @export
PCA.BayesHilbert = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  torg = seq(0,1,length.out=m)
  
  mu = FrechetMean.BayesHilbert(X)
  
  X.centered = RieLog.BayesHilbert(mu,X)
  Cov = (t(X.centered) %*% X.centered) / n
  eig.decomp = eigen(Cov)
  
  W = norm.BayesHilbert(eig.decomp$vectors)
  if (length(W)==1){
    Wmat = W^(-1)
  } else{
    Wmat = diag(W^(-1))
  }
  
  values = eig.decomp$values * W^2
  vectors = t(eig.decomp$vectors) %*% Wmat
  mu = inv.clr.BayesHilbert(vec.to.mat(mu))[1,]
  
  if (alpha<0){alpha=0} else if (alpha>1){alpha=1}
  
  values.cumsum = cumsum(values)
  idx = min(which(values.cumsum/values.cumsum[length(values.cumsum)]>=alpha))
  idx = max(idx,2)
  
  values = values[1:idx]
  vectors = vectors[1:idx,]
  
  result = list(values=values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.BayesHilbert = function(object,Xnew){
  object[['space']] = 'BayesHilbert'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


