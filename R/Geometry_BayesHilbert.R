
#' @title Inverse centered log-transform for density data
#' 
#' @description 
#' Inverse transform R^m to density data
#' inv.clr(x)_i = exp(x_i)/(sum_j exp(x_i)) * normalizer
#' 
#' @param x clr data, (n,m) matrix
#' @param normalizer normalized constant, default=1 
#' 
#' @return density data, (n,m) matrix
#' @export
inv.clr.BayesHilbert = function(x,normalizer=1){
  x = vec.to.mat(x)
  x = x - rowMeans(x)
  y = exp(x)
  y = y/rowMeans(y)*normalizer
  return(y)
}


#' Inner product on tangent space at p for Bayes-Hilbert data

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


#' norm on tangent space at p for Bayes-Hilbert data

#' @describeIn norm.manifold Method
#' @export
norm.BayesHilbert = function(u,p=NULL){
  z = sqrt(inner.BayesHilbert(u,u,p))
  return(z)
}


#' Geodesic distance for Bayes-Hilbert data

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


#' Riemannian exponential map for Bayes-Hilbert data

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


#' Riemannian logarithmic map for Bayes-Hilbert data

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


#' Basis on the tangent space at the point for Bayes-Hilbert data

#' @describeIn Basis.manifold Method
#' @export
basis.BayesHilbert = function(p,dim=50){
  ngrid = length(p)
  z = sapply(1:dim,function(j){basis.functional_(2*j,ngrid)})
  z = t(z)/apply(z,2,norm.BayesHilbert) # normalize
  return(z)
}


#' Frechet mean for Bayes-Hilbert data

#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.BayesHilbert = function(X){
  X = vec.to.mat(X)
  
  X = clr(X)
  mu = colMeans(X)
  mu = inv.clr.BayesHilbert(mu)[1,]
  return(mu)
}


#' Principal component analysis for Bayes-Hilbert data

#' @describeIn PCA.manifold Method
#' @export
PCA.BayesHilbert = function(X){
  X = vec.to.mat(X)
  X = clr(X)
  n = nrow(X)
  m = ncol(X)
  torg = seq(0,1,length.out=m)
  
  L = fdapace::MakeFPCAInputs(IDs=rep(1:n,each=m),tVec=rep(torg,n),t(X))
  fpca = fdapace::FPCA(L$Ly,L$Lt)
  
  W = apply(fpca$phi,2,function(x){norm.BayesHilbert(x)})
  
  fpca$lambda = fpca$lambda * W^2
  fpca$phi = fpca$phi %*% diag(W^(-1))
  
  values = fpca$lambda
  vectors = t(fpca$phi)
  mu = inv.clr.BayesHilbert(vec.to.mat(fpca$mu))[1,]
  
  result = list(fpca=fpca,values=values,vectors=vectors,mu=mu,dim=nrow(vectors))
  return(result)
}


#' Prediction score matrix for Bayes-Hilbert data

#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.BayesHilbert = function(object,Xnew){
  object[['space']] = 'BayesHilbert'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


# predict.PCA.BayesHilbert = function(object,Xnew){
#   Xnew = clr(Xnew)
#   n2 = nrow(Xnew)
#   m = ncol(Xnew)
#   tnew = seq(0,1,length.out=m)
# 
#   Lnew = fdapace::MakeFPCAInputs(IDs=rep(1:n2,each=m),tVec=rep(tnew,n2),t(Xnew))
#   scores = predict(object$fpca,Lnew$Ly,Lnew$Lt)$scores
#   return (scores)
# }

