


#' @describeIn inner.manifold Method
#' @export
inner.sphere = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowSums(u * v)
  return(z)
}

inner.each.sphere = function(u,v,p=NULL){
  z = u %*% v
  return(z)
}


#' @describeIn norm.manifold Method
#' @export
norm.sphere = function(u,p=NULL){
  z = sqrt(inner.sphere(u,u,p))
  return(z)
}


#' @describeIn dist.manifold Method
#' @export
dist.sphere = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = rowSums(p*q)
  z[z > 1] <- 1
  z[z < -1] <- -1
  z = acos(z)
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieExp.sphere = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  z = sapply(1:nrow(p),function(i){Exp1(p[i,],u[i,])})
  z = t(z)
  return(z)
}

# Use a function in the library "manifold"
Exp1 <- function(mu,v,tol=1e-10) {
  vNorm <- as.numeric(sqrt(crossprod(v)))
  if (!is.na(vNorm) && vNorm <= tol) {
    mu
  } else {
    cos(vNorm) * mu + sin(vNorm) * v / vNorm
  }
}


#' @describeIn RieExp.manifold Method
#' @export
RieLog.sphere = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = Log2(t(q),t(p))
  z = t(z)
  return(z)
}

# Use a function in the library "manifold"
Log2 <- function(X, Mu, tol=1e-10) {
  dimAmbient <- nrow(X)
  N <- ncol(X)
  cprod <- colSums(X * Mu)
  U <- X - matrix(cprod, dimAmbient, N, byrow=TRUE) * Mu
  uNorm <- sqrt(colSums(U^2))
  res <- matrix(0, dimAmbient, N)
  ind <- uNorm > tol
  distS <- acos(ifelse(cprod > 1, 1, ifelse(cprod < -1, -1, cprod)))
  res[, ind] <- (U * matrix(distS / uNorm, dimAmbient, N, byrow=TRUE))[, ind, drop=FALSE]
  res
}


#' @describeIn basis.manifold Method
#' @export
basis.sphere = function(p,dim=NULL){
  z = manifold::basisTan(manifold::createM('Sphere'),p)
  z = t(z)
  return(z)
}


#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.sphere = function(X){
  X = vec.to.mat(X)
  
  mu = manifold::frechetMean(manifold::createM('Sphere'),t(X))
  mu = t(mu)[1,]
  return(mu)
}


#' @describeIn PCA.manifold Method
#' @export
PCA.sphere = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  
  mu = FrechetMean.sphere(X)
  X.tangent = RieLog.sphere(mu,X)
  
  Cov = t(X.tangent) %*% X.tangent / n
  C = eigen(Cov)
  
  values = C$values[1:(m-1)]
  vectors = t(C$vectors[,1:(m-1)])
  
  result = list(values=values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.sphere = function(object,Xnew){
  object[['space']] = 'sphere'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


