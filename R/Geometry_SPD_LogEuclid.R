

# matrix logarithmic and exponential maps for SPD matrices
logMvec = function(M){
  M = vec.to.mat(M)
  
  n = nrow(M)
  m = ncol(M)
  dim = sqrt(m)
  
  logM = sapply(1:n,function(i){x=matrix(M[i,],dim,dim);as.vector(manifold::LogMSPD(x))})
  logM = t(logM)
  
  return(logM)
}

expMvec = function(M){
  M = vec.to.mat(M)
  
  n = nrow(M)
  m = ncol(M)
  dim = sqrt(m)
  
  expM = sapply(1:n,function(i){x=matrix(M[i,],dim,dim);as.vector(manifold::ExpM(x))})
  expM = t(expM)
  
  return(expM)
}



#' @describeIn inner.manifold Method
#' @export
inner.SPD.LogEuclid = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowSums(u * v)
  return(z)
}

inner.each.SPD.LogEuclid = function(u,v,p=NULL){
  z = u %*% v
  return(z)
}


#' @describeIn norm.manifold Method
#' @export
norm.SPD.LogEuclid = function(u,p=NULL){
  z = sqrt(inner.SPD.LogEuclid(u,u,p))
  return(z)
}


#' @describeIn dist.manifold Method
#' @export
dist.SPD.LogEuclid = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  logP = logMvec(p)
  logQ = logMvec(q)
  
  z = sqrt(rowSums((logP-logQ)^2))
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieExp.SPD.LogEuclid = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  logP = logMvec(p)
  z = expMvec(logP+u)
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieLog.SPD.LogEuclid = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  logP = logMvec(p)
  logQ = logMvec(q)
  
  z = logQ-logP
  return(z)
}


#' @describeIn Basis.manifold Method
#' @export
basis.SPD.LogEuclid = function(p,dim=50){
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


#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.SPD.LogEuclid = function(X){
  X = vec.to.mat(X)
  m = ncol(X)
  dim = sqrt(m)
  
  logX = logMvec(X)
  logmu = colMeans(logX)
  mu = manifold::ExpM(matrix(logmu,dim,dim))
  mu = as.vector(mu)
  return(mu)
}


#' @describeIn PCA.manifold Method
#' @export
PCA.SPD.LogEuclid = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  dim = sqrt(m)
  k = dim*(dim+1)/2
  
  mu = FrechetMean.SPD.LogEuclid(X)
  X.tangent = RieLog.SPD.LogEuclid(mu,X)
  
  Cov = t(X.tangent) %*% X.tangent / n
  C = eigen(Cov)
  
  values = C$values[1:k]
  vectors = t(C$vectors)[1:k,]
  
  result = list(values=values,vectors=vectors,mu=mu,Cov=Cov,dim=nrow(vectors))
  return(result)
}


#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.SPD.LogEuclid = function(object,Xnew){
  object[['space']] = 'SPD.LogEuclid'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}


