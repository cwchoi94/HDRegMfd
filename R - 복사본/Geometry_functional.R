



#' @describeIn inner.manifold Method
#' @export
inner.functional = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowMeans(u * v)
  return(z)
}

inner.each.functional = function(u,v,p=NULL){
  m = length(u)
  z = u %*% v / m
  return(z)
}


#' @describeIn norm.manifold Method
#' @export
norm.functional = function(u,p=NULL){
  z = sqrt(inner.functional(u,u,p))
  return(z)
}


#' @describeIn dist.manifold Method
#' @export
dist.functional = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = sqrt(rowMeans((p-q)^2))
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieExp.functional = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  z = p+u
  return(z)
}


#' @describeIn RieExp.manifold Method
#' @export
RieLog.functional = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = q-p
  return(z)
}


basis.functional_ = function(l,ngrid=100,t=NULL){
  if (is.null(t)){t = seq(0,1,length.out=ngrid)}
  if (l%%2==1){
    phi = 2^(1/2)*cos((l-1)*pi*t)
  } else{
    phi = 2^(1/2)*sin(l*pi*t)
  }
  return(phi)
}


#' @describeIn Basis.manifold Method
#' @export
basis.functional = function(p,dim=50){
  ngrid = length(p)
  z = sapply(1:dim,function(j){basis.functional_(j,ngrid)})
  z = t(z)/norm.functional(t(z)) # normalize
  return(z)
}


#' @describeIn FrechetMean.manifold Method
#' @export
FrechetMean.functional = function(X){
  X = vec.to.mat(X)
  mu = colMeans(X)
  return(mu)
}


#' @describeIn PCA.manifold Method
#' @export
PCA.functional = function(X,alpha=0.95){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  torg = seq(0,1,length.out=m)
  
  mu = FrechetMean.functional(X)
  
  X.centered = RieLog.functional(mu,X)
  Cov = (t(X.centered) %*% X.centered) / n
  eig.decomp = eigen(Cov)
  
  W = norm.functional(eig.decomp$vectors)
  if (length(W)==1){
    Wmat = W^(-1)
  } else{
    Wmat = diag(W^(-1))
  }
  
  values = eig.decomp$values * W^2
  vectors = t(eig.decomp$vectors) %*% Wmat
  
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
predict.PCA.functional = function(object,Xnew){
  object[['space']] = 'functional'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}

