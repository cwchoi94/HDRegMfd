


#' Inner product on tangent space at p for Wasserstein data
#' 
#' @describIn inner.manifold Method
#' @export
inner.Wasserstein = function(u,v,p=NULL){
  u = vec.to.mat(u)
  v = vec.to.mat(v)
  
  u = vec.duplicate(u,nrow(v))
  v = vec.duplicate(v,nrow(u))
  
  z = rowMeans(u * v)
  return(z)
}

inner.each.Wasserstein = function(u,v,p=NULL){
  m = length(u)
  z = u %*% v / m
  return(z)
}


#' norm on tangent space at p for Wasserstein data
#' 
#' @describIn norm.manifold Method
#' @export
norm.Wasserstein = function(u,p=NULL){
  z = sqrt(inner.Wasserstein(u,u,p))
  return(z)
}


#' Geodesic distance for Wasserstein data
#' 
#' @describIn dist.manifold Method
#' @export
dist.Wasserstein = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = sqrt(rowMeans((p-q)^2))
  return(z)
}



#' Riemannian exponential map for Wasserstein data
#' 
#' @describIn RieExp.manifold Method
#' @export
RieExp.Wasserstein = function(p,u){
  p = vec.to.mat(p)
  u = vec.to.mat(u)
  
  p = vec.duplicate(p,nrow(u))
  u = vec.duplicate(u,nrow(p))
  
  z = p+u
  return(z)
}


#' Riemannian logarithmic map for Wasserstein data
#' 
#' @describIn RieExp.manifold Method
#' @export
RieLog.Wasserstein = function(p,q){
  p = vec.to.mat(p)
  q = vec.to.mat(q)
  
  p = vec.duplicate(p,nrow(q))
  q = vec.duplicate(q,nrow(p))
  
  z = q-p
  return(z)
}


basis.Wasserstein_ = function(l,ngrid=100,t=NULL){
  if (is.null(t)){t = seq(0,1,length.out=ngrid)}
  phi = t^l
  return(phi)
}

#' Basis on the tangent space at the point for Wasserstein data
#' 
#' @describIn Basis.manifold Method
#' @export
basis.Wasserstein = function(p,dim=50){
  ngrid = length(p)
  z = sapply(1:dim,function(j){basis.functional_(j,ngrid)})
  z = t(z)/apply(z,2,norm.Wasserstein) # normalize
  return(z)
}


#' Frechet mean for Wasserstein data
#' 
#' @describIn FrechetMean.manifold Method
#' @export
FrechetMean.Wasserstein = function(X){
  X = vec.to.mat(X)
  mu = colMeans(X)
  return(mu)
}


#' Principal component analysis for Wasserstein data
#' 
#' @describeIn PCA.manifold Method
#' @export
PCA.Wasserstein = function(X){
  X = vec.to.mat(X)
  n = nrow(X)
  m = ncol(X)
  torg = seq(0,1,length.out=m)
  
  L = fdapace::MakeFPCAInputs(IDs=rep(1:n,each=m),tVec=rep(torg,n),t(X))
  fpca = fdapace::FPCA(L$Ly,L$Lt)
  
  W = apply(fpca$phi,2,norm.Wasserstein)
  if (length(W)==1){
    Wmat = W^(-1)
  } else{
    Wmat = diag(W^(-1))
  }
  
  fpca$lambda = fpca$lambda * W^2
  fpca$phi = fpca$phi %*% Wmat
  
  values = fpca$lambda
  vectors = t(fpca$phi)
  mu = fpca$mu
  
  result = list(fpca=fpca,values=values,vectors=vectors,mu=mu,dim=nrow(vectors))
  return(result)
}


#' Prediction score matrix for Wasserstein data
#' 
#' @describeIn predict.PCA.manifold Method
#' @export
predict.PCA.Wasserstein = function(object,Xnew){
  object[['space']] = 'Wasserstein'
  scores = predict.PCA.manifold(object,Xnew)
  return (scores)
}

# predict.PCA.Wasserstein = function(object,Xnew){
#   n2 = nrow(Xnew)
#   m = ncol(Xnew)
#   tnew = seq(0,1,length.out=m)
#   
#   Lnew = fdapace::MakeFPCAInputs(IDs=rep(1:n2,each=m),tVec=rep(tnew,n2),t(Xnew))
#   scores = predict(object$fpca,Lnew$Ly,Lnew$Lt)$scores
#   return (scores)
# }

