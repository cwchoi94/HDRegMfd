



#' @title Minimal Bandwidths
#' 
#' @description
#' Computes the minimal bandwidths to ensure that one-dimensional KDEs are invertible on the interval \eqn{[0,1]}.
#' The function returns a vector of minimal bandwidths, each multiplied by an adjustment constant.
#' 
#' @param X an \eqn{n\times p} score matrix supported \eqn{[0,1]}.
#' @param degree the degree of the local polynomial fitting. Options are: 0 for Nadaraya-Watson (default), 1 for local-linear, etc.
#' @param const a constant to adjust the minimal bandwidths.
#' 
#' @return a \eqn{p} vector of adjusted minimal bandwidths.
#' @export
get.min.bandwidths = function(X,degree=0,const=1.001){
  n = nrow(X)
  p = ncol(X)
  r = degree+1
  
  h.min.vec = sapply(1:p,function(j){
    x = sort(X[,j])
    h.min = max(x[r],diff(x,lag=r)/2,1-x[n-r])
    return(h.min)
  })
  
  h.min.vec = h.min.vec * const
  
  return(h.min.vec)
}


#' @title Minimal Bandwidths for K-fold CV
#' 
#' @description
#' Computes the minimal bandwidths to ensure that one-dimensional KDEs are invertible on the interval \eqn{[0,1]}.
#' The function returns a vector of minimal bandwidths, each multiplied by an adjustment constant.
#' 
#' @param X.list a list of \eqn{n\times p} score matrices supported \eqn{[0,1]}.
#' @param degree the degree of the local polynomial fitting. Options are: 0 for Nadaraya-Watson (default), 1 for local-linear, etc.
#' @param const a constant to adjust the minimal bandwidths.
#' 
#' @return a \eqn{p} vector of adjusted minimal bandwidths.
#' @export
get.min.bandwidths.kfold = function(X.list,degree=0,const=1.001){
  
  kfold = length(X.list)
  p = ncol(X.list[[1]])
  r = degree+1
  h.min.mat = matrix(0,kfold,p)
  for (i in 1:kfold){
    X = X.list[[i]]
    n = nrow(X)
    
    h.min.vec = sapply(1:p,function(j){
      x = sort(X[,j])
      h.min = max(x[r],diff(x,lag=r)/2,1-x[n-r])
      return(h.min)
    })
    h.min.mat[i,] = h.min.vec
  }
  
  h.min = apply(h.min.mat,2,max) * const
  
  return(h.min)
}




#' @title Minimal Bandwidths
#' 
#' @description
#' Computes the minimal bandwidths to ensure that one-dimensional KDEs are invertible on the interval \eqn{[0,1]}.
#' The function returns a vector of minimal bandwidths, each multiplied by an adjustment constant.
#' 
#' @param X an \eqn{n\times p} score matrix supported \eqn{[0,1]}.
#' @param degree the degree of the local polynomial fitting. Options are: 0 for Nadaraya-Watson (default), 1 for local-linear, etc.
#' @param const a constant to adjust the minimal bandwidths.
#' 
#' @return a \eqn{p} vector of adjusted minimal bandwidths.
#' @export
get.rule.of.thumbs.bandwidths = function(X,degree=0,const=1.001){
  n = nrow(X)
  p = ncol(X)
  
  h.min = get.min.bandwidths(X,degree,const)
  h.vec = apply(X,2,std)
  if (degree==0){
    h.vec = h.vec * n**(-1/5)
  }else{
    h.vec = h.vec * n**(-1/(2*degree+3))
  }
  h.vec = pmax(h.vec,h.min)
  
  return(h.vec)
}






#' @title Kernel Density Estimator
#' 
#' @description
#' Computes one- and two-dimensional kernel density estimators for \eqn{X_j} on the interval \eqn{[0,1]} for \eqn{1\le j\le p}.
#' This function also computes the projection operator \eqn{\Pi_{jk}:L^2(f_k) \to L_2(f_j)} for \eqn{1\le j\neq k\le p}.
#' It uses the normalized Epanechnikov kernel for computing the KDEs.
#' 
#' @param X an \eqn{n\times p} score matrix. 
#' @param bandwidths a \eqn{p} vector of bandwidths.
#' @param degree the degree of the local polynomial fitting. Options are: 0 for Nadaraya-Watson (default), 1 for local-linear, etc.
#' @param ngrid the number of grid points for evaluating the KDE (default: 51).
#' @param Kdenom_method a method used for the denominator in the normalized kernel calculation.
#' @param is.proj logical. If \code{TRUE}, this function gives the projection operators.
#' 
#' @return a \code{KDE} object with the following components: 
#' \describe{
#'       \item{kde.1d}{a \eqn{p} list of one-dimensional KDE where each element is (ngrid,degree+1,degree+1) array.}
#'       \item{kde.2d}{a \eqn{(p*p*ngrid*ngrid,degree+1,degree+1)} array of two-dimensional KDE. This only exists when \code{is.proj=FALSE}.}
#'       \item{proj}{a \eqn{(p\times p\times ngrid\times ngrid) \times (degree+1) \times (degree+1)}-dimensional array. This only exists when \code{is.proj=TRUE}.}
#'       \item{kde.1d.inv}{a \eqn{p} list of inverse matrices of one-dimensional KDE where each element is (ngrid,degree+1,degree+1) array.}
#'       \item{grids}{grid points.}
#'       \item{weights}{a vector of weights for fast computing numerical integration.}
#' }
#' @export
KDE = function(X,bandwidths,degree=0,ngrid=51,Kdenom_method='numeric',is.proj=FALSE){
  
  if (ngrid<=1){ngrid = 2}
  p = ncol(X)
  grids = seq(0,1,length.out=ngrid)
  
  weights = rep(1/(ngrid-1),ngrid)
  weights[1] <- weights[ngrid] <- 1/(2*(ngrid-1))
  
  kde = KDE_(X,bandwidths,grids,weights,degree,Kdenom_method,is.proj)
  kde[['kde.1d']] = lapply(1:p,function(j){kde[['kde.1d']][(ngrid*(j-1)+1):(ngrid*j),,,drop=FALSE]})
  kde[['kde.1d.inv']] = lapply(1:p,function(j){kde[['kde.1d.inv']][(ngrid*(j-1)+1):(ngrid*j),,,drop=FALSE]})
  
  gc()
  class(kde) = 'KDE'
  
  return(kde)
}








#' @title SBF preprocessing
#' 
#' @description
#' Computes the pre-requisite components for smooth backfitting (SBF).
#' It uses \code{\link{KDE}} for density estimation.
#' 
#' @inheritParams KDE
#' 
#' @param LogY an \eqn{n\times m} response matrix. 
#' 
#' @return a \code{SBF.comp} object with the following components: 
#' \describe{
#'       \item{tildem}{a \eqn{p} list of one-dimensional KDE where each element is (ngrid,degree+1,m) array.}
#'       \item{kde.1d}{a \eqn{p} list of \eqn{\tilde{m}_j} where each element is (one-dimensional KDE where each element is (ngrid,degree+1,degree+1) array.}
#'       \item{proj}{a \eqn{(p\times p\times ngrid\times ngrid) \times (degree+1) \times (degree+1)}-dimensional array.}
#'       \item{bandwidths}{a \eqn{p} vector of bandwidths.}
#'       \item{grids}{grid points.}
#'       \item{weights}{a vector of weights for fast computing numerical integration.}
#' }
#' @export
SBF.preprocessing = function(X,LogY,bandwidths,degree=0,ngrid=51,Kdenom_method='numeric'){
  
  if (ngrid<=1){ngrid = 2}
  grids = seq(0,1,length.out=ngrid)
  
  weights = rep(1/(ngrid-1),ngrid)
  weights[1] <- weights[ngrid] <- 1/(2*(ngrid-1))
  
  SBF.comp = SBF_preprocessing(X,LogY,bandwidths,grids,weights,degree,Kdenom_method)
  class(SBF.comp) = 'SBF.comp'
  
  return(SBF.comp)
}











# for fast computation in simulation
compute.SBF.comp = function(Xorg,Yorg,Yspace,degree=0,Xdim.max.list=NULL,transform='Gaussian',normalize=FALSE,ngrid=51,Kdenom_method='numeric'){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.penalty(penalty)
  Check.manifold(Yspace)
  
  if ((penalty=='SCAD') & (gamma<2)){
    gamma = 3.7
  } else if ((penalty=='MCP') & (gamma<1)){
    gamma = 3
  }
  
  if (ngrid<=1){
    ngrid = 2
  }
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  
  # preprocessing for X
  ## PCA
  pca = PCA.manifold.list(Xorg)
  X_ = predict(pca,Xorg)
  
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(X_,ncol)))}
  Xdim.max.max = max(Xdim.max.list)
  X = reduce.dimension(X_,Xdim.max.max)
  Xdims = sapply(X,ncol)
  
  ## transformation
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  
  # define some functions for SBF from the SBF.preprocessing function
  # kde.1d: p list - (g,r,r) cube
  # proj: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2,r1,r2) cube
  # tildem: p list - (g,r,m) cube
  bandwidths = get.rule.of.thumbs.bandwidths(X,degree)
  
  SBF.comp = SBF.preprocessing(X,LogY,bandwidths,degree,ngrid,Kdenom_method)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  SBF.comp[['runtime']] = runtime
  
  return(SBF.comp) 
}




















