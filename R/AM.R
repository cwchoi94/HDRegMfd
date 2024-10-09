



#' @title High-dimensional additive regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Estimate additive component maps using an ADMM-based algorithm.
#' This function supports 'LASSO', 'SCAD', or 'MCP' penalty functions.
#' An \eqn{\ell^1}-type constrained bound \code{R} is multiplied by \eqn{\hat\sigma_Y = \sqrt{n^{-1}\sum_{i=1}^n d^2(Y_i,\hat\mu_Y)}}.
#' 
#' @inheritParams LM
#' @inheritParams Transform.Score
#' @inheritParams KDE
#' 
#' @param R an \eqn{\ell^1}-type constrained bound, multiplied by \eqn{\hat\sigma_Y} (default: 100).
#' @param bandwidth.list a \eqn{p} list of bandwidths. If \code{bandwidth.list==NULL}, then we compute the bandwidths by \code{\link{get.bandwidths}}.
#' @param max.iter a maximum number of iterations (default: 200).
#' @param threshold a convergence threshold for the algorithm (default: 1e-6).
#' @param SBF.compe a \code{\link{SBF.comp}} object (default: \code{NULL}). If \code{SBF.comp==NULL}, this function computes the SBF.comp using \code{\link{SBF.comp}}.
#'
#' @return a \code{AM} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{Yspace}{the underlying space of of \eqn{Y}.}
#'       \item{kde.1d}{a \eqn{p} list of one-dimensional KDE where each element is (ngrid,degree+1,degree+1) array.}
#'       \item{X.vector}{a \eqn{p} list of the corresponding eigenvectors of \eqn{X_j}.}
#'       \item{transform}{a \code{\link{Transform.Score}} object.}
#'       \item{mhat}{a \eqn{p} list of estimated \eqn{\hat m_{jk}}, where each element is an \eqn{ngrid \times (degree+1) \times m} cube.}
#'       \item{mhat.norm}{a \eqn{p\times 4} matrix of norms of \eqn{\hat m_{jk}}, with columns c('index','j','k','mhat.norm').}
#'       \item{all.indices}{an index set used in the estimation.}
#'       \item{proper.ind.mat}{a subset of \code{mhat.norm} with \eqn{\mathcal{S}=\{(j,k): \hat{m_{jk}}\neq0\}}.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
AM = function(Xorg,Yorg,Yspace,degree=0,penalty='LASSO',gamma=0,lambda=0.1,Xdim.max=100,R=100,bandwidths.list=NULL,
              transform='Gaussian',normalize=FALSE,ngrid=51,Kdenom_method='numeric',phi=1,eta=1e-3,max.iter=200,threshold=1e-6,SBF.comp=NULL){
  
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
  m = ncol(LogY)
  
  Ysigma = sqrt(mean(dist.manifold(Yorg,Ymu,Yspace)**2))
  R = Ysigma * R
  
  # preprocessing for X
  ## PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  Xdims = sapply(X,ncol)
  
  ## transformation
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  P = ncol(X)
  
  # define some functions for SBF from the SBF.preprocessing function
  # kde.1d: p list - (g,r,r) cube
  # proj: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2,r1,r2) cube
  # tildem: p list - (g,r,m) cube
  if (is.null(SBF.comp)){
    if (is.null(bandwidths.list)){
      bandwidths = get.rule.of.thumbs.bandwidths(X,degree)
    } else{
      bandwidths = do.call(c,sapply(1:p,function(j){bandwidths.list[[j]][1:Xdims[j]]}))
    }
    min.bandwidths = get.min.bandwidths(X,degree)
    bandwidths = pmax(bandwidths,min.bandwidths)
    
    SBF.comp = SBF.preprocessing(X,LogY,bandwidths,degree,ngrid,Kdenom_method)
  }
  
  
  # apply AM_each function in cpp
  object = AM_each(SBF.comp,Ymu,Yspace,lambda,R,penalty,gamma,phi,eta,max.iter,threshold)
  
  
  # compute other parameters
  object[['mhat.norm']] = cbind(index.mat,object[['mhat.norm']])
  colnames(object[['mhat.norm']]) = c('index','j','k','mhat.norm')
  
  X.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  X.vectors = reduce.dimension(X.vectors,Xdim.max,margin=2)
  proper.ind.mat.all = object[['mhat.norm']][which(object[['mhat.norm']][,'mhat.norm']!=0),]
  proper.ind.mat = proper.ind.mat.all[,2:3]
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['P']] = P
  object[['m']] = m
  object[['pca']] = pca
  object[['Ymu']] = Ymu
  object[['Yspace']] = Yspace
  object[['kde.1d']] = SBF.comp[['kde.1d']]
  object[['X.vectors']] = X.vectors
  object[['Xdim.max']] = Xdim.max
  object[['transform']] = object.transform
  object[['bandwidths']] = bandwidths
  object[['all.indices']] = index.mat[,1]
  object[['proper.ind.mat.all']] = proper.ind.mat.all
  object[['proper.ind.mat']] = proper.ind.mat
  object[['runtime']] = runtime
  class(object) = 'AM'
  
  return(object)
}



#' @title Prediction for Additive Models.
#' 
#' @description 
#' Predict \eqn{\hat{Y}_{new}} for the given \eqn{X_{new}} using an \code{\link{AM}} object.
#' 
#' @inheritParams predict.LM
#' 
#' @param object an \code{\link{AM}} object.
#'
#' @return an \eqn{n'\times m} matrix of predicted values \eqn{\hat{Y}_{new}}.
#' @export
predict.AM = function(object,Xorgnew){
  P = object[['P']]
  m = object[['m']]
  grids = object[['grids']]
  all.indices = object[['all.indices']]
  mhat = object[['mhat']]
  
  Xnew_ = predict.PCA.manifold.list(object$pca,Xorgnew)
  Xnew = reduce.dimension(Xnew_,object$Xdim.max)
  Xnew = predict(object$transform,Xnew)
  Xnew = Xnew[,all.indices]
  
  # linear interpolation
  mhatnew = lapply(object[['proper.ind.mat.all']][,1],function(j){
    sapply(1:m,function(k){approx(grids,mhat[[j]][,1,k],Xnew[,j])$y})
  })
  
  LogYhat = Reduce('+',mhatnew)
  if (is.null(LogYhat)){
    LogYhat = matrix(0,nrow(Xnew),m)
  }
  Yhat = RieExp.manifold(object[['Ymu']],LogYhat,object[['Yspace']])
  return(Yhat)
}









