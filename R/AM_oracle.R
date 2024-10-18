


#' @title Oracle additive regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Estimate additive component maps using an ADMM-based algorithm.
#' This function supports 'LASSO', 'SCAD', or 'MCP' penalty functions.
#' An \eqn{\ell^1}-type constrained bound \code{R} is multiplied by \eqn{\hat\sigma_Y = \sqrt{n^{-1}\sum_{i=1}^n d^2(Y_i,\hat\mu_Y)}}.
#' 
#' @inheritParams AM
#' 
#' @param proper.ind.mat a \eqn{s\times 2} index matrix such that \eqn{\mathcal{S}=\{(j,k) : m_{jk}^*\neq0\}}.
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
AM.oracle = function(Xorg,Yorg,Yspace,proper.ind.mat=NULL,degree=0,bandwidths.list=NULL,
                     transform='Gaussian',normalize=FALSE,ngrid=51,Kdenom_method='numeric',phi=1,max.iter=200,threshold=1e-6,SBF.comp=NULL){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.manifold(Yspace)
  
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
  
  # preprocessing for X
  Xdim.max = max(proper.ind.mat[,2])
  
  ## PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  Xdims = sapply(X,ncol)
  
  ## transformation
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  
  ## slice X_jk for (j,k) in proper.ind.mat
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  X = X[,all.indices]
  P = ncol(X)
  
  # define some functions for SBF from the SBF.preprocessing function
  # kde.1d: p list - (g,r,r) cube
  # proj: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2,r1,r2) cube
  # tildem: p list - (g,r,m) cube
  if (is.null(SBF.comp)){
    if (is.null(bandwidths.list)){
      bandwidths = get.rule.of.thumbs.bandwidths(X,degree)
    } else{
      bandwidths = sapply(1:P,function(j){ind=proper.ind.mat[j,];bandwidths.list[[ind[1]]][ind[2]]})
    }
    min.bandwidths = get.min.bandwidths(X,degree)
    bandwidths = pmax(bandwidths,min.bandwidths)
    
    SBF.comp = SBF.preprocessing(X,LogY,bandwidths,degree,ngrid,Kdenom_method)
  }
  bandwidths = SBF.comp[['bandwidths']]
  
  
  # apply AM_each function in cpp
  object = AM_each(SBF.comp,Ymu,Yspace,0,1e10,'LASSO',0,phi,0,max.iter,threshold)
  
  
  # compute other parameters
  object[['mhat.norm']] = cbind(1:P,proper.ind.mat,object[['mhat.norm']])
  colnames(object[['mhat.norm']]) = c('index','j','k','mhat.norm')
  
  X.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  X.vectors = reduce.dimension(X.vectors,Xdim.max,margin=2)
  proper.ind.mat.all = object[['mhat.norm']][which(object[['mhat.norm']][,'mhat.norm']!=0),,drop=FALSE]
  proper.ind.mat = proper.ind.mat.all[,2:3,drop=FALSE]
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
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
  object[['all.indices']] = all.indices
  object[['proper.ind.mat.all']] = proper.ind.mat.all
  object[['proper.ind.mat']] = proper.ind.mat
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  class(object) = 'AM'
  
  return(object)
}


