


#' @title Generalized Cross-Validation for High-Dimensional Additive Models
#' 
#' @description 
#' Implements generalized cross-validation (GCV) for high-dimensional additive models.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function \code{AM_CV} in 'AM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{AM}}.
#' 
#' @inheritParams AM.CV
#' @inheritParams LM.GCV
#' 
#' @param h.grid grid length of h_list for each \eqn{h}th bandwidth.
#' @param h.count a number of h_list for each \eqn{h}th bandwidth.
#' @param R.list a vector of \eqn{\ell^1}-type constrained bounds, multiplied by \eqn{\hat\sigma_Y} (default: c(1e10)).
#' 
#' @return a \code{AM} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{kde.1d}{a \eqn{p} list of one-dimensional KDE where each element is (ngrid,degree+1,degree+1) array.}
#'       \item{transform}{a \code{\link{Transform.Score}} object.}
#'       \item{mhat}{a \eqn{p} list of estimated \eqn{\hat m_{jk}}, where each element is an \eqn{ngrid \times (degree+1) \times m} cube.}
#'       \item{mhat.norm}{a \eqn{p} vector of norms of \eqn{\hat m_{jk}}.}
#'       \item{col.indices}{an index set used in the estimation.}
#'       \item{proper.indices}{an estimated index set \eqn{\mathcal{S}=\{(j,k): \hat{m_{jk}}\neq0\}}.}
#'       \item{proper.indices.1d}{an estimated index set with a version of sub-vector of \code{col.indices}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the CV running time (HH:MM:SS).}
#'       \item{runtime.second}{the CV running time (second).}
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
AM.CBS.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,Yspace,proper.ind.mat=NULL,degree=0,h.grid=0.01,h.count=10,penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,R.list=NULL,
                      max.cv.iter=20,cv.threshold=1e-6,transform='Gaussian',normalize=FALSE,ngrid=51,Kdenom_method='numeric',phi=1,eta=1e-3,max.iter=200,threshold=1e-6){
  
  start.time = Sys.time()
  
  # check validity of inputs
  Check.penalty(penalty)
  Check.manifold(Yspace)
  
  if ((penalty=='SCAD') & (gamma<2)){
    gamma = 3.7
  } else if ((penalty=='MCP') & (gamma<1)){
    gamma = 3
  }
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  if(is.null(lambda.list)){lambda.list=c(0)}
  if(is.null(R.list)){R.list=c(1e10)}
  
  if (ngrid<=1){ngrid = 2}
  grids = seq(0,1,length.out=ngrid)
  
  weights = rep(1/(ngrid-1),ngrid)
  weights[1] <- weights[ngrid] <- 1/(2*(ngrid-1))
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  LogYnew = RieLog.manifold(Ymu,Yorgnew,Yspace)
  m = ncol(LogY)
  
  Ysigma = sqrt(mean(dist.manifold(Yorg,Ymu,Yspace)**2))
  R.list = Ysigma * R.list
  
  # preprocessing for X
  ## PCA
  pca = PCA.manifold.list(Xorg)
  X_ = predict(pca,Xorg)
  Xnew_ = predict(pca,Xorgnew)
  
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(X_,ncol)))}
  Xdim.max.max = max(Xdim.max.list)
  X = reduce.dimension(X_,Xdim.max.max)
  Xnew = reduce.dimension(Xnew_,Xdim.max.max)
  Xdims = sapply(X,ncol)
  
  if (is.null(proper.ind.mat)){
    proper.ind.mat = do.call(rbind,lapply(1:p,function(j){cbind(rep(j,Xdims[j]),1:Xdims[j])}))
  }
  
  ## transformation
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  Xnew = predict(object.transform,Xnew)
  
  ## slice X_jk for (j,k) in proper.ind.mat
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  X = X[,all.indices]
  Xnew = Xnew[,all.indices]
  P = ncol(X)
  
  index.mat = cbind(1:P,proper.ind.mat)
  
  # compute a matrix of bandwidth candidates, with (i,j) elements being ith candidate of h_j 
  min.bandwidths = get.min.bandwidths(X,degree)
  bandwidths.mat = sapply(min.bandwidths,function(x){x+h.grid*(1:h.count)})
  
  
  # Use AM_CBS function to obtain the optimal parameters
  # not compute AIC or BIC loss
  cv.type = 'NONE' 
  result = AM_CBS_GCV(X,LogY,Xnew,LogYnew,Ymu,Yspace,bandwidths.mat,grids,weights,lambda.list,Xdim.max.list,R.list,index.mat,cv.type,
                      penalty,gamma,degree,Kdenom_method,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c(sapply(1:P,function(i){paste0('h',i)}),'lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply the AM function with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.bandwidths = result$opt.bandwidths
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  ## preprocessing
  opt.SBF.comp = SBF.preprocessing(X,LogY,opt.bandwidths,degree,ngrid,Kdenom_method)
  opt.SBF.comp = SBF_preprocessing_reduce_dim(opt.SBF.comp,opt.Xdim.max,index.mat)
  
  ## apply AM_each
  object = AM_each(opt.SBF.comp,Ymu,Yspace,opt.lambda,opt.R,penalty,gamma,phi,eta,max.iter,threshold)
  
  
  # compute other parameters
  X = reduce.dimension(X_,opt.Xdim.max)
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  P = ncol(X)
  
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  X = X[,all.indices]
  P = ncol(X)
  
  object[['mhat.norm']] = cbind(1:P,proper.ind.mat,object[['mhat.norm']])
  colnames(object[['mhat.norm']]) = c('index','j','k','mhat.norm')
  
  X.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  X.vectors = reduce.dimension(X.vectors,opt.Xdim.max,margin=2)
  proper.ind.mat.all = object[['mhat.norm']][which(object[['mhat.norm']][,'mhat.norm']!=0),,drop=FALSE]
  proper.ind.mat = proper.ind.mat.all[,2:3,drop=FALSE]
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime.opt.second = as.numeric(difftime(Sys.time(),opt.start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['P']] = P
  object[['m']] = m
  object[['pca']] = pca
  object[['Ymu']] = Ymu
  object[['Yspace']] = Yspace
  object[['kde.1d']] = opt.SBF.comp[['kde.1d']]
  object[['X.vectors']] = X.vectors
  object[['Xdim.max']] = opt.Xdim.max
  object[['transform']] = object.transform
  object[['bandwidths']] = opt.bandwidths
  object[['bandwidths.mat']] = bandwidths.mat
  object[['all.indices']] = all.indices
  object[['proper.ind.mat.all']] = proper.ind.mat.all
  object[['proper.ind.mat']] = proper.ind.mat
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  object[['runtime.opt.second']] = runtime.opt.second
  class(object) = 'AM'
  
  return(object)
}





