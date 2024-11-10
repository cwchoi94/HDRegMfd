


#' @title Cross-Validation for High-Dimensional Additive Models
#' 
#' @description 
#' Implements cross-validation (CV) for high-dimensional additive models based on 'AIC' or 'BIC'.
#' The CV process follows the coordinate-wise variable selection and is implemented using a function \code{AM_CV} in 'AM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{AM}}.
#' 
#' @inheritParams AM
#' @inheritParams LM.CV
#' 
#' @param R.list a vector of \eqn{\ell^1}-type constrained bounds, multiplied by \eqn{\hat\sigma_Y} (default: c(100)).
#' @param max.cv.iter a maximum number of CV iterations (default: 20).
#' @param cv.threshold a convergence threshold for the CV (default: 1e-6).
#' @param loss.type the type of loss function. Options are 'average' or 'integral' (default).
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
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the CV running time (HH:MM:SS).}
#'       \item{runtime.second}{the CV running time (second).}
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
AM.CV = function(Xorg,Yorg,Yspace,degree=0,cv.type='AIC',penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,R.list=NULL,bandwidths.list=NULL,
                 max.cv.iter=20,cv.threshold=1e-6,transform='Gaussian',normalize=FALSE,ngrid=51,Kdenom_method='numeric',phi=1,eta=1e-3,max.iter=200,threshold=1e-6,loss.type='integral',SBF.comp=NULL){
  
  start.time = Sys.time()
  
  # check validity of inputs
  Check.penalty(penalty)
  Check.manifold(Yspace)
  Check.loss.type(loss.type)
  
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
  if(is.null(lambda.list)){lambda.list=c(0)}
  if(is.null(R.list)){R.list=c(100)}
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  m = ncol(LogY)
  
  Ysigma = sqrt(mean(dist.manifold(Yorg,Ymu,Yspace)**2))
  R.list = Ysigma * R.list
  
  # preprocessing for X
  ## PCA
  pca = PCA.manifold.list(Xorg)
  X_ = predict(pca,Xorg)
  
  if(is.null(Xdim.max.list)){Xdim.max.list = c(min(max(sapply(X_,ncol)),ceiling(n**(1/3))))}
  Xdim.max.max = max(Xdim.max.list)
  X = reduce.dimension(X_,Xdim.max.max)
  Xdims = sapply(X,ncol)
  
  ## transformation
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  P = ncol(X)
  
  # define some functions for SBF
  # kde.1d: p list - (g,r,r) cube
  # projection: (p1,p2,g1,g2,r1,r2) = (p1*p2*g1*g2,r1,r2) cube
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
  bandwidths = SBF.comp[['bandwidths']]
  
  
  # Use AM_CV function to obtain the optimal parameters
  if (loss.type=='average'){
    result = AM_CV_average(SBF.comp,X,LogY,Ymu,Yspace,lambda.list,Xdim.max.list,R.list,index.mat,
                           cv.type,penalty,gamma,max.cv.iter,cv.threshold)
  }else{
    result = AM_CV_integral(SBF.comp,X,LogY,Ymu,Yspace,lambda.list,Xdim.max.list,R.list,index.mat,
                            cv.type,penalty,gamma,max.cv.iter,cv.threshold)
  }
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply the AM function with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  ## preprocessing
  opt.SBF.comp = SBF_preprocessing_reduce_dim(SBF.comp,opt.Xdim.max,index.mat)
  
  ## apply AM_each
  object = AM_each(opt.SBF.comp,Ymu,Yspace,opt.lambda,opt.R,penalty,gamma,phi,eta,max.iter,threshold)
  
  
  # compute other parameters
  X = reduce.dimension(X_,opt.Xdim.max)
  object.transform = Transform.Score(X,transform,normalize)
  index.mat = object.transform$index.mat
  X = predict(object.transform,X)
  P = ncol(X)
  
  object[['mhat.norm']] = cbind(index.mat,object[['mhat.norm']])
  colnames(object[['mhat.norm']]) = c('index','j','k','mhat.norm')
  
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
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
  object[['Xdim.max']] = opt.Xdim.max
  object[['transform']] = object.transform
  object[['bandwidths']] = bandwidths
  object[['all.indices']] = index.mat[,1]
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




