






#' @title Kfold Cross-Validation for High-Dimensional Additive Models
#' 
#' @description 
#' Implements Kfold cross-validation (Kfold-CV) for high-dimensional additive models.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function \code{AM_CV} in 'AM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{AM}}.
#' 
#' @inheritParams AM.CBS.GCV
#' @inheritParams LM.kfold
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
AM.CBS.kfold = function(Xorg,Yorg,Yspace,proper.ind.mat=NULL,degree=0,h.grid=0.01,h.count=10,kfold=5,seed=NULL,penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,R.list=NULL,
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
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  if(is.null(lambda.list)){lambda.list=c(0)}
  if(is.null(R.list)){R.list=c(1e10)}
  
  if (ngrid<=1){ngrid = 2}
  grids = seq(0,1,length.out=ngrid)
  
  weights = rep(1/(ngrid-1),ngrid)
  weights[1] <- weights[ngrid] <- 1/(2*(ngrid-1))
  
  
  ## compute LogYall and R.list
  Ymu = FrechetMean.manifold(Yall,Yspace)
  LogYall = RieLog.manifold(Ymu,Yall,Yspace)
  m = ncol(LogYall)
  
  Ysigma = sqrt(mean(dist.manifold(Yall,Ymu,Yspace)**2))
  R.list = Ysigma * R.list
  
  ## compute Xdim.max.max
  pca = PCA.manifold.list(Xall)
  X_ = predict(pca,Xall)
  if(is.null(Xdim.max.list)){Xdim.max.list = c(min(max(sapply(X_,ncol)),ceiling(n**(1/3))))}
  Xdim.max.max = max(Xdim.max.list)
  
  # data split and preprocessing
  test.indices.list = get.kfold.test.indices(n,kfold,seed)
  
  Xorg.list = list()
  Xnew.list = list()
  LogY.list = list()
  LogYnew.list = list()
  Ymu.list = list()
  Xdims.mat = matrix(0,kfold,p)
  index.mat.list = list()
  for (i in 1:kfold){
    test.indices = test.indices.list[[i]]
    data.split = split.data.org(Xall,Yall,test.indices)
    
    # compute Xorg, Xnew, LogY, LogYnew, and Ymu for each kfold data.
    Ymu = FrechetMean.manifold(data.split$Yorg,Yspace)
    
    pca = PCA.manifold.list(data.split$Xorg)
    X_ = predict(pca,data.split$Xorg)
    Xnew_ = predict(pca,data.split$Xnew)
    X = reduce.dimension(X_,Xdim.max.max)
    Xnew = reduce.dimension(Xnew_,Xdim.max.max)
    
    Xdims = sapply(X,ncol)
    
    object.transform = Transform.Score(X,transform,normalize)
    index.mat = object.transform$index.mat
    X = predict(object.transform,X)
    Xnew = predict(object.transform,Xnew)
    
    Xorg.list[[i]] = X
    Xnew.list[[i]] = Xnew
    LogY.list[[i]] = RieLog.manifold(Ymu,data.split$Yorg,Yspace)
    LogYnew.list[[i]] = RieLog.manifold(Ymu,data.split$Ynew,Yspace)
    Ymu.list[[i]] = Ymu
    Xdims.mat[i,] = Xdims
    index.mat.list[[i]] = index.mat
  }
  
  # compute index.mat and Xdims
  ## index.mat: intersection of all index.mat
  ## Xdims: minimum dimension of all X_j
  tmp = lapply(index.mat.list,function(y){which(apply(y[,-1], 1, function(y_row) {
    any(apply(index.mat.list[[1]][,-1], 1, function(x_row){all(x_row == y_row)}))
  }))})
  index.mat = index.mat.list[[1]][Reduce(intersect,tmp),]
  Xdims = apply(Xdims.mat,2,min)
  
  if (is.null(proper.ind.mat)){
    proper.ind.mat = do.call(rbind,lapply(1:p,function(j){cbind(rep(j,Xdims[j]),1:Xdims[j])}))
  }
  
  ## slice X_jk for (j,k) in proper.ind.mat
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  for (i in 1:kfold){
    Xorg.list[[i]] = Xorg.list[[i]][,all.indices]
    Xnew.list[[i]] = Xnew.list[[i]][,all.indices]
  }
  P = length(all.indices)
  
  index.mat = cbind(1:P,proper.ind.mat)
  
  # compute a matrix of bandwidth candidates, with (i,j) elements being ith candidate of h_j 
  min.bandwidths = get.min.bandwidths.kfold(Xorg.list,degree)
  bandwidths.mat = sapply(min.bandwidths,function(x){x+h.grid*(1:h.count)})
  
  
  # Use AM_CBS function to obtain the optimal parameters
  # not compute AIC or BIC loss
  cv.type = 'NONE' 
  result = AM_CBS_kfold(Xorg.list,LogY.list,Xnew.list,LogYnew.list,Ymu.list,Yspace,kfold,bandwidths.mat,grids,weights,lambda.list,Xdim.max.list,R.list,index.mat,cv.type,
                        penalty,gamma,degree,Kdenom_method,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c(sapply(1:P,function(i){paste0('h',i)}),'lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply the AM function with the optimal parameters
  opt.bandwidths = result$opt.bandwidths
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  
  # compute LogY and X with the optimal parameters
  opt.start.time = Sys.time()
  
  # preprocessing for X
  ## PCA
  pca = PCA.manifold.list(Xall)
  X_ = predict(pca,Xall)
  X = reduce.dimension(X_,opt.Xdim.max)
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
  
  ## preprocessing
  opt.SBF.comp = SBF.preprocessing(X,LogYall,opt.bandwidths,degree,ngrid,Kdenom_method)
  
  ## apply AM_each
  object = AM_each(opt.SBF.comp,Ymu,Yspace,opt.lambda,opt.R,penalty,gamma,phi,eta,max.iter,threshold)
  
  
  # compute other parameters
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





