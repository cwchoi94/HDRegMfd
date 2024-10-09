


#' @title Generalized Cross-Validation for High-Dimensional Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements a generalized cross-validation (GCV) for high-dimensional Hilbert-Schmidt linear models.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function 'LM_CV' in 'LM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{LM}}.
#' 
#' @inheritParams LM.CV
#' 
#' @param Xorgnew a list of new covariates.
#' @param Yorgnew an \eqn{n'\times m} new responses.
#'
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,Yspace,lambda.list,Xdim.max.list,R.list,penalty='LASSO',gamma=0,
                  phi=1,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
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
  
  # PCA for X
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  Xnew = predict(pca,Xorgnew)
  
  # projection of Yorg and Yorgnew to the tangent space
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  LogYnew = RieLog.manifold(Ymu,Yorgnew,Yspace)
  
  # Use LM_GCV function to obtain the optimal parameters
  result = LM_GCV(X,LogY,Xnew,LogYnew,Ymu,Yspace,lambda.list,Xdim.max.list,R.list,penalty,gamma,phi,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = LM_each(X,LogY,Ymu,Yspace,opt.lambda,opt.Xdim.max,opt.R,penalty,gamma,phi,eta,max.iter,threshold)
  
  # compute other parameters
  Xdims = object$Xdims
  Xdims_cumul = c(0,cumsum(Xdims))
  
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['pca']] = pca
  object[['Ymu']] = Ymu
  object[['Yspace']] = Yspace
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  class(object) = 'LM'
  
  return(object)
}




