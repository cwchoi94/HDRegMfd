


#' @title Cross-Validation for High-Dimensional Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements cross-validation (CV) for high-dimensional Hilbert-Schmidt linear models based on 'AIC' or 'BIC'.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function 'LM_CV' in 'LM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{LM}}.
#' 
#' @inheritParams LM
#' 
#' @param lambda.list a vector of non-negative penalty constants.
#' @param Xdim.max.list a vector of the maximum dimension to which \eqn{X_j} will be reduced.
#' @param R.list a vector of \eqn{\ell^1}-type constrained bounds.
#' @param cv.type a CV method, which must be one of 'AIC', 'BIC', and 'ABIC' (default: 'AIC').
#' @param max.cv.iter a maximum number of CV iterations (default 20).
#' @param cv.threshold a convergence threshold for the CV (default 1e-10).
#'
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
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
LM.CV = function(Xorg,Yorg,Yspace,lambda.list,Xdim.max.list,R.list,cv.type='AIC',penalty='LASSO',phi=1,gamma=0,
                 max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validity of inputs
  Check.penalty(penalty)
  Check.cv.type(cv.type)
  
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
  
  # projection of Yorg and Yorgnew onto the tangent space
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  
  # Use LM_CV function to obtain the optimal parameters
  result = LM_CV(X,LogY,Ymu,Yspace,lambda.list,Xdim.max.list,R.list,cv.type,
                 penalty,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply an LM function with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = LM_each(X,LogY,Ymu,Yspace,opt.lambda,opt.Xdim.max,opt.R,penalty,phi,gamma,eta,max.iter,threshold)
  
  # compute other parameters
  Xdims = object$Xdims
  Xdims_cumul = c(0,cumsum(Xdims))
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
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




