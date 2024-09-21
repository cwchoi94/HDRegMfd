### ADMM|MM algorithm update for high-dimensional linear regression
### developed only for LASSO, SCAD, MCP




#' @title Demension reduction for a list of covariates
#' 
#' @description 
#' Reduces the dimension of each \eqn{X_j} to \eqn{\min\{D_j,\text{Xdim.max}\}}.
#' 
#' @param X a \eqn{p} list of covariates, where each \eqn{X_j} is an \eqn{n\times D_j} matrix.
#' @param Xdim.max the maximum dimension to which \eqn{X_j} will be reduced, with a default value of 100.
#' @param margin the direction of the dimension reduction, with a default value of 1 (1: column-wise, 2: row-wise).
#' 
#' @return a \eqn{p} list of reduced covariates \eqn{X_j'}.
#' @export
reduce.dimension = function(X,Xdim.max=100,margin=1){
  if (!(margin %in% c(1,2))){stop("margin should be 1 or 2.")}
  
  if (margin==1){
    Xdims = sapply(X,ncol)
    p = length(Xdims)
    X = lapply(1:p,function(j){matrix(X[[j]][,1:min(Xdims[[j]],Xdim.max)],nrow=nrow(X[[j]]))})
  } else{
    Xdims = sapply(X,nrow)
    p = length(Xdims)
    X = lapply(1:p,function(j){matrix(X[[j]][1:min(Xdims[[j]],Xdim.max),],ncol=ncol(X[[j]]))})
  }
  
  return(X)
}


#' @title High-dimensional Hilbert-Schmidt linear regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Estimate Hilbert-Schmidt operators using an ADMM-based algorithm.
#' This function supports 'LASSO', 'SCAD', or 'MCP' penalty functions.
#' 
#' @param Xorg a list of covariates with the following components (see also \code{\link{covariates.generate}}):
#' \describe{
#'       \item{j}{a \eqn{p} list of manifold-valued covariates, where each \eqn{j}th element is an \eqn{n\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}, see \code{\link{Check.manifold}}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
#' @param Yorg an \eqn{n\times m} matrix of manifold-valued responses.
#' @param Yspace the name of the underlying space \eqn{\mathcal{M}_Y} of \eqn{Y}.
#' @param lambda a non-negative penalty constant (default: 0.1).
#' @param Xdim.max the maximum dimension to which \eqn{X_j} will be reduced (default: 100).
#' @param R an \eqn{\ell^1}-type constrained bound (default: 100).
#' @param penalty the name of a penalty function. This must be one of 'LASSO', 'SCAD', or 'MCP' (default: 'LASSO').
#' @param phi a parameter for computing the ADMM-based algorithm for the majorized objective function (default: 1).
#' @param gamma a parameter for SCAD (default: 3.7) or MCP (default: 3).
#' @param eta a parameter for computing the ADMM-based algorithm for the proximal norm square (default: 1e-3).
#' @param max.iter a maximum number of iterations (default: 500).
#' @param threshold a convergence threshold for the algorithm (default: 1e-10).
#'
#' @return an \code{LM} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM = function(Xorg,Yorg,Yspace,lambda=0.1,Xdim.max=100,R=100,penalty='LASSO',
              phi=1,gamma=0,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  inner = eval(parse(text=paste0('inner.each.',Yspace)))
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  
  # apply LM_each function in cpp
  object = LM_each(X,LogY,Ymu,Yspace,lambda,Xdim.max,R,penalty,phi,gamma,eta,max.iter,threshold)
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,Xdim.max,margin=2)
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
  object[['runtime']] = runtime
  class(object) = 'LM'
  
  return(object)
}



#' @title Prediction for Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Predict \eqn{\hat{Y}_{new}} for the given \eqn{X_{new}} using an \code{\link{LM}} object.
#' 
#' @param object an \code{\link{LM}} object.
#' @param Xnew a new list of covariates with the following components (see also \code{\link{covariates.generate}}):
#' \describe{
#'       \item{j}{a \eqn{p} list of manifold-valued covariates, where each \eqn{j}th element is an \eqn{n'\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}, see \code{\link{Check.manifold}}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
#'
#' @return an \eqn{n'\times m} matrix of predicted values \eqn{\hat{Y}_{new}}.
#' @export
predict.LM = function(object,Xnew){
  Xnew = predict.PCA.manifold.list(object$pca,Xnew)
  Xnew = reduce.dimension(Xnew,object$Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  LogYhat = Xnew %*% object$beta
  Yhat = RieExp.manifold(object[['Ymu']],LogYhat,object[['Yspace']])
  return(Yhat)
}



