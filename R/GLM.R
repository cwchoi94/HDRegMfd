### ADMM|MM algorithm update for high-dimensional generalized linear regression
### developed only for LASSO, SCAD, MCP





#' @title High-dimensional generalized Linear Regression for manifold-valued covariates.
#' 
#' @description 
#' Estimate Hilbert-Schmidt operators using an ADMM-based algorithm.
#' This function supports 'LASSO', 'SCAD', or 'MCP' penalty functions.
#' 
#' @inheritParams LM
#' 
#' @param Yorg an \eqn{n\times m} matrix of responses.
#' @param link a link function, see \code{\link{Check.link}}.
#'
#' @return a 'GLM' object with the following compnents:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{link}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM = function(Xorg,Yorg,lambda=0.1,Xdim.max=100,R=100,penalty='LASSO',link='binomial',
               phi=1,gamma=0,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validity of inputs
  Check.penalty(penalty)
  Check.link(link)
  
  if ((penalty=='SCAD') & (gamma<2)){
    gamma = 3.7
  } else if ((penalty=='MCP') & (gamma<1)){
    gamma = 3
  }
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  # apply GLM_each function in cpp
  object = GLM_each(X,Yorg,lambda,Xdim.max,R,penalty,link,phi,gamma,eta,max.iter,threshold)
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['pca']] = pca
  object[['link']] = link
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['runtime']] = runtime
  class(object) = 'GLM'
  
  return(object)
}



#' @title Prediction for generalized Linear Models
#' 
#' @description 
#' Predict \eqn{\hat{Y}_{new}} for the given \eqn{X_{new}} using an \code{\link{GLM}} object.
#' If 'is.inv.link' is FALSE, this function returns \eqn{\hat{\theta}_{new}}.
#' 
#' @param object an \code{\link{GLM}} object.
#' @param Xnew a new list of covariates with the following components (see also \code{\link{covariates.generate}}):
#' \describe{
#'       \item{j}{a \eqn{p} list of manifold-valued covariates, where each \eqn{j}th element is an \eqn{n'\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}, see \code{\link{Check.manifold}}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
#' @param is.inv.link whether to apply the inverse link function (default=TRUE).
#'
#' @return an \eqn{n'\times m} matrix of predicted values \eqn{\hat{Y}_{new}}.
#' @export
predict.GLM = function(object,Xnew,is.inv.link=TRUE){
  Xnew = predict.PCA.manifold.list(object$pca,Xnew)
  Xnew = reduce.dimension(Xnew,object$Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  n2 = nrow(Xnew)
  theta = Xnew %*% object$beta + matrix(rep(object$beta0,n2),nrow=n2,byrow=TRUE)
  
  if(is.inv.link){
    Ymu = Inv_Link(theta,object$link)
    return (Ymu)
  }else{
    return(theta)
  }
  
}



