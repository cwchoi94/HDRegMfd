### ADMM|MM algorithm update for high-dimensional generalized linear regression
### developed only for LASSO, SCAD, MCP




#' @title High-dimensional Quantile Linear Regression for manifold-valued covariates.
#' 
#' @description 
#' Estimate Hilbert-Schmidt operators using an I-LAMM algorithm.
#' This function supports 'LASSO', 'SCAD', or 'MCP' penalty functions.
#' 
#' @inheritParams LM
#' 
#' @param Yorg an \eqn{n\times m} matrix of responses.
#' @param tau the quantile of Y given X (default: 0.5).
#' @param h the bandwidth for smoothing a check function.
#' @param kernel the kernel function for smoothing a check function.
#' @param phi0 the initial \eqn{phi} for updating \eqn{\hat{\beta}}.
#' @param c.phi the multiplying constant for \eqn{\phi}.
#'
#' @return a '\code{QM}' object with the following compnents:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{tau}{the quantile.}
#'       \item{kernel}{the kernel function for smoothing a check function.}
#'       \item{h}{the bandwidth for smoothing a check function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{iter.inner}{a vector of iteration numbers for each sub-iteration.}
#'       \item{runtime}{the running time (HH:MM:SS).}
#'       \item{runtime.second}{the running time (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
QM = function(Xorg,Yorg,tau=0.5,h=NULL,kernel='Gaussian',penalty='LASSO',gamma=0,lambda=0.1,Xdim.max=100,
              phi0=1e-4,c.phi=1.1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validity of inputs
  Check.penalty(penalty)
  Check.kernel.QM(kernel)
  
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
  
  # compute the default value of h
  if (is.null(h)){
    h = max(0.05,tau * (1-tau) * (log(sum(Xdims)) / n)^(1/4))
  }
  
  # apply QM_each function in cpp
  object = QM_each(X,Yorg,lambda,Xdim.max,tau,h,kernel,penalty,gamma,phi0,c.phi,max.iter,threshold)
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['pca']] = pca
  object[['tau']] = tau
  object[['kernel']] = kernel
  object[['h']] = h
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['iter.inner']] = object$iter.inner[1:object$iter]
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  class(object) = 'QM'
  
  return(object)
}



#' @title Prediction for quantile Linear Models
#' 
#' @description 
#' Predict \eqn{\hat{Q}_{Y|X}(X_{new})} for the given \eqn{X_{new}} using an \code{\link{QM}} object.
#' 
#' @param object an \code{\link{QM}} object.
#' @param Xnew a new list of covariates with the following components (see also \code{\link{covariates.generate}}):
#' \describe{
#'       \item{j}{a \eqn{p} list of manifold-valued covariates, where each \eqn{j}th element is an \eqn{n'\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}, see \code{\link{Check.manifold}}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
#'
#' @return an \eqn{n'\times m} matrix of predicted values \eqn{\hat{Q}_{Y|X}(X_{new})}.
#' @export
predict.QM = function(object,Xnew){
  Xnew = predict.PCA.manifold.list(object$pca,Xnew)
  Xnew = reduce.dimension(Xnew,object$Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  n2 = nrow(Xnew)
  theta = Xnew %*% object$beta + matrix(rep(object$beta0,n2),nrow=n2,byrow=TRUE)
  return(theta)
}



