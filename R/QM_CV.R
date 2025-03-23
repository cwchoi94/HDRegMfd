



#' @title Cross-Validation for High-Dimensional Quantile Linear Models
#' 
#' @description 
#' Implements cross-validation (CV) for high-dimensional quantile linear models based on 'AIC' or 'BIC'.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function 'QM_CV' in 'QM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{QM}}.
#' 
#' @inheritParams LM.CV
#' @inheritParams QM
#'
#' @return an \code{\link{QM}} object with the following components:
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
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
QM.CV = function(Xorg,Yorg,tau=0.5,h=NULL,kernel='Gaussian',cv.type='AIC',penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,
                 max.cv.iter=20,cv.threshold=1e-10,cv.const=2,phi0=1e-4,c.phi=1.1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.penalty(penalty)
  Check.kernel.QM(kernel)
  Check.cv.type(cv.type)
  
  if ((penalty=='SCAD') & (gamma<2)){
    gamma = 3.7
  } else if ((penalty=='MCP') & (gamma<1)){
    gamma = 3
  }
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  if(is.null(lambda.list)){lambda.list=c(0)}
  Ymu = 0
  Yspace = 'Euclid'
  
  # PCA for X
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  Xdims = sapply(X,ncol)
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(X,ncol)))}
  
  # compute the default value of h
  # In the CV function, the default value of h is computed in the 'QM_CV' function.
  if (is.null(h)){h = c(-1.0)}
  
  # Use QM_CV function to obtain the optimal parameters
  result = QM_CV(X,Yorg,lambda.list,Xdim.max.list,cv.type,tau,h,kernel,penalty,gamma,cv.const,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply QM with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  
  object = QM_each(X,Yorg,opt.lambda,opt.Xdim.max,tau,h,kernel,penalty,gamma,phi0,c.phi,max.iter,threshold)
  
  # compute other parameters
  Xdims = object$Xdims
  Xdims_cumul = c(0,cumsum(Xdims))
  
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime.opt.second = as.numeric(difftime(Sys.time(),opt.start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['pca']] = pca
  object[['tau']] = tau
  object[['kernel']] = kernel
  object[['h']] = object$h
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['iter.inner']] = object$iter.inner[1:object$iter]
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  object[['runtime.opt.second']] = runtime.opt.second
  class(object) = 'QM'
  
  return(object)
}




