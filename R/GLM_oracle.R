


#' @title Oracle generalized linear regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Estimate Hilbert-Schmidt operators using an ADMM-based algorithm with the knowledge of the nonzero index set \eqn{\mathcal{S}=\{j:\mathfrak{B}_j\neq0}}.
#' 
#' @inheritParams GLM
#' 
#' @param proper.indices an index set \eqn{\mathcal{S}=\{1\le j\le p : \mathfrak{B}_j\neq0\}}.
#'
#' @return a \code{\link{GLM}} object with the following compnents:
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
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM.oracle = function(Xorg,Yorg,Xdim.max=100,proper.indices=NULL,link='binomial',phi=1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.link(link)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  Xdims.oracle = sapply(Xoracle,ncol)
  Xdims_cumul.oracle = c(0,cumsum(Xdims.oracle))
  
  # apply GLM_each function in cpp
  object = GLM_each(Xoracle,Yorg,0,Xdim.max,1e10,'LASSO',link,phi,0,1e-3,max.iter,threshold)
  
  # compute oracle estimator
  beta.oracle = object$beta
  beta = matrix(0,sum(Xdims),ncol(Yorg))
  for (i in 1:length(proper.indices)){
    j = proper.indices[i]
    beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],] = beta.oracle[(Xdims_cumul.oracle[i]+1):Xdims_cumul.oracle[i+1],]
  }
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['beta']] = beta
  object[['pca']] = pca
  object[['link']] = link
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['Xdim.max']] = Xdim.max
  object[['runtime']] = runtime
  class(object) = 'GLM'
  
  return(object)
}






