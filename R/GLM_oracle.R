### Implementation of linear regression using the knowledge of the nonzero index set.



#' @title High-dimensional linear regression for manifold-valued responses and covariates.
#' 
#' @description 
#' This function is based on a high-dimensional linear regression models with the knowledge of the nonzero indices.
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param Xdim.max a max dimension of \eqn{X_j}, real>0.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#'
#' @return an \code{\link{LM}} object.
#'    \describe{
#'       \item{pca}{an \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{an \eqn{m} vector of the Frechet mean of \eqn{Y}.}
#'       \item{beta}{a \eqn{P\times m} matrix of estimated beta, where \eqn{P=\sum_{j=1}^p K_j}.}
#'       \item{beta.each}{a \eqn{p} list of each \eqn{beta_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norm of each \eqn{beta_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of corresponding bases of \eqn{X_j}. Each basis is a \eqn{K_j}-by-\eqn{T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of tensor operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an indices of nonzero \code{beta_j}.}
#'       \item{runtime}{a running time.}
#'       \item{...}{other input parameters.}
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






