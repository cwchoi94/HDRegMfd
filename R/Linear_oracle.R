### Implementation of linear regression using the knowledge of the nonzero index set.



get.proper.indices = function(proper.indices,p){
  if (is.null(proper.indices)){
    proper.indices = 1:p
  } else{
    proper.indices = proper.indices[proper.indices<=p]
  }
  
  return(proper.indices)
}


# compute LSE
compute_beta = function(X,Y,proper.indices){
  
  n = nrow(Y)
  m = ncol(Y)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  X = do.call(cbind,X)
  P = ncol(X)
  W = sapply(1:P,function(j){(t(X[,j]) %*% X[,j] / n)})
  X = X %*% diag(W^(-1/2))
  
  proper.indices.all = unlist(sapply(proper.indices,function(j){(Xdims_cumul[j]+1):Xdims_cumul[j+1]}))
  X = X[,proper.indices.all]
  XX = t(X) %*% X / n
  XY = t(X) %*% Y / n
  
  beta = matrix(0,P,m)
  beta[proper.indices.all,] = solve(XX) %*% XY
  beta = diag(W^(-1/2)) %*% beta
  
  return(beta)
}



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
LM.oracle = function(Xorg,Yorg,Yspace,Xdim.max=100,proper.indices=NULL){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.manifold(Yspace)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  X = reduce.dimension(X,Xdim.max)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  
  
  # compute oracle LSE
  beta = compute_beta(X,LogY,proper.indices)
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object = list()
  
  object[['beta']] = beta
  object[['pca']] = pca
  object[['Ymu']] = Ymu
  object[['Yspace']] = Yspace
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['Xdim.max']] = Xdim.max
  object[['runtime']] = runtime
  class(object) = 'LM'
  
  return(object)
}






