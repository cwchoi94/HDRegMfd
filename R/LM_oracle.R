

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



#' @title Oracle Hilbert-Schmidt linear regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Estimate Hilbert-Schmidt operators using the knowledge of the nonzero index set \eqn{\mathcal{S}=\{j:\mathfrak{B}_j\neq0\}}.
#' 
#' @inheritParams LM
#' 
#' @param proper.indices an index set \eqn{\mathcal{S}=\{1\le j\le p : \mathfrak{B}_j\neq0\}}.
#'
#' @return a \code{\link{LM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.oracle = function(Xorg,Yorg,Yspace,proper.indices=NULL,Xdim.max=100){
  
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
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
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
  object[['runtime.second']] = runtime.second
  class(object) = 'LM'
  
  return(object)
}






