### ADMM|MM algorithm update for high-dimensional generalized linear regression
### developed only for LASSO, SCAD, MCP





#' @title High-dimensional linear regression for manifold-valued responses and covariates.
#' 
#' @description 
#' Compute tensor product operators by ADMM-MM algorithm.
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param lambda a penalty constant, real>0.
#' @param Xdim.max a max dimension of \eqn{X_j}, real>0.
#' @param R a constrained bound, real>0.
#' @param penalty a method of penalty. It should be one of 'LASSO', 'SCAD', or 'MCP'.
#' @param phi a parameter in computing ADMM-MM algorithm for the majorized objective function, default 1.
#' @param gamma a parameter for SCAD (3.7) or MCP (3), parentheses: default value.
#' @param eta a parameter in computing ADMM-MM algorithm for the proximal norm square, default 1e-3.
#' @param max.iter a maximum iteration, default 500.
#' @param threshold an algorihtm convergence threshold, default 1e-10.
#'
#' @return a 'LM' object.
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{Ymu}{an \eqn{m} vector of the Frechet mean of \eqn{Y}.}
#'       \item{beta}{a \eqn{P\times m} matrix of estimated beta, where \eqn{P=\sum_{j=1}^p K_j}.}
#'       \item{beta.each}{a \eqn{p} list of each \eqn{beta_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norm of each \eqn{beta_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of corresponding bases of \eqn{X_j}. Each basis is a \eqn{K_j\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of tensor operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{a indices of nonzero \code{beta_j}.}
#'       \item{runtime}{a running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM = function(Xorg,Yorg,lambda=0.1,Xdim.max=100,R=100,penalty='LASSO',link='binomial',phi=1,gamma=0,
               eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
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



#' @title Prediction function for a LM object
#' 
#' @description 
#' Predict \eqn{\hat{Y}_{new}} for given \eqn{X_{new}}.
#' 
#' @param object a \code{\link{LM}} class object.
#' @param Xnew a \eqn{p} list of new observations.
#'
#' @return an \eqn{n'\times m} matrix of \eqn{\hat{Y}_{new}}.
#' @export
predict.GLM = function(object,Xnew,is.inv.link=TRUE){
  Xnew = predict.PCA.manifold.list(object$pca,Xnew)
  Xnew = reduce.dimension(Xnew,object$Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  theta = Xnew %*% object$beta + matrix(rep(object$beta0),nrow=nrow(theta))
  
  if(is.inv.link){
    Ymu = Inv_Link(theta,object$link)
    return (Ymu)
  }else{
    return(theta)
  }
  
}



