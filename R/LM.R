### ADMM|MM algorithm update for high-dimensional linear regression
### developed only for LASSO, SCAD, MCP




#' @title Demension reduction function for a list of X.
#' 
#' @description 
#' Change the dimension of \eqn{X_j} at most \eqn{min(d_j,Xdim.max)}.
#' 
#' @param X a \eqn{p} list of matrices. Each \eqn{X_j} is a \eqn{n\times d_j} matrix.
#' @param Xdim.max a max dimension of \eqn{X_j}, int>0.
#' @param margin the direction which the dimension reduction applied over.
#' 
#' @return a \eqn{p} list of \eqn{X_j'}.
#' @export
reduce.dimension = function(X,Xdim.max=200,margin=1){
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
LM = function(Xorg,Yorg,Yspace,lambda=0.1,Xdim.max=100,R=100,penalty='LASSO',phi=1,gamma=0,
              eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  object = LM_each(X,LogY,Ymu,inner,lambda,Xdim.max,R,penalty,phi,gamma,eta,max.iter,threshold)
  
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
predict.LM = function(object,Xnew){
  Xnew = predict.PCA.manifold.list(object$pca,Xnew)
  Xnew = reduce.dimension(Xnew,object$Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  LogYhat = Xnew %*% object$beta
  Yhat = RieExp.manifold(object[['Ymu']],LogYhat,object[['Yspace']])
  return(Yhat)
}



