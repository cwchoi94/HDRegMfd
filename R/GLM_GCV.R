### Cross validation for high-dimensional linear regression


#' @title Generalized cross validation for high-dimensional linear regression
#' 
#' @description 
#' GCV for a LM object
#' It is based on CBS algorithm and uses a function LM_GCV in LM_GCV.cpp.
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Xnew a list of new covariate observations.
#' @param Ynew an \eqn{n'\times} new response observations.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param lambda.list a vector of lambda.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param R.list a vector of constrained bound.
#' @param penalty a method of penalty. It should be one of 'LASSO', 'SCAD', or 'MCP'.
#' @param phi a parameter in computing ADMM-MM algorithm for the majorized objective function, default 1.
#' @param gamma a parameter for SCAD (3.7) or MCP (3), parentheses: default value.
#' @param max.cv.iter a number of maximum CV iterations, default 20.
#' @param cv.threshold a parameter to modify the computation error in CV, default 1e-10.
#' @param eta a parameter in computing ADMM-MM algorithm for the proximal norm square, default 1e-3.
#' @param max.iter a maximum iteration, default 500.
#' @param threshold an algorihtm convergence threshold, default 1e-10.
#'
#' @return an \code{\link{LM}} object.
#'    \describe{
#'       \item{parameter.list}{a list of optimal parameters for each CV step.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time}
#'       \item{...}{see \code{\link{LM}}.}
#' }
#' @export
GLM.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,lambda.list,Xdim.max.list,R.list,penalty='LASSO',link='binomial',phi=1,gamma=0,
                  max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  
  # PCA for X
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  Xnew = predict(pca,Xorgnew)
  
  # Use GLM_GCV function to obtain the optimal parameters
  result = GLM_GCV(X,Yorg,Xnew,Yorgnew,lambda.list,Xdim.max.list,R.list,
                   penalty,link,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = GLM_each(X,Yorg,opt.lambda,opt.Xdim.max,opt.R,penalty,link,phi,gamma,eta,max.iter,threshold)
  
  # compute other parameters
  Xdims = object$Xdims
  Xdims_cumul = c(0,cumsum(Xdims))
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  beta.each = lapply(1:p,function(j){object$beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
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
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  class(object) = 'GLM'
  
  return(object)
}




