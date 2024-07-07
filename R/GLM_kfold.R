


#' @title Kfold cross validation for high-dimensional linear regression for manifold-valued data
#' 
#' @description 
#' Kfold CV for an \code{\link{LM}} function.
#' It is based on CBS algorithm and uses a function LM_kfold in LM_kfold.cpp.
#' 
#' @param Xall a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yall an \eqn{n\times m} response matrix.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param kfold a number of kfold CV, int>0.
#' @param lambda.list a vector of lambda.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param R.list a vector of constrained bound.
#' @param penalty a method of penalty. It should be one of 'LASSO', 'SCAD', or 'MCP'.
#' @param phi a parameter in computing ADMM-MM algorithm for the majorized objective function, default 1.
#' @param gamma a parameter for SCAD (3.7) or MCP (3), parentheses: default value.
#' @param seed a random seed, int>0, default: non-random (NULL).
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
GLM.kfold = function(Xall,Yall,kfold,lambda.list,Xdim.max.list,R.list,penalty='LASSO',link='binomial',phi=1,gamma=0,seed=NULL,
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
  n = nrow(Yall)
  p = Xall[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  # data split and preprocessing
  test.indices.list = get.kfold.test.indices(n,kfold,seed)
  
  Xorg.list = list()
  Xnew.list = list()
  Y.list = list()
  Ynew.list = list()
  for (i in 1:kfold){
    test.indices = test.indices.list[[i]]
    data.split = split.data.org(Xall,Yall,test.indices)
    
    # compute Xorg, Xnew, LogY and LogYnew for each kfold data.
    pca = PCA.manifold.list(data.split$Xorg)
    Ymu = FrechetMean.manifold(data.split$Yorg,Yspace)
    
    Xorg.list[[i]] = predict(pca,data.split$Xorg)
    Xnew.list[[i]] = predict(pca,data.split$Xnew)
    Y.list[[i]] = RieLog.manifold(Ymu,data.split$Yorg,Yspace)
    Ynew.list[[i]] = RieLog.manifold(Ymu,data.split$Ynew,Yspace)
  }
  
  # Use LM_kfold function defined in cpp
  result = GLM_Kfold(Xorg.list,Y.list,Xnew.list,Ynew.list,kfold,lambda.list,Xdim.max.list,R.list,
                     penalty,link,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = GLM(Xall,Yall,opt.lambda,opt.Xdim.max,opt.R,penalty,link,phi,gamma,eta,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}


