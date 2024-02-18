### Cross validation for high-dimensional linear regression


# partition all indices into kfold indices
get.kfold.test.indices = function(n,kfold,seed=NULL){
  
  all.indices = 1:n
  if (!is.null(seed)){
    set.seed(seed)
    all.indices = sample(all.indices,n)
  }
  
  test.indices.list = list()
  for (i in 1:kfold){
    idx1 = floor(n/kfold*(i-1))+1
    idx2 = floor(n/kfold*i)
    
    test.indices.list[[i]] = all.indices[idx1:idx2]
  }
  
  return(test.indices.list)
}


# split data by test.indices
split.data.org = function(Xall,Yall,test.indices){
  p = Xall[['p']]
  spaces = Xall[['spaces']]
  
  Xorg = list()
  Xnew = list()
  for (j in 1:p){
    Xorg[[j]] = as.matrix(Xall[[j]][-test.indices,])
    Xnew[[j]] = as.matrix(Xall[[j]][test.indices,])
  }
  Xorg[['p']] = p
  Xnew[['p']] = p
  Xorg[['spaces']] = spaces
  Xnew[['spaces']] = spaces
  
  Yorg = vec.to.mat(Yall[-test.indices,])
  Ynew = vec.to.mat(Yall[test.indices,])
  data.split = list(Xorg=Xorg,Yorg=Yorg,Xnew=Xnew,Ynew=Ynew)
  return(data.split)
}



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
#' @param phi a parameter in computing ADMM-MM algorithm for the majorized objective function, default 1.
#' @param penalty a method of penalty. It should be one of 'LASSO', 'SCAD', or 'MCP'.
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
LM.kfold = function(Xall,Yall,Yspace,kfold,lambda.list,Xdim.max.list,R.list,phi=1,penalty,gamma=0,seed=NULL,
                    max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  n = nrow(Yall)
  p = Xall[['p']]
  inner = eval(parse(text=paste0('inner.each.',Yspace)))
  
  # data split and preprocessing
  test.indices.list = get.kfold.test.indices(n,kfold,seed)
  
  Xorg.list = list()
  Xnew.list = list()
  LogY.list = list()
  LogYnew.list = list()
  Ymu.list = list()
  for (i in 1:kfold){
    test.indices = test.indices.list[[i]]
    data.split = split.data.org(Xall,Yall,test.indices)
    
    # compute Xorg, Xnew, LogY, LogYnew, and Ymu for each kfold data.
    pca = PCA.manifold.list(data.split$Xorg)
    Ymu = FrechetMean.manifold(data.split$Yorg,Yspace)
    
    Xorg.list[[i]] = predict(pca,data.split$Xorg)
    Xnew.list[[i]] = predict(pca,data.split$Xnew)
    LogY.list[[i]] = RieLog.manifold(Ymu,data.split$Yorg,Yspace)
    LogYnew.list[[i]] = RieLog.manifold(Ymu,data.split$Ynew,Yspace)
    Ymu.list[[i]] = Ymu
  }
  
  # Use LM_kfold function defined in cpp
  result = LM_Kfold(Xorg.list,LogY.list,Xnew.list,LogYnew.list,Ymu.list,inner,kfold,
                    lambda.list,Xdim.max.list,R.list,phi,penalty,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = LM(Xall,Yall,Yspace,opt.lambda,opt.Xdim.max,opt.R,phi,penalty,gamma,eta,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}


