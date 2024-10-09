

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
  
  Yorg = vec.to.mat(Yall[-test.indices,,drop=FALSE])
  Ynew = vec.to.mat(Yall[test.indices,,drop=FALSE])
  data.split = list(Xorg=Xorg,Yorg=Yorg,Xnew=Xnew,Ynew=Ynew)
  return(data.split)
}



#' @title Kfold Cross-Validation for High-Dimensional Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements Kfold cross-validation (Kfold-CV) for high-dimensional Hilbert-Schmidt linear models.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function 'LM_CV' in 'LM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{LM}}.
#' 
#' @inheritParams LM.CV
#' 
#' @param kfold the number of kfolds (default: 5).
#' @param seed a random seed for kfolds. When it is not provided, the data will be partitioned in index order.
#' 
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.kfold = function(Xorg,Yorg,Yspace,kfold=5,seed=NULL,penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,R.list=NULL,
                    phi=1,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  
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
  result = LM_Kfold(Xorg.list,LogY.list,Xnew.list,LogYnew.list,Ymu.list,Yspace,kfold,
                    lambda.list,Xdim.max.list,R.list,penalty,gamma,phi,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = LM(Xall,Yall,Yspace,opt.lambda,opt.Xdim.max,opt.R,penalty,gamma,phi,eta,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}


