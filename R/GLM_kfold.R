


#' @title Kfold Cross-Validation for High-Dimensional generalized Linear Models
#' 
#' @description 
#' Implements Kfold cross-validation (Kfold-CV) for high-dimensional generalized linear models.
#' The CV process is based on the coordinate-wise variable selection and is implemented using a function 'GLM_CV' in 'GLM_CV.cpp'.
#' For a more detailed description of parameters, see \code{\link{GLM}}.
#' 
#' @inheritParams LM.kfold
#' @inheritParams GLM
#'
#' @return an \code{\link{GLM}} object with the following components:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{link}{the link function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta0.norm}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time (HH:MM:SS).}
#'       \item{runtime.second}{the running time (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM.kfold = function(Xorg,Yorg,link='binomial',kfold=5,seed=NULL,penalty='LASSO',gamma=0,lambda.list=NULL,Xdim.max.list=NULL,R.list=NULL,
                     phi=1,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  if(is.null(lambda.list)){lambda.list=c(0)}
  if(is.null(R.list)){R.list=c(100)}
  Ymu = 0
  Yspace = 'Euclid'
  
  ## compute Xdim.max.list when it doesn't exist
  if (is.null(Xdim.max.list)){
    pca = PCA.manifold.list(Xall)
    X = predict(pca,Xall)
    Xdim.max.list = c(max(sapply(X,ncol)))
  }
  
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
                     penalty,link,gamma,phi,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.lambda = result$opt.lambda
  opt.Xdim.max = result$opt.Xdim.max
  opt.R = result$opt.R
  
  object = GLM(Xall,Yall,link,penalty,gamma,opt.lambda,opt.Xdim.max,opt.R,phi,eta,max.iter,threshold)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  
  return(object)
}


