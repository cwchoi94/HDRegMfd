




#' @title Cross-Validation for Oracle generalized Linear Models
#' 
#' @description 
#' Implements cross-validation (CV) for oracle generalized linear models based on 'AIC' or 'BIC'.
#' For a more detailed description of parameters, see \code{\link{GLM.oracle}}.
#' 
#' @inheritParams GLM.CV
#' @inheritParams GLM.oracle
#' 
#' @return a \code{\link{GLM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{link}{the link function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM.oracle.CV = function(Xorg,Yorg,link='binomial',proper.indices=NULL,cv.type='AIC',Xdim.max.list=NULL,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.link(link)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  R.list = c(100000)
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(Xoracle,ncol)))}
  
  
  # Use GLM_GCV function to obtain the optimal parameters
  result = GLM_CV(Xoracle,Yorg,lambda.list,Xdim.max.list,R.list,cv.type,
                  'LASSO',link,1,0,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM.oracle(Xorg,Yorg,link,proper.indices,opt.Xdim.max,1,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['cv.type']] = cv.type
  object[['runtime']] = runtime
  
  return(object)
}




#' @title Generalized Cross-Validation for Oracle Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements a generalized cross-validation (GCV) for oracle generalized linear models.
#' For a more detailed description of parameters, see \code{\link{GLM.oracle}}.
#' 
#' @inheritParams GLM.GCV
#' @inheritParams GLM.oracle.CV
#' 
#' @return a \code{\link{GLM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{link}{the link function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM.oracle.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,link='binomial',proper.indices=NULL,Xdim.max.list=NULL,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.link(link)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  R.list = c(1000)
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  Xnew = predict(pca,Xorgnew)
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  Xoracle.new = lapply(proper.indices,function(j){Xnew[[j]]})
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(Xoracle,ncol)))}
  
  
  # Use GLM_GCV function to obtain the optimal parameters
  result = GLM_GCV(Xoracle,Yorg,Xoracle.new,Yorgnew,lambda.list,Xdim.max.list,R.list,
                   'LASSO',link,1,0,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM.oracle(Xorg,Yorg,link,proper.indices,opt.Xdim.max,1,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}





#' @title Kfold Cross-Validation for Oracle generalized Linear Models
#' 
#' @description 
#' Implements a Kfold cross-validation (Kfold-CV) for oracle generalized linear models.
#' For a more detailed description of parameters, see \code{\link{GLM.oracle}}.
#' 
#' @inheritParams GLM.kfold
#' @inheritParams GLM.oracle.CV
#' 
#' @return a \code{\link{GLM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a \code{\link{PCA.manifold.list}} object.}
#'       \item{link}{the link function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
GLM.oracle.kfold = function(Xorg,Yorg,link='binomial',proper.indices=NULL,kfold=5,seed=NULL,Xdim.max.list=NULL,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.link(link)
  
  # define basic parameters
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  R.list = c(1000)
  
  ## compute Xdim.max.list when it doesn't exist
  if (is.null(Xdim.max.list)){
    pca = PCA.manifold.list(Xall)
    X = predict(pca,Xall)
    X.oracle = lapply(proper.indices,function(j){X[[j]]})
    Xdim.max.list = c(max(sapply(X,ncol)))
  }
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # data split and preprocessing
  test.indices.list = get.kfold.test.indices(n,kfold,seed)
  
  Xoracle.list = list()
  Xoracle.new.list = list()
  Y.list = list()
  Ynew.list = list()
  for (idx in 1:kfold){
    test.indices = test.indices.list[[idx]]
    data.split = split.data.org(Xall,Yall,test.indices)
    
    # compute Xorg, Xnew, LogY, LogYnew, and Ymu for each kfold data.
    pca = PCA.manifold.list(data.split$Xorg)
    Ymu = FrechetMean.manifold(data.split$Yorg,Yspace)
    
    X = predict(pca,data.split$Xorg)
    Xnew = predict(pca,data.split$Xnew)
    Xoracle.list[[idx]] = lapply(proper.indices,function(j){X[[j]]})
    Xoracle.new.list[[idx]] = lapply(proper.indices,function(j){Xnew[[j]]})
    Y.list[[idx]] = RieLog.manifold(Ymu,data.split$Yorg,Yspace)
    Ynew.list[[idx]] = RieLog.manifold(Ymu,data.split$Ynew,Yspace)
  }
  
  # Use GLM_kfold function defined in cpp
  result = GLM_Kfold(Xoracle.list,Y.list,Xoracle.new.list,Ynew.list,kfold,lambda.list,Xdim.max.list,R.list,
                     'LASSO',link,0,1,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max','R')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM.oracle(Xorg,Yorg,link,proper.indices,opt.Xdim.max,1,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}





