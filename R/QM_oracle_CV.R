




#' @title Cross-Validation for Oracle Quantile Linear Models
#' 
#' @description 
#' Implements cross-validation (CV) for oracle quantile linear models based on 'AIC' or 'BIC'.
#' For a more detailed description of parameters, see \code{\link{QM.oracle}}.
#' 
#' @inheritParams QM.CV
#' @inheritParams QM.oracle
#' 
#' @return a \code{\link{QM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{tau}{the quantile.}
#'       \item{kernel}{the kernel function for smoothing a check function.}
#'       \item{h}{the bandwidth for smoothing a check function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{iter.inner}{a vector of iteration numbers for each sub-iteration.}
#'       \item{runtime}{the running time (HH:MM:SS).}
#'       \item{runtime.second}{the running time (second).}
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
QM.oracle.CV = function(Xorg,Yorg,tau=0.5,h=NULL,kernel='Gaussian',proper.indices=NULL,cv.type='AIC',Xdim.max.list=NULL,
                        max.cv.iter=20,cv.threshold=1e-10,cv.const=2,phi0=1e-4,c.phi=1.1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.kernel.QM(kernel)
  Check.cv.type(cv.type)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # compute the default value of h
  # In the CV function, the default value of h is computed in the 'QM_CV' function.
  if (is.null(h)){h = c(-1.0)}
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(Xoracle,ncol)))}
  
  
  # Use QM_CV function to obtain the optimal parameters
  result = QM_CV(Xoracle,Yorg,lambda.list,Xdim.max.list,cv.type,tau,h,kernel,
                  'LASSO',0,cv.const,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply QM with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.Xdim.max = result$opt.Xdim.max
  
  object = QM.oracle(Xorg,Yorg,proper.indices,tau,h,kernel,opt.Xdim.max,phi0,c.phi,max.iter,threshold)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime.opt.second = as.numeric(difftime(Sys.time(),opt.start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['cv.type']] = cv.type
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  object[['runtime.opt.second']] = runtime.opt.second
  
  return(object)
}




#' @title Generalized Cross-Validation for Oracle Quantile Linear Models
#' 
#' @description 
#' Implements a generalized cross-validation (GCV) for oracle quantile linear models.
#' For a more detailed description of parameters, see \code{\link{QM.oracle}}.
#' 
#' @inheritParams QM.GCV
#' @inheritParams QM.oracle.CV
#' 
#' @return a \code{\link{QM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{tau}{the quantile.}
#'       \item{kernel}{the kernel function for smoothing a check function.}
#'       \item{h}{the bandwidth for smoothing a check function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{iter.inner}{a vector of iteration numbers for each sub-iteration.}
#'       \item{runtime}{the running time (HH:MM:SS).}
#'       \item{runtime.second}{the running time (second).}
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
QM.oracle.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,proper.indices=NULL,tau=0.5,h=NULL,kernel='Gaussian',Xdim.max.list=NULL,
                         max.cv.iter=20,cv.threshold=1e-10,phi0=1e-4,c.phi=1.1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.kernel.QM(kernel)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  Xnew = predict(pca,Xorgnew)
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  Xoracle.new = lapply(proper.indices,function(j){Xnew[[j]]})
  if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(Xoracle,ncol)))}
  
  # compute the default value of h
  # In the CV function, the default value of h is computed in the 'QM_CV' function.
  if (is.null(h)){h = c(-1.0)}
  
  # Use QM_GCV function to obtain the optimal parameters
  result = QM_GCV(Xoracle,Yorg,Xoracle.new,Yorgnew,lambda.list,Xdim.max.list,
                  tau,h,kernel,'LASSO',0,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.Xdim.max = result$opt.Xdim.max
  
  object = QM.oracle(Xorg,Yorg,proper.indices,tau,h,kernel,opt.Xdim.max,phi0,c.phi,max.iter,threshold)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime.opt.second = as.numeric(difftime(Sys.time(),opt.start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  object[['runtime.opt.second']] = runtime.opt.second
  
  return(object)
}





#' @title Kfold Cross-Validation for Oracle quantile Linear Models
#' 
#' @description 
#' Implements a Kfold cross-validation (Kfold-CV) for oracle quantile linear models.
#' For a more detailed description of parameters, see \code{\link{QM.oracle}}.
#' 
#' @inheritParams QM.kfold
#' @inheritParams QM.oracle.CV
#' 
#' @return a \code{\link{QM}} object with the following compnents:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{tau}{the quantile.}
#'       \item{kernel}{the kernel function for smoothing a check function.}
#'       \item{h}{the bandwidth for smoothing a check function.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta0}{an \eqn{m} vector of the intercept constant.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm0}{the norm of \eqn{{\beta}_0^*}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an estimated index set an index set \eqn{\mathcal{S}=\{1\le j\le p : \hat{\mathfrak{B}}_j\neq0\}}.}
#'       \item{iter.inner}{a vector of iteration numbers for each sub-iteration.}
#'       \item{runtime}{the running time (HH:MM:SS).}
#'       \item{runtime.second}{the running time (second).}
#'       \item{runtime.opt.second}{the running time with the optimal parmaters (second).}
#'       \item{...}{other parameters.}
#' }
#' @export
QM.oracle.kfold = function(Xorg,Yorg,proper.indices=NULL,kfold=5,seed=NULL,tau=0.5,h=NULL,kernel='Gaussian',Xdim.max.list=NULL,
                           max.cv.iter=20,cv.threshold=1e-10,phi0=1e-4,c.phi=1.1,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.kernel.QM(kernel)
  
  # define basic parameters
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  
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
  
  # compute the default value of h
  # In the CV function, the default value of h is computed in the 'QM_CV' function.
  if (is.null(h)){h = c(-1.0)}
  
  # Use QM_kfold function defined in cpp
  result = QM_Kfold(Xoracle.list,Y.list,Xoracle.new.list,Ynew.list,kfold,lambda.list,Xdim.max.list,
                    tau,h,kernel,'LASSO',0,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
  colnames(parameter.list) = c('lambda','Xdim.max')
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  # apply GLM with the optimal parameters
  opt.start.time = Sys.time()
  
  opt.Xdim.max = result$opt.Xdim.max
  
  object = QM.oracle(Xorg,Yorg,proper.indices,tau,h,kernel,opt.Xdim.max,phi0,c.phi,max.iter,threshold)
  
  runtime.second = as.numeric(difftime(Sys.time(),start.time,units='secs'))
  runtime.opt.second = as.numeric(difftime(Sys.time(),opt.start.time,units='secs'))
  runtime = hms::hms(round(runtime.second))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  object[['runtime.second']] = runtime.second
  object[['runtime.opt.second']] = runtime.opt.second
  
  return(object)
}





