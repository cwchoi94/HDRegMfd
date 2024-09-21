


# Compute loss for each Xdim.max
get.loss.LM.oracle = function(X,LogY,Xnew,LogYnew,Ymu,Yspace,Xdim.max,proper.indices){
  X = reduce.dimension(X,Xdim.max)
  Xnew = reduce.dimension(Xnew,Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  # compute oracle LSE
  beta = compute_beta(X,LogY,proper.indices)
  
  LogYhat = Xnew %*% beta
  rmse = vector.norm(LogYnew-LogYhat,Ymu,Yspace,'L2')/sqrt(n2)
  
  return(rmse)
}



#' @title Cross-Validation for Oracle Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements cross-validation (CV) for oracle Hilbert-Schmidt linear models based on 'AIC' or 'BIC'.
#' For a more detailed description of parameters, see \code{\link{LM.oracle}}.
#' 
#' @inheritParams LM.CV
#' @inheritParams LM.oracle
#' 
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.oracle.CV = function(Xorg,Yorg,Yspace,Xdim.max.list,proper.indices=NULL,cv.type='AIC',max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.cv.type(cv.type)
  
  # define basic parameters
  n = nrow(Yorg)
  p = Xorg[['p']]
  inner = eval(parse(text=paste0('inner.each.',Yspace)))
  
  lambda.list = c(0)
  R.list = c(100000)
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # PCA for X
  pca = PCA.manifold.list(Xorg)
  X = predict(pca,Xorg)
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  
  # projection of Yorg and Yorgnew to the tangent space
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  
  # Use LM_CV function to obtain the optimal parameters
  result = LM_CV(Xoracle,LogY,Ymu,Yspace,lambda.list,Xdim.max.list,R.list,cv.type,
                 'LASSO',1,0,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply LM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  X = reduce.dimension(X,opt.Xdim.max)
  beta = compute_beta(X,LogY,proper.indices)
  
  
  # compute other parameters
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
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
  object[['parameter.list']] = parameter.list
  object[['loss.list']] = loss.list
  object[['cv.type']] = cv.type
  object[['Xdim.max']] = opt.Xdim.max
  object[['runtime']] = runtime
  class(object) = 'LM'
  
  return(object)
}



#' @title Generalized Cross-Validation for Oracle Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements a generalized cross-validation (GCV) for oracle Hilbert-Schmidt linear models.
#' For a more detailed description of parameters, see \code{\link{LM.oracle}}.
#' 
#' @inheritParams LM.GCV
#' @inheritParams LM.oracle.CV
#'
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.oracle.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,Yspace,Xdim.max.list,proper.indices=NULL){
  
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
  Xnew = predict(pca,Xorgnew)
  
  # compute LogY
  Ymu = FrechetMean.manifold(Yorg,Yspace)
  LogY = RieLog.manifold(Ymu,Yorg,Yspace)
  LogYnew = RieLog.manifold(Ymu,Yorgnew,Yspace)
  
  
  # GCV for Xdim.max
  max.Xdim.max = max(sapply(X,ncol))
  Xdim.max.list[Xdim.max.list>max.Xdim.max] = max.Xdim.max
  Xdim.max.list = unique(Xdim.max.list)
  r2 = length(Xdim.max.list)
  
  loss.list = sapply(1:r2,function(i){get.loss.LM.oracle(X,LogY,Xnew,LogYnew,Ymu,Yspace,Xdim.max.list[i],proper.indices)})
  
  # compute optimal beta
  opt.Xdim.max = Xdim.max.list[which.min(loss.list)]
  X = reduce.dimension(X,opt.Xdim.max)
  beta = compute_beta(X,LogY,proper.indices)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
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
  object[['Xdim.max']] = opt.Xdim.max
  object[['runtime']] = runtime
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  class(object) = 'LM'
  
  return(object)
}




#' @title Kfold Cross-Validation for Oracle Hilbert-Schmidt Linear Models
#' 
#' @description 
#' Implements a Kfold cross-validation (Kfold-CV) for oracle Hilbert-Schmidt linear models.
#' For a more detailed description of parameters, see \code{\link{LM.oracle}}.
#' 
#' @inheritParams LM.kfold
#' @inheritParams LM.oracle.CV
#' 
#' @return an \code{\link{LM}} object with the following components:
#'    \describe{
#'       \item{pca}{a 'PCA.manifold.list' object, see \code{\link{PCA.manifold.list}}.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{beta}{a \eqn{L_+^{*} \times m} matrix of estimated \eqn{\bm{\beta}}, where \eqn{L_+^{*}=\sum_{j=1}^p L_j^*} and \eqn{m} is the intrinsic dimension of \eqn{T_{\mu_Y}\mathcal{M}_Y}.}
#'       \item{beta.each}{a \eqn{p} list of \eqn{L_j^*\times m} matrices of \eqn{\bm{\beta}_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of norms of \eqn{\bm{\beta}_j}.}
#'       \item{beta.vectors}{a \eqn{p} list of orthonormal bases of \eqn{X_j} obtained by \code{\link{PCA.manifold.list}}. Each basis is an \eqn{L_j^*\times T_j} matrix.}
#'       \item{beta.tensor}{a \eqn{p} list of estimated Hilbert-Schmidt operators, see \code{\link{make.tensor}}.}
#'       \item{proper.indices}{an index set an index set \eqn{\mathcal{S}=\{1\le j\le p : {\mathfrak{B}}_j\neq0\}}.}
#'       \item{parameter.list}{a list of optimal parameters for each CV update.}
#'       \item{loss.list}{a list of loss for each CV update.}
#'       \item{runtime}{the running time.}
#'       \item{...}{other parameters.}
#' }
#' @export
LM.oracle.kfold = function(Xorg,Yorg,Yspace,kfold,Xdim.max.list,proper.indices=NULL,seed=NULL){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.manifold(Yspace)
  
  # define basic parameters
  Xall = Xorg
  Yall = Yorg
  
  n = nrow(Yall)
  p = Xall[['p']]
  
  proper.indices = get.proper.indices(proper.indices,p)
  
  # data split and preprocessing
  test.indices.list = get.kfold.test.indices(n,kfold,seed)
  
  Xorg.list = list()
  Xnew.list = list()
  LogY.list = list()
  LogYnew.list = list()
  Ymu.list = list()
  for (idx in 1:kfold){
    test.indices = test.indices.list[[idx]]
    data.split = split.data.org(Xall,Yall,test.indices)
    
    # compute Xorg, Xnew, LogY, LogYnew, and Ymu for each kfold data.
    pca = PCA.manifold.list(data.split$Xorg)
    Ymu = FrechetMean.manifold(data.split$Yorg,Yspace)
    
    Xorg.list[[idx]] = predict(pca,data.split$Xorg)
    Xnew.list[[idx]] = predict(pca,data.split$Xnew)
    LogY.list[[idx]] = RieLog.manifold(Ymu,data.split$Yorg,Yspace)
    LogYnew.list[[idx]] = RieLog.manifold(Ymu,data.split$Ynew,Yspace)
    Ymu.list[[idx]] = Ymu
  }
  
  
  # kfold for Xdim.max
  max.Xdim.max = max(sapply(1:kfold,function(i){max(sapply(Xorg.list[[i]],ncol))}))
  Xdim.max.list[Xdim.max.list>max.Xdim.max] = max.Xdim.max
  Xdim.max.list = unique(Xdim.max.list)
  r2 = length(Xdim.max.list)
  
  loss.list = matrix(0,r2,kfold)
  for (idx in 1:kfold){
    Xorg = Xorg.list[[idx]]
    Xnew = Xnew.list[[idx]]
    LogY = LogY.list[[idx]]
    LogYnew = LogYnew.list[[idx]]
    Ymu = Ymu.list[[idx]]
    
    loss = sapply(1:r2,function(i){get.loss.LM.oracle(Xorg,LogY,Xnew,LogYnew,Ymu,Yspace,Xdim.max.list[i],proper.indices)})
    loss.list[,idx] = unlist(loss)
  }
    
  loss = rowMeans(loss.list)
  opt.Xdim.max = Xdim.max.list[which.min(loss)]
  
  # fit opt model
  object = LM.oracle(Xall,Yall,Yspace,opt.Xdim.max,proper.indices)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}




