### Cross validation for LM.oracle


#' GCV for an LM.oracle function
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Xnew a new list of manifold-valued covariates.
#' @param Ynew a new \eqn{n'\times m} response matrix.
#' @param Yspace an underlying space of Yorg and Ynew.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#'
#' @return an \eqn{\code{LM}} object
#'    \describe{
#'       \item{Xdim.max.list}{a vector of Xdim.max.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time.}
#'       \item{...}{see \code{\link{LM.oracle}}.}
#' }
#' @export
GLM.oracle.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,Xdim.max.list,proper.indices=NULL,link='binomial',max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  
  
  # Use GLM_GCV function to obtain the optimal parameters
  result = GLM_GCV(Xoracle,Yorg,Xoracle.new,Yorgnew,lambda.list,Xdim.max.list,R.list,
                   penalty,link,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM_each(Xoracle,Yorg,0,opt.Xdim.max,1000,'LASSO',link,phi,0,1e-3,max.iter,threshold)
  
  
  # compute oracle estimator
  X = reduce.dimension(X,opt.Xdim.max)
  Xoracle = reduce.dimension(Xoracle,opt.Xdim.max)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  Xdims.oracle = sapply(Xoracle,ncol)
  Xdims_cumul.oracle = c(0,cumsum(Xdims.oracle))
  
  beta.oracle = object$beta
  beta = matrix(0,sum(Xdims),ncol(Yorg))
  for (i in 1:length(proper.indices)){
    j = proper.indices[i]
    beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],] = beta.oracle[(Xdims_cumul.oracle[i]+1):Xdims_cumul.oracle[i+1],]
  }
  
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['beta']] = beta
  object[['pca']] = pca
  object[['link']] = link
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['Xdim.max']] = opt.Xdim.max
  object[['runtime']] = runtime
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  class(object) = 'GLM'
  
  return(object)
}





#' Kfold cross validation for a GLM.oracle function
#' 
#' @param Xall a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yall an \eqn{n\times m} response matrix.
#' @param Yspace an underlying space of Yorg and Ynew.
#' @param kfold a number of kfold CV, int>0.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#' @param seed a random seed, int>0, default: non-random (NULL).
#'
#' @return an \eqn{\code{LM}} object
#'    \describe{
#'       \item{Xdim.max.list}{a vector of Xdim.max.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time.}
#'       \item{...}{see \code{\link{LM.oracle}}.}
#' }
#' @export
GLM.oracle.kfold = function(Xall,Yall,kfold,Xdim.max.list,proper.indices=NULL,link='binomial',seed=NULL,max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.link(link)
  
  # define basic parameters
  n = nrow(Yall)
  p = Xall[['p']]
  Ymu = 0
  Yspace = 'Euclid'
  
  lambda.list = c(0)
  R.list = c(1000)
  
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
                     penalty,link,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM.oracle(Xall,Yall,opt.Xdim.max,proper.indices,link,phi,max.iter,threshold)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}



#' GCV for an LM.oracle function
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Xnew a new list of manifold-valued covariates.
#' @param Ynew a new \eqn{n'\times m} response matrix.
#' @param Yspace an underlying space of Yorg and Ynew.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#'
#' @return an \eqn{\code{LM}} object
#'    \describe{
#'       \item{Xdim.max.list}{a vector of Xdim.max.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time.}
#'       \item{...}{see \code{\link{LM.oracle}}.}
#' }
#' @export
GLM.oracle.CV = function(Xorg,Yorg,Xdim.max.list,proper.indices=NULL,cv.type='AIC',link='binomial',max.cv.iter=20,cv.threshold=1e-10,eta=1e-3,max.iter=500,threshold=1e-10){
  
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
  
  Xoracle = lapply(proper.indices,function(j){X[[j]]})
  
  
  # Use GLM_GCV function to obtain the optimal parameters
  result = GLM_CV(Xoracle,Yorg,lambda.list,Xdim.max.list,R.list,cv.type,
                  penalty,link,phi,gamma,max.cv.iter,cv.threshold)
  
  parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),]
  loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]
  
  
  # apply GLM with the optimal parameters
  opt.Xdim.max = result$opt.Xdim.max
  
  object = GLM_each(Xoracle,Yorg,0,opt.Xdim.max,1000,'LASSO',link,phi,0,1e-3,max.iter,threshold)
  
  
  # compute oracle estimator
  X = reduce.dimension(X,opt.Xdim.max)
  Xoracle = reduce.dimension(Xoracle,opt.Xdim.max)
  
  Xdims = sapply(X,ncol)
  Xdims_cumul = c(0,cumsum(Xdims))
  
  Xdims.oracle = sapply(Xoracle,ncol)
  Xdims_cumul.oracle = c(0,cumsum(Xdims.oracle))
  
  beta.oracle = object$beta
  beta = matrix(0,sum(Xdims),ncol(Yorg))
  for (i in 1:length(proper.indices)){
    j = proper.indices[i]
    beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],] = beta.oracle[(Xdims_cumul.oracle[i]+1):Xdims_cumul.oracle[i+1],]
  }
  
  
  # compute other parameters
  beta.each = lapply(1:p,function(j){beta[(Xdims_cumul[j]+1):Xdims_cumul[j+1],]})
  beta.norm = sapply(1:p,function(j){vector.norm(beta.each[[j]],Ymu,Yspace,'L2')})
  beta.vectors = lapply(1:p,function(j){pca[[j]]$vectors})
  beta.vectors = reduce.dimension(beta.vectors,opt.Xdim.max,margin=2)
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,pca[[j]]$mu,Ymu)})
  proper.indices = which(beta.norm!=0)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['beta']] = beta
  object[['pca']] = pca
  object[['link']] = link
  object[['beta.each']] = beta.each
  object[['beta.norm']] = beta.norm
  object[['beta.vectors']] = beta.vectors
  object[['beta.tensor']] = beta.tensor
  object[['proper.indices']] = proper.indices
  object[['Xdim.max']] = opt.Xdim.max
  object[['runtime']] = runtime
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  class(object) = 'GLM'
  
  return(object)
}



