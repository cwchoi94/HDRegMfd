### Cross validation for LR.oracle



# Compute loss for each Xdim.max
get.loss.LR.oracle = function(X,LogY,Xnew,LogYnew,Ymu,Xdim.max,proper.indices){
  X = reduce.dimension(X,Xdim.max)
  Xnew = reduce.dimension(Xnew,Xdim.max)
  Xnew = do.call(cbind,Xnew)
  
  # compute oracle LSE
  beta = compute_beta(X,LogY,proper.indices)
  
  LogYhat = Xnew %*% beta
  rmse = vector.norm(LogYnew-LogYhat,Ymu,Yspace,'L2')/sqrt(n2)
  
  return(rmse)
}


#' GCV for an LR.oracle function
#' 
#' @param Xorg a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yorg an \eqn{n\times m} response matrix.
#' @param Xnew a new list of manifold-valued covariates.
#' @param Ynew a new \eqn{n'\times m} response matrix.
#' @param Yspace an underlying space of Yorg and Ynew.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#'
#' @return an \eqn{\code{LR}} object
#'    \describe{
#'       \item{Xdim.max.list}{a vector of Xdim.max.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time.}
#'       \item{...}{see \code{\link{LR.oracle}}.}
#' }
#' @export
LR.oracle.GCV = function(Xorg,Yorg,Xorgnew,Yorgnew,Yspace,Xdim.max.list,proper.indices){
  
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
  
  loss.list = sapply(1:r2,function(i){get.loss.LR.oracle(X,LogY,Xnew,LogYnew,Ymu,Xdim.max.list[i],proper.indices)})
  
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
  beta.tensor = lapply(1:p,function(j){make.tensor(beta.vectors[[j]],beta.each[[j]],pca$spaces[j],Yspace,Ymu,pca[[j]]$mu)})
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
  class(object) = 'LR'
  
  return(object)
}




#' Kfold cross validation for an LR.oracle function
#' 
#' @param Xall a list of manifold-valued covariates, see \code{\link{PCA.manifold.list}}.
#' @param Yall an \eqn{n\times m} response matrix.
#' @param Yspace an underlying space of Yorg and Ynew.
#' @param kfold a number of kfold CV, int>0.
#' @param Xdim.max.list a vector of max dimension of \eqn{X_j}.
#' @param proper.indices a vector of indices of relevant \eqn{X_j}.
#' @param seed a random seed, int>0, default: non-random (NULL).
#'
#' @return an \eqn{\code{LR}} object
#'    \describe{
#'       \item{Xdim.max.list}{a vector of Xdim.max.}
#'       \item{loss.list}{a list of loss for each CV step.}
#'       \item{runtime}{running time.}
#'       \item{...}{see \code{\link{LR.oracle}}.}
#' }
#' @export
LR.oracle.kfold = function(Xall,Yall,Yspace,kfold,Xdim.max.list,proper.indices,seed=NULL){
  
  start.time = Sys.time()
  
  # check validility of inputs
  Check.manifold(Yspace)
  
  # define basic parameters
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
    
    loss = sapply(1:r2,function(i){get.loss.LR.oracle(Xorg,LogY,Xnew,LogYnew,Ymu,Xdim.max.list[i],proper.indices)})
    loss.list[,idx] = unlist(loss)
  }
    
  loss = rowMeans(loss.list)
  opt.Xdim.max = Xdim.max.list[which.min(loss)]
  
  # fit opt model
  object = LR.oracle(Xall,Yall,Yspace,opt.Xdim.max,proper.indices)
  
  runtime = hms::hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  
  object[['Xdim.max.list']] = Xdim.max.list
  object[['loss.list']] = loss.list
  object[['runtime']] = runtime
  
  return(object)
}




