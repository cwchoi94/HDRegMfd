


# Check a transform function.
Check.transform = function(transform){
  if (!(transform %in% c('linear','Gaussian','trigonometric'))){
    stop("The transformation must be one of 'linear','Gaussian','trigonometric'")
  }
}



#' @title Transform Score Matrices
#' 
#' @description
#' Transforms score matrices into the interval \eqn{[0,1]}.
#' This function is implemented for three types of transformations:
#' - linear: Min-max scaling.
#' - Gaussian: Cumulative distribution function of the standard Gaussian distribution.
#' - trigonometric: Inverse tangent function.
#' If \code{normalize} is \code{TRUE}, the function standardizes each column of the score matrices.
#' 
#' @param score.list a \eqn{p} list of score matrices. It is only used for computing \code{min.list}, \code{max.list} and \code{omega.sq.list}.
#' @param transform a method of transformation. Must be one of 'linear', 'Gaussian' or 'trigonometric' (default: 'linear').
#' @param normalize a parameter indicating whether to standardize each column of score matrices (default: \code{TRUE}).
#' 
#' @return a \code{Transform.Score} object with the following components:
#' \describe{
#'       \item{p}{the length of list of score matrices.}
#'       \item{index.mat}{a \eqn{L_+^* \times 3} matrix, with each row being \eqn{(i,j,k)} such that \eqn{i\equiv (j,k)}.}
#'       \item{min.list}{a \eqn{p} list of minimum values of columns of (normalized) score matrices.}
#'       \item{max.list}{a \eqn{p} list of maximum values of columns of (normalized) score matrices.}
#'       \item{transform}{the method of transformation.}
#'       \item{normalize}{the parameter indicating whether standardize each column of score matrices.}
#'       \item{omega.sq.list}{a \eqn{p} list of standard deviations of each column of score matrices.}
#' }
#' @export
Transform.Score = function(score.list,transform='linear',normalize=TRUE){
  Check.transform(transform)
  
  p = length(score.list)
  dims = sapply(score.list,ncol)
  
  # make an index matrix identifying (j,k) = l
  index.mat = do.call(rbind,lapply(1:p,function(j){cbind(j,1:dims[j])}))
  index.mat = cbind(1:sum(dims),index.mat)
  
  # if normalize==TRUE, we normalize each column of score matrices
  omega.sq.list = lapply(1:p,function(j){sqrt(colMeans(score.list[[j]]^2))})
  if (normalize==TRUE){
    score.list = lapply(1:p,function(j){sweep(score.list[[j]],2,omega.sq.list[[j]],FUN='/')})
  }
  
  min.list = lapply(1:p,function(j){apply(score.list[[j]],2,min)})
  max.list = lapply(1:p,function(j){apply(score.list[[j]],2,max)})
  
  object = list()
  object[['p']] = p
  object[['index.mat']] = index.mat
  object[['min.list']] = min.list
  object[['max.list']] = max.list
  object[['transform']] = transform
  object[['normalize']] = normalize
  object[['omega.sq.list']] = omega.sq.list
  class(object) = 'Transform.Score'
  
  return(object)
}



#' @title Prediction of Transform Score Matrices
#' 
#' @description
#' Transforms given score matrices into the interval \eqn{[0,1]} using a \code{\link{Transform.Score}} object.
#' 
#' @param a \code{\link{Transform.Score}} object.
#' @param score.list a \eqn{p} list of score matrices to be transformed.
#' 
#' @return a \eqn{p} list of transformed score matrices.
#' @export
predict.Transform.Score = function(object,score.list){
  p = object[['p']]
  min.list = object[['min.list']]
  max.list = object[['max.list']]
  
  # normalize
  omega.sq.list = object[['omega.sq.list']]
  if (object[['normalize']]==TRUE){
    score.list = lapply(1:p,function(j){sweep(score.list[[j]],2,omega.sq.list[[j]],FUN='/')})
  }
  
  # transform
  transform = object[['transform']]
  if (transform=='linear'){
    score.transform.list = lapply(1:p,function(j){
      min.vec = min.list[[j]]
      range.vec = max.list[[j]] - min.vec
      range.vec[range.vec==0] = 1 # if min==max, it may cause division error
      score = score.list[[j]]
      
      score.transform = sweep(sweep(score,2,min.vec,'-'),2,range.vec,'/')
      return(score.transform)
    })
  }else if (transform=='Gaussian'){
    score.transform.list = lapply(1:p,function(j){
      return(pnorm(score.list[[j]]))
    })
  }else if (transform=='trigonometric'){
    score.transform.list = lapply(1:p,function(j){
      return(atan(score.list[[j]])/pi + 0.5)
    })
  }
  
  score.transform = do.call(cbind,score.transform.list)
  
  return(score.transform)
}






