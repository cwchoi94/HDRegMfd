# Functions for PCA.manifold.list


#' @title Principal component analysis for manifold-valued data,
#' 
#' @description 
#' PCA fitting for a list of manifold-valued data.
#' 
#' @param Xdata a list of manifold-valued data.
#' \describe{
#'       \item{j}{a \eqn{p} list of generated data. Each \eqn{j}th element is an \eqn{n\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of underlying spaces of \eqn{X_j}.}
#'       \item{p}{a number of \eqn{X_j}.}
#' }
#' @param alpha a truncation parameter of the number of vectors, only used for infinite dimensional manifolds. Select the first index where the sum of the variances is equal to or greater than the alpha of the total.
#' 
#' @return a PCA.manifold.list object.
#' \describe{
#'       \item{j}{a \code{\link{PCA.manifold}} object for \eqn{X_j}.}
#'       \item{spaces}{a \eqn{p} vector of underlying spaces of \eqn{X_j}.}
#'       \item{p}{a number of \eqn{X_j}.}
#' }
#' @export
PCA.manifold.list = function(Xdata,alpha=0.95){
  pca.list = PCA_list(Xdata,alpha)
  class(pca.list) = 'PCA.manifold.list'
  return(pca.list)
}
# PCA.manifold.list2 = function(Xdata){
#   spaces = Xdata[['spaces']]
#   p = Xdata[['p']]
# 
#   pca.list = lapply(1:p,function(j){PCA.manifold(Xdata[[j]],spaces[j])})
#   pca.list[['p']] = p
#   pca.list[['spaces']] = spaces
#   class(pca.list) = 'PCA.manifold.list2'
#   return(pca.list)
# }


#' @title Prediction of principal component score for manifold-valued data
#' 
#' @description 
#' Prediction of principal component scores for a list of manifold-valued data.
#' 
#' @param object a \code{\link{PCA.manifold.list}} object.
#' @param Xdatanew a new \eqn{p} list of manifold-valued data.
#' 
#' @return a \eqn{p} list of score matrices, see \code{\link{predict.PCA.manifold}}.
#' @export
predict.PCA.manifold.list = function(object,Xdatanew){
  scores = predict_PCA_list(object,Xdatanew)
  p = object$p
  scores = lapply(1:p,function(j){vec.to.mat(scores[[j]])})
  return(scores)
}
# predict.PCA.manifold.list2 = function(object,Xdatanew){
#   p = object[['p']]
#   scores = lapply(1:p,function(j){predict(object[[j]],Xdatanew[[j]])})
#   return(scores)
# }


