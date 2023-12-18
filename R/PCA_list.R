# Functions for PCA.manifold.list


#' @title Principal component analysis for manifold-valued data,
#' 
#' @description 
#' PCA fitting for a list of manifold-valued data.
#' 
#' @param Xdata a list of manifold-valued data.
#' \describe{
#'       \item{j}{A \eqn{p} list of generated data. Each jth element is an \eqn{n}-by-\eqn{dim_j} matrix.}
#'       \item{Xspaces}{A \eqn{p} vector of underlying spaces of \eqn{X_j}.}
#'       \item{p}{A number of \eqn{X_j}.}
#' }
#' 
#' @return a \eqn{p} list of PCA.manifold object, see \code{\link{PCA.manifold}}.
#' @export
PCA.manifold.list = function(Xdata){
  pca.list = PCA_list(Xdata)
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
#' @param object a object obtained from \code{\link{PCA.manifold.list}}.
#' @param Xdatanew a \eqn{p} list of new manifold-valued data.
#' 
#' @return a \eqn{p} list of score matrices, see \code{\link{predict.PCA.manifold}}.
#' @export
predict.PCA.manifold.list = function(object,Xdatanew){
  scores = predict_PCA_list(object,Xdatanew)
  return(scores)
}
# predict.PCA.manifold.list2 = function(object,Xdatanew){
#   p = object[['p']]
#   scores = lapply(1:p,function(j){predict(object[[j]],Xdatanew[[j]])})
#   return(scores)
# }


