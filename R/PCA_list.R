


#' @title Principal Component Analysis for a List of Manifold-valued Data,
#' 
#' @description 
#' Performs the principal component analysis (spectral decomposition) for a list of manifold-valued covariates.
#' 
#' @param Xdata a list of manifold-valued covariates with the following arguments:
#' \describe{
#'       \item{j}{a \eqn{p} list of manifold-valued covariates, where each \eqn{j}th element is an \eqn{n\times T_j} matrix.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
#' @param alpha a truncation parameter of the number of basis vectors, used only for infinite dimensional \eqn{\mathcal{M}}. Selects the first index where the cumulative variance is equal to or greater than \eqn{\alpha} of the total variance.
#' 
#' @return a 'PCA.manifold.list' object with the following arguments:
#' \describe{
#'       \item{j}{a \code{\link{PCA.manifold}} object for \eqn{X_j}.}
#'       \item{spaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.}
#'       \item{p}{the number of \eqn{X_j}.}
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


#' @title Prediction of Score Matrices for a List of Manifold-valued Data
#' 
#' @description
#' Computes a list of score matrices for a list of the manifold-valued data based on \code{\link{PCA.manifold.list}}.
#' 
#' @param object a \code{\link{PCA.manifold.list}} object.
#' @param Xnew a new list of manifold-valued covariates, see the "Xdata" argument in \code{\link{PCA.manifold.list}}.
#' 
#' @return a \eqn{p} list of score matrices.
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


