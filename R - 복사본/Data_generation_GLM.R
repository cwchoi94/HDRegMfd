



#' @title Generate generalized linear regression data
#' 
#' @description
#' Generate generalized linear regression data.
#' The response \eqn{Y} given the covariates \eqn{X} follows an exponential family distribution, and thus \eqn{Y} is real-valued.
#' The seed for generating Hilbert-Schmidt operators \eqn{\mathfrak{B}_j} is fixed as zero, ensuring that each \eqn{\mathfrak{B}_j} is the same across simulations.
#' 
#' @inheritParams LM.data.generate
#' 
#' @param link a link function, see \code{\link{Check.link}}.
#' @param Ydim the intrinsic dimension of \eqn{Y}, typically 1.
#' @param beta0.norm an integer specifying the norm of \eqn{\beta_0^*}, with a default value of 1.
#' @param c.beta a \eqn{p} vector of pre-computed parameters to ensure that \eqn{\mathbb{E}\| \mathfrak{B_j}(Log_{\mu_j}X_j) \|^2 = (\text{beta.norm[j]})^2}, only computed if not provided.
#' 
#' @return a list of data containing:
#'    \describe{
#'       \item{X}{a list of manifold-valued covariates, see \code{\link{covariates.generate}}.}
#'       \item{Y}{an \eqn{n\times m} matrix of responses.}
#'       \item{p}{the number of \eqn{X_j}.}
#'       \item{link}{a link function, see \code{\link{Check.link}}.}
#'       \item{theta}{an \eqn{n\times m} matrix of canonical parameters for \eqn{Y}.}
#'       \item{Ymu}{an \eqn{n\times m} matrix of the means of \eqn{Y_i}, which is an inverse link of theta.}
#'       \item{beta}{a \eqn{p} list of Hilbert-Schmidt operators, see \code{\link{tensor.beta.generate}}.}
#'       \item{beta0}{a \eqn{m} vector, with \eqn{\|\beta_0^*\| = \text{beta0.norm}}.}
#'       \item{Xmu}{a \eqn{p} list of the Frechet means \eqn{\mu_j} of \eqn{X_j}.}
#'       \item{LogX}{a \eqn{p} list of \eqn{Log_{\mu_j}X_j}, the Riemannian logarithmic transformations of \eqn{X_j}.}
#'       \item{Xbeta.each}{a \eqn{p} list of \eqn{\mathfrak{B}_j(Log_{\mu_j}X_j)}.}
#'       \item{Xbeta}{the sum of Xbeta.each, i.e., \eqn{\sum_{j=1}^p\mathfrak{B}_j(Log_{\mu_j}X_j)}.}
#'       \item{...}{other input parameters.}
#' }
#' @export
GLM.data.generate = function(n,Xspaces,Xdims,link='binomial',Ydim=1,proper.indices=NULL,beta.norm=1,beta0.norm=1,Xrho=0.5,Xsigma=1,ngrid=100,seed=1,c.beta=NULL){
  
  Yspace = 'Euclid'
  Check.link(link)
  if (link!='multinomial'){Ydim = 1}
  if(is.null(proper.indices)){proper.indices = seq(length(Xdims))}
  if(length(beta.norm)==1){beta.norm = rep(beta.norm,length(Xdims))}
  p = length(Xdims)
  s = length(proper.indices)
  
  # generate Xmu 
  Xmu.list = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j],ngrid)})
  Ymu = Ymu.generate(Ydim,Yspace,ngrid) # pseudo Ymu
  
  # generate same beta for each simulation
  set.seed(0)
  beta = tensor.beta.generate(Xspaces,Yspace,Xmu.list,Ymu,Xdims,Ydim,proper.indices,1)
  beta.oracle = lapply(proper.indices,function(j){beta[[j]]})
  
  # generate nuisance X to set |B_j(LogX_j)|=beta.norm[j]
  if (is.null(c.beta)){
    X.base = covariates.generate(10000,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)$X
    LogX.base = lapply(proper.indices,function(j){RieLog.manifold(Xmu.list[[j]],X.base[[j]],Xspaces[j])})
    Xbeta.each.base = lapply(1:s,function(j){operator.tensor(beta.oracle[[j]],LogX.base[[j]])})
    c.beta.tmp = sapply(1:s,function(j){sqrt(mean(norm.manifold(Xbeta.each.base[[j]])^2))})
    c.beta = rep(0,p)
    c.beta[proper.indices] = c.beta.tmp
  }
  
  for (j in proper.indices){
    beta[[j]]$element2 = beta[[j]]$element2 / c.beta[j] * beta.norm[j]
  }
  
  # generate X
  set.seed(seed)
  X = covariates.generate(n,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)$X
  
  # compute Xbeta
  LogX = lapply(1:length(Xdims),function(j){RieLog.manifold(Xmu.list[[j]],X[[j]],Xspaces[j])})
  Xbeta.each = lapply(1:length(Xdims),function(j){operator.tensor(beta[[j]],LogX[[j]])})
  Xbeta = Reduce('+',Xbeta.each)
  
  # make Y
  beta0 = rep(1,Ydim)/Ydim * beta0.norm
  theta = beta0 + Xbeta
  Ymu = Inv_Link(theta,link)
  if (link=='binomial'){
    Y = matrix(rbinom(n,1,prob=Ymu),nrow=n,ncol=1)
  }else if (link=='poisson'){
    Y = matrix(rpois(n,Ymu),nrow=n,ncol=1)
  }
  
  data = list(X=X,Y=Y,link=link,theta=theta,Ymu=Ymu,beta=beta,beta0=beta0,c.beta=c.beta,
              Xmu.list=Xmu.list,LogX=LogX,Xbeta=Xbeta,Xbeta.each=Xbeta.each,
              n=n,p=length(Xdims),Xspaces=Xspaces,Xdims=Xdims,Ydim=Ydim,
              proper.indices=proper.indices,beta.norm=beta.norm,beta0.norm=beta0.norm,
              Xrho=Xrho,Xsigma=Xsigma,seed=seed)
  return(data)
}


