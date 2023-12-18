### functions to generate Riemannian metric space-valued data.


#####################################################
#####################################################
### Frechet means generate


#' Generate Frechet means of \eqn{X_j}
#' 
#' @param m length of each observation.
#' @param space a name of space, see \code{\link{Check.manifold}}.
#' 
#' @return a \eqn{m} vector.
Xmu.generate = function(m,space){
  Check.manifold(space)
  if (space=='Euclid'){
    mu = rep(0,m)
  } else if (space=='functional'){
    m = 100
    mu = rep(0,m)
  } else if (space=='simplex'){
    mu = inv.clr.simplex(rep(0,m))[1,]
  } else if (space=='sphere'){
    mu = c(rep(0,m-1),1)
  } else if (space=='BayesHilbert'){
    m = 100
    mu = inv.clr.BayesHilbert(rep(0,m))[1,]
  } else if (space=='Wasserstein'){
    m = 100
    t = seq(0,1,length.out=m+2)[2:(m+1)]
    mu = qnorm(t)
  }
  return(mu)
}


#' Generate the Frechet means for \eqn{Y}
#' 
#' @param m length of each observation.
#' @param space a name of space, see \code{\link{Check.manifold}}.
#' 
#' @return an \eqn{m} vector.
Ymu.generate = function(m,space){
  Check.manifold(space)
  if (space=='Euclid'){
    mu = rep(0,m)
  } else if (space=='functional'){
    m = 100
    mu = rep(0,m)
  } else if (space=='simplex'){
    mu = inv.clr.simplex(rep(0,m))[1,]
  } else if (space=='sphere'){
    mu = c(rep(0,m-1),1)
  } else if (space=='BayesHilbert'){
    m = 100
    mu = inv.clr.BayesHilbert(rep(0,m))[1,]
  } else if (space=='Wasserstein'){
    m = 100
    t = seq(0,1,length.out=m+2)[2:(m+1)]
    mu = qnorm(t)
  }
  return(mu)
}



#####################################################
#####################################################
### Covariate generate


#' Generate real covariates
#' 
#' @param n a number of data.
#' @param p a number of covariates.
#' @param Zrho a correlation parameter which ranges from -1 to 1.
#' @param Zsigma a variance parameter, default 1.
#' 
#' @return A \eqn{n}-by-\eqn{p} covariate matrix.
covariates.generate.real = function(n,p,Zrho,Zsigma){
  Zcov = outer(1:p,1:p,function(i,j){Zrho^(abs(i-j))})*Zsigma
  Z = MASS::mvrnorm(n,rep(0,p),Zcov)
  return(Z)
}


#' Generate each manifold-valued covariates
#' 
#' @param Xi an \eqn{n}-by\eqn{m'} score matrix.
#' @param dim a dimension of \eqn{X}.
#' @param space the underlying spaces, see \code{\link{Check.manifold}}
#' 
#' @return an \eqn{n}-by-\eqn{m} matrix with dimension. Each row is an element of the manifold.
#' @export
covariates.generate.each = function(Xi,dim,space='Euclid'){
  Check.manifold(space)
  mu = Xmu.generate(dim,space)
  basis = basis.manifold(mu,ncol(Xi),space)
  LogX = Xi %*% basis
  X = RieExp.manifold(mu,LogX,space)
  return(X)
}


#' @title Generate a list of Riemannian metric space-valued covariates.
#' 
#' @description
#' Generate the covariate list.
#' If the underlying space of \eqn{X_j} is finite-dimensional, the dimension of \eqn{X_j} is equal to \eqn{dim_j}.
#' If the underlying space of \eqn{X_j} is inifite-dimensional, the dimension of \eqn{X_j} is set 100.
#' 
#' @param n a number of data
#' @param Xspaces a \eqn{p} vector of underlying spaces of \eqn{X_j}.
#' @param dims a \eqn{p} vector of dimension of \eqn{X_j}.
#' @param Xrho a correlation parameter which ranges from -1 to 1.
#' @param Xsigma a correlation parameter which ranges from -1 to 1.
#' 
#' @return a list of generated data
#'    \describe{
#'       \item{j}{A \eqn{p} list of generated data. Each jth element is an \eqn{n}-by-\eqn{dim_j} matrix.}
#'       \item{Xspaces}{A \eqn{p} vector of underlying spaces of \eqn{X_j}.}
#'       \item{p}{A number of \eqn{X_j}.}
#' }
covariates.generate = function(n,Xspaces,dims,Xrho=0.5,Xsigma=1){
  # compute intrinsic dimension
  p = length(dims)
  dims_ = sapply(1:p,function(i){if(Xspaces[i] %in% c('simplex','sphere')){dims[i]-1}else{dims[i]}})
  
  # generate scores
  Z = covariates.generate.real(n,sum(dims_),Xrho,Xsigma)
  
  # weight bound for w_jk
  w.ftn = function(l){
    if (l<=4){
      w = l^(-1)
    } else if (l>4){
      w = l^(-1)
    }
    return(w)
  }
  
  # generate X_j
  Xdata = list()
  dims.cumul = c(0,cumsum(dims_))
  for (j in 1:p){
    a = dims.cumul[j]+1
    b = dims.cumul[j+1]
    
    Xi = as.matrix(Z[,a:b])
    Xi = Xi %*% diag(sapply(1:ncol(Xi),function(l){w.ftn(l)}))
    Xdata[[j]] = covariates.generate.each(Xi,dims[j],Xspaces[j])
  }
  Xdata[['spaces']] = Xspaces
  Xdata[['p']] = p
  return(Xdata)
}



######################################################
######################################################
### Operator generate


#' @title Generate operators from \eqn{H_1} to \eqn{H_2}
#' 
#' @description
#' Generate an operator \eqn{B_j} from \eqn{H_1} to \eqn{H_2}. 
#' It can be identified as an element in the tensor product space \eqn{H_1\otimes H_2}. 
#' 
#' @param Xspace an underlying space of \eqn{H_1}.
#' @param Yspace an underlying space of \eqn{H_2}.
#' @param Xdim an intrinsic dimension of \eqn{H_1}.
#' @param Ydim an intrinsic dimension of \eqn{H_2}.
#' @param beta.norm the tensor norm.
#' 
#' @return a tensor product element in \eqn{H_1\otimes H_2} generated by \code{\link{make.tensor}}.
tensor.beta.generate.each = function(Xspace,Yspace,Xdim,Ydim,beta.norm=1){
  Xmu = Xmu.generate(Xdim,Xspace)
  Ymu = Ymu.generate(Ydim,Yspace)
  
  Xbasis = basis.manifold(Xmu,Xdim,Xspace)
  Ybasis = basis.manifold(Ymu,Ydim,Yspace)
  
  # weight bound for b_l1l2
  b.ftn = function(l1,l2){
    if (l1<=5 & l2<=5){
      b = 1
    } else if (l1<=5 & l2>5){
      b = (l2-2)^(-1)
    } else if (l1>5 & l2<=5){
      b = (l1-2)^(-2)
    } else if (l1>5 & l2>5){
      b = (l1-2)^(-2)*(l2-2)^(-1)
    }
    return(b)
  }
  
  b = sapply(1:nrow(Xbasis),function(l1){sapply(1:nrow(Ybasis),function(l2){runif(1,-1,1)*b.ftn(l1,l2)})})
  b = t(b) %*% Ybasis
  beta = make.tensor(Xbasis,b,Xspace,Yspace,Xmu,Ymu)
  beta$element2 = beta$element2 / norm.tensor(beta) * beta.norm
  return(beta)
}


#' @title Generate a list of operators from \eqn{H_j} to \eqn{H_Y}.
#' 
#' @description
#' Generate a list of operators \eqn{B_j} from \eqn{H_j} to \eqn{H_Y}.
#' Each operator can be identified as an element in the tensor product space \eqn{H_j\otimes H_Y}.
#' 
#' @param Xspaces a \eqn{p} vector of underlying spaces of \eqn{X_j}.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param Xdims a \eqn{p} vector of intrinsic dimensions of \eqn{H_j}.
#' @param Ydim an intrinsic dimension of \eqn{H_Y}.
#' @param proper.indices an indices of nonzero operators.
#' @param beta.norm the tensor norm.
#' 
#' @return a list of operators,
#'    \describe{
#'       \item{j}{a \eqn{p} list of generated data. Each jth element is an \eqn{n}-by-\eqn{Xdim_j} matrix.}
#'       \item{Xspaces}{a \eqn{p} vector of spaces of \eqn{X_j}.}
#'       \item{p}{a number of \eqn{X_j}.}
#' } 
tensor.beta.generate = function(Xspaces,Yspace,Xdims,Ydim,proper.indices=NULL,beta.norm=1,...){
  p = length(Xspaces)
  if(is.null(proper.indices)){proper.indices = seq(p)}
  
  beta.list = list()
  for (j in 1:p){
    beta.list[[j]] = tensor.beta.generate.each(Xspaces[[j]],Yspace,Xdims[j],Ydim,beta.norm*(j%in%proper.indices))
  }
  
  beta.list[['Xspaces']] = Xspaces
  beta.list[['p']] = p
  return(beta.list)
}




############################################################
############################################################
### error generate


#' Gerenate random errors on the tangent space at mu.
#' 
#' @param n a number of data.
#' @param space an underlying space.
#' @param dim an intrinsic dimension of the tangent space.
#' @param error.rho a correlation parameter.
#' @param error.std a standard deviation parameter.
#' 
#' @return an \eqn{n}-by-\eqn{m} matrix of generated random error.
error.generate = function(n,space,dim,error.rho=0.5,error.std=1){
  mu = Ymu.generate(dim,space)
  
  if (space %in% c('Euclid','simplex','sphere')){
    basis = basis.manifold(mu,dim,space)
    dim = nrow(basis)
    
    error.cov = matrix(error.rho,dim,dim) + diag(1-error.rho,dim)
    error.cov = error.cov * error.std^2 / dim
    error = mvrnorm(n,rep(0,dim),error.cov)
    error = error %*% basis
  } else{
    m = length(mu)
    const = 5/(m-1)
    error.cov = outer(1:m,1:m,function(x,y){error.rho^(const*abs(x-y))})
    error.cov = error.cov * error.std^2
    error = mvrnorm(n,rep(0,m),error.cov)
    
    if (space=='BayesHilbert'){
      rho2 = error.rho^(const)
      error.std.correction = sqrt(1-1/(m^2 * (1-rho2))*(m*(1+rho2) - 2*(rho2-rho2^(m+1))/(1-rho2)))
      error = error/error.std.correction
      error = error - rowMeans(error)
    }
  }
  
  return(error)
}


#' @title Generate linear regression data
#' 
#' @description
#' Generate linear regression data.
#' The seed for generating beta is fixed as zero, so each beta is equal in simulations.
#' 
#' @param n a number of observation.
#' @param Xspaces a \eqn{p} vector of underlying spaces of \eqn{X_j}.
#' @param Yspace an underlying space of \eqn{Y}.
#' @param Xdims a \eqn{p} vector of intrinsic dimensions of \eqn{H_j}.
#' @param Ydim an intrinsic dimension of \eqn{H_Y}.
#' @param proper.indices an indices of nonzero operators.
#' @param beta.norm the tensor norm of nonzero operators.
#' @param Xrho a correlation parameter for \eqn{X_j}.
#' @param Xsigma a variance parameter for \eqn{X_j}.
#' @param error.rho a correlation parameter for error.
#' @param error.std a standard deviation parameter for error.
#' @param seed a random seed (1<=seed)
#' 
#' @return a list of generated data
#'    \describe{
#'       \item{X}{a \eqn{p} list of generated covariates. Each \eqn{X_j} is an \eqn{n}-by-\eqn{K_j} matrix.}
#'       \item{Y}{an \eqn{n}-by-\eqn{m} matrix of response.}
#'       \item{beta}{a \eqn{p} list of generated operators.}
#'       \item{error}{a generated random error.}
#'       \item{Xbeta.each}{a \eqn{p} list of \eqn{B_j(X_j)}.}
#'       \item{Xbeta}{a sum of Xbeta.each.}
#'       \item{Ymu}{the Frechet mean of \eqn{Y}.}
#'       \item{...}{the other input parameters.}
#' }
linear.data.generate = function(n,Xspaces,Yspace,Xdims,Ydim,proper.indices=NULL,beta.norm=1,Xrho=0.5,Xsigma=1,
                                error.rho=0.5,error.std=1,seed=1){
  
  # generate same beta for each simulation
  set.seed(0)
  beta = tensor.beta.generate(Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm)
  
  # generate X and error
  set.seed(seed)
  X = covariates.generate(n,Xspaces,Xdims,Xrho,Xsigma)
  error = error.generate(n,Yspace,Ydim,error.rho,error.std)
  
  # compute Xbeta
  Xbeta.each = lapply(1:length(Xdims),function(j){operator.tensor(beta[[j]],X[[j]])})
  Xbeta = Reduce('+',Xbeta.each)
  
  # make Y
  Ymu = Ymu.generate(Ydim,Yspace)
  LogY = Xbeta + error
  ExpXbeta = RieExp.manifold(Ymu,Xbeta,Yspace)
  Y = RieExp.manifold(Ymu,LogY,Yspace)
  
  data = list(X=X,Y=Y,beta=beta,error=error,ExpXbeta=ExpXbeta,Ymu=Ymu,
              LogY=LogY,Xbeta=Xbeta,Xbeta.each=Xbeta.each,
              n=n,p=length(Xdims),Xspaces=Xspaces,Yspace=Yspace,Xdims=Xdims,Ydim=Ydim,
              proper.indices=proper.indices,beta.norm=beta.norm,Xrho=Xrho,Xsigma=Xsigma,
              error.rho=error.rho,error.std=error.std,seed=seed)
  return(data)
}



