### data generation codes for simulations.



### basic functions


#' @title Generate a non-symmetric identity matrix
#' 
#' @description
#' Creates a \eqn{p1 \times p2} matrix \eqn{I}, where each element \eqn{I_{jk}} is equal to \eqn{\mathbb{I}(j=k)}.
#' 
#' @param p1 the number of rows in the matrix.
#' @param p2 the number of columns in the matrix.
#' 
#' @return a \eqn{p1 \times p2} matrix where each element \eqn{I_{jk}} is \eqn{\mathbb{I}(j=k)}.
make.nonsym.Idmat = function(p1,p2){
  if (p1>=p2){
    mat = rbind(diag(p2),matrix(0,p1-p2,p2))
  } else{
    mat = cbind(diag(p1),matrix(0,p1,p2-p1))
  }
  return(mat)
}



Transform.Xi = function(Xi,space,a){
  Xi = Xi %*% diag(sapply(1:ncol(Xi),function(k){w.ftn(k,a)}))
  return(Xi)
}

# weight bound for w_jk
w.ftn = function(k,a=1){
  if (k<=4){
    w = k^(-a)
  } else if (k>4){
    w = k^(-a)
  }
  return(w)
}


# weight bound for b_l1l2
b.ftn = function(l1,l2){
  if (l1<=5 & l2<=5){
    z = 1
  } else if (l1<=5 & l2>5){
    z = (l2-2)^(-1)
  } else if (l1>5 & l2<=5){
    z = (l1-2)^(-2)
  } else if (l1>5 & l2>5){
    z = (l1-2)^(-2)*(l2-2)^(-1)
  }
  return(z)
}



#####################################################
#####################################################
### Generate the Frechet means 


#' @title Generate Frechet means of \eqn{X_j}
#' 
#' @description 
#' Generate Frechet means \eqn{\mu_j} of \eqn{X_j}.
#' 
#' @param m the length of each observation. Typically, this equals the dimension, but it can differ in some cases.
#' @param space the name of space, see \code{\link{Check.manifold}}.
#' @param ngrid the number of grids. This is only used for infinite-dimensional spaces, see \code{\link{Check.manifold}}: "functional", "BayesHilbert" and "Wasserstein" spaces.
#' 
#' @return the Frechet mean \eqn{\mu_j} of \eqn{X_j}, an \eqn{m} vector.
Xmu.generate = function(m,space,ngrid=100){
  Check.manifold(space)
  if (space=='Euclid'){
    mu = rep(0,m)
  } else if (space=='simplex'){
    mu = inv.clr.simplex(rep(0,m))[1,]
  } else if (space=='SPD.LogEuclid'){
    m = sqrt(m)
    mu = as.vector(0.8*diag(m)+0.2)
  } else if (space=='SPD.AffInv'){
    m = sqrt(m)
    mu = as.vector(0.8*diag(m)+0.2)
  } else if (space=='sphere'){
    mu = c(rep(0,m-1),1)
  } else if (space=='functional'){
    mu = rep(0,ngrid)
  } else if (space=='BayesHilbert'){
    mu = inv.clr.BayesHilbert(rep(0,ngrid))[1,]
  } else if (space=='Wasserstein'){
    t = seq(0,1,length.out=ngrid+2)[2:(ngrid+1)]
    mu = 10*t
  }
  return(mu)
}


#' @title Generate the Frechet mean for \eqn{Y}
#' 
#' @description
#' This code is the same as \code{\link{Xmu.generate}}.
#' 
#' @inheritParams Xmu.generate
#' 
#' @return The Frechet mean \eqn{\mu_Y} of \eqn{Y}, an \eqn{m} vector.
Ymu.generate = function(m,space,ngrid=100){
  mu = Xmu.generate(m,space,ngrid)
  return(mu)
}



#####################################################
#####################################################
### Generate covariates


#' @title Generate real covariates
#' 
#' @description
#' Generate an \eqn{n\times p} matrix \eqn{\bm{X}=(X_{ij})}, where each \eqn{X_{ij}} represents a covariate
#' 
#' @param n the number of data points.
#' @param p the number of covariates.
#' @param Zrho a correlation parameter ranging from -1 to 1.
#' @param Zsigma a variance parameter, with a default value of 1.
#' 
#' @return an \eqn{n\times p} covariate matrix.
covariates.generate.real = function(n,p,Zrho,Zsigma){
  Zcov = outer(1:p,1:p,function(i,j){Zrho^(abs(i-j))})*Zsigma
  Z = MASS::mvrnorm(n,rep(0,p),Zcov)
  return(Z)
}



#' @title Generate a manifold-valued covariate
#' 
#' @description
#' Given \eqn{j}, generate manifold-valued covariates \eqn{X_{ij}} for \eqn{1\le i\le n}.
#' 
#' @param Xi an \eqn{n\times m'} score matrix, where \eqn{m'} is the intrinsic dimension (the number of scores).
#' @param mu the Frechet mean of \eqn{X_j}.
#' @param space the underlying space, see \code{\link{Check.manifold}}.
#' 
#' @return an \eqn{n\times m} matrix. Each row represents a manifold-valued covariate.
#' @export
covariates.generate.each = function(Xi,mu,space='Euclid'){
  Check.manifold(space)
  basis = basis.manifold(mu,ncol(Xi),space)
  
  LogX = Xi %*% basis
  X = RieExp.manifold(mu,LogX,space)
  return(X)
}




#' @title Generate a list of manifold-valued covariates.
#' 
#' @description
#' Generate a list of manifold-valued covariates.
#' 
#' @param n the number of data points.
#' @param Xspaces a \eqn{p} vector specifying the underlying spaces of \eqn{X_j}.
#' @param Xmu.list a \eqn{p} list of the Frechet means of \eqn{X_j}.
#' @param Xdims a \eqn{p} vector specifying the dimensions of \eqn{X_j}.
#' @param Xrho a correlation parameter ranging from -1 to 1, with a default value of 0.5.
#' @param Xsigma a common standard deviation parameter, with a default value of 1.
#' @param a a parameter such that \eqn{\text{Var}(\xi_{jk})\asymp k^{-2a}}.
#' 
#' @return a list of data containing:
#'    \describe{
#'       \item{j}{a \eqn{p}-list of generated data, where each \eqn{j}th element is an \eqn{n\times D_j} matrix.}
#'       \item{Xspaces}{a \eqn{p} vector specifying the underlying spaces of \eqn{X_j}, see \code{\link{Check.manifold}}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
covariates.generate = function(n,Xspaces,Xmu.list,Xdims,Xrho=0.5,Xsigma=1,a=1){
  # compute intrinsic dimension
  p = length(Xdims)
  Xdims_ = sapply(1:p,function(i){
    if(Xspaces[i] %in% c('simplex','sphere')){
      Xdims[i]-1
    }else if (Xspaces[i]=='SPD.LogEuclid'){
      (Xdims[i]+sqrt(Xdims[i]))/2
    }else{
      Xdims[i]
    }})
  
  # generate scores
  Zeta = lapply(1:p,function(j){covariates.generate.real(n,Xdims_[j],0,Xsigma)})
  Xi = lapply(1:p,function(j){
    if (j==1){
      I2 = make.nonsym.Idmat(Xdims_[j+1],Xdims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j+1]] %*% I2)
    } else if (j==p){
      I1 = make.nonsym.Idmat(Xdims_[j-1],Xdims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j-1]] %*% I1)
    } else{
      I1 = make.nonsym.Idmat(Xdims_[j-1],Xdims_[j])
      I2 = make.nonsym.Idmat(Xdims_[j+1],Xdims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j-1]] %*% I1 + Zeta[[j+1]] %*% I2)
    }
    if (Xspaces[j]=='Wasserstein'){
      Xi.each = (2*pnorm(Xi.each) - 1)/sqrt(2)
    }
    return(Xi.each)
  })
  
  Xi.scaled = lapply(1:p,function(j){Transform.Xi(Xi[[j]],Xspaces[j],a)})
  
  # generate X_j
  X = lapply(1:p,function(j){covariates.generate.each(Xi.scaled[[j]],Xmu.list[[j]],Xspaces[j])})
  X[['spaces']] = Xspaces
  X[['p']] = p
  
  Xdata = list(X=X,Xi=Xi,Xi.scaled=Xi.scaled)
  return(Xdata)
}


######################################################
######################################################
### Generate Hilbert-Schmidt operators


#' @title Generate a Hilbert-Schmidt operator between tangent spaces of \eqn{\mathcal{M}_X} and \eqn{\mathcal{M}_Y}
#' 
#' @description
#' Generate a Hilbert-Schmidt operator \eqn{\mathfrak{B}: T_{\mu_X}\mathcal{M}_X \to T_{\mu_Y}\mathcal{M}_Y}, where \eqn{\mu_X} and \eqn{\mu_Y} are the Frechet means of \eqn{X} and \eqn{Y}, respectively.
#' This operator can be identified as an element in the tensor product space \eqn{T_{\mu_X}\mathcal{M}_X \otimes T_{\mu_Y}\mathcal{M}_Y}. 
#' 
#' @param Xspace the name of the underlying space \eqn{\mathcal{M}_X} of \eqn{X}.
#' @param Yspace the name of the underlying space \eqn{\mathcal{M}_Y} of \eqn{Y}.
#' @param Xmu the Frechet mean of \eqn{X_j}.
#' @param Ymu the Frechet mean of \eqn{Y}.
#' @param Xdim the intrinsic dimension of \eqn{\mathcal{M}_X}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param beta.norm the tensor norm, with a default value of 1.
#' 
#' @return a Hilbert-Schmidt operators (an tensor element) in \eqn{T_{\mu_X}\mathcal{M}_X \otimes T_{\mu_Y}\mathcal{M}_Y} generated by \code{\link{make.tensor}}.
tensor.beta.generate.each = function(Xspace,Yspace,Xmu,Ymu,Xdim,Ydim,beta.norm=1){
  
  Xbasis = basis.manifold(Xmu,Xdim,Xspace)
  Ybasis = basis.manifold(Ymu,Ydim,Yspace)
  
  b = sapply(1:nrow(Xbasis),function(l1){sapply(1:nrow(Ybasis),function(l2){runif(1,-1,1)*b.ftn(l1,l2)})})
  b = vec.to.mat(b)
  b = t(b) %*% Ybasis
  beta = make.tensor(Xbasis,b,Xspace,Yspace,Xmu,Ymu)
  beta$element2 = beta$element2 / norm.tensor(beta) * beta.norm
  return(beta)
}


#' @title Generate a list of Hilbert-SChmidt operators
#' 
#' @description
#' Generate a list of Hilbert-Schmidt operators \eqn{\mathfrak{B}_j: T_{\mu_j}\mathcal{M}_j \to T_{\mu_Y}\mathcal{M}_Y}, where \eqn{\mu_j} and \eqn{\mu_Y} are the Frechet means of \eqn{X_j} and \eqn{Y}, respectively.
#' Each operator is created by \code{\link{tensor.beta.generate.each}} and can be identified as an element in \eqn{T_{\mu_j}\mathcal{M}_j \otimes T_{\mu_Y}\mathcal{M}_Y}.
#' 
#' @param Xspaces a \eqn{p} vector of the names of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.
#' @param Yspace the name of the underlying space \eqn{\mathcal{M}_Y} of \eqn{Y}.
#' @param Xmu.list a \eqn{p} list of the Frechet means of \eqn{X_j}.
#' @param Ymu the Frechet mean of \eqn{Y}.
#' @param Xdims a \eqn{p} vector of the intrinsic dimension of \eqn{\mathcal{M}_j}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param proper.indices an index set \eqn{\mathcal{S}=\{1\le j\le p : \mathfrak{B}_j\neq0\}}.
#' @param beta.norm a \eqn{p} vector (or a single integer) of the tensor norms, used only for \eqn{j\notin\mathcal{S}}. If an integer is provided, the same value is applied for all \eqn{X_j}.
#' 
#' @return a list of Hilbert-Schmidt operators (tensor elements) containing:
#'    \describe{
#'       \item{j}{For \eqn{1\le j\le p}, the \eqn{j}th element is an \eqn{n\times D_j} matrix, with each \eqn{i}th row of \eqn{X_{ij}}.}
#'       \item{Xspaces}{a \eqn{p} vector of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.}
#'       \item{p}{the number of \eqn{X_j}.}
#'       \item{beta.norm}{a \eqn{p} vector of the tensor norms}
#' } 
tensor.beta.generate = function(Xspaces,Yspace,Xmu.list,Ymu,Xdims,Ydim,proper.indices=NULL,beta.norm=1){
  p = length(Xspaces)
  if(is.null(proper.indices)){proper.indices = seq(p)}
  if(length(beta.norm)==1){beta.norm = rep(beta.norm,p)}
  beta.norm[-proper.indices] = 0
  beta.norm = beta.norm[1:p]
  
  beta.list = list()
  for (j in 1:p){
    beta.list[[j]] = tensor.beta.generate.each(Xspaces[[j]],Yspace,Xmu.list[[j]],Ymu,Xdims[j],Ydim,beta.norm[j])
  }
  
  beta.list[['Xspaces']] = Xspaces
  beta.list[['p']] = p
  beta.list[['beta.norm']] = beta.norm
  return(beta.list)
}




############################################################
############################################################
### Generate random errors


#' @title Generate a random error \eqn{\varepsilon} of \eqn{Y}
#' 
#' @description
#' For a random element \eqn{Y} taking values in a manifold \eqn{\mathcal{M}_Y} with the Frechet mean \eqn{\mu_Y}, 
#' generate a random error \eqn{\varepsilon} taking values in \eqn{T_{\mu_Y}\mathcal{M}_Y}.
#' 
#' @param n the number of data points.
#' @param Yspace the underlying space of \eqn{\mathcal{M}_Y}.
#' @param Ymu the Frechet mean of \eqn{Y}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param error.rho a correlation parameter for \eqn{\varepsilon}.
#' @param error.std a standard deviation parameter for \eqn{\varepsilon}.
#' 
#' @return an \eqn{n\times m} matrix of random errors.
error.generate = function(n,Yspace,Ymu,Ydim,error.rho=0.5,error.std=1){
  
  if (Yspace %in% c('Euclid','simplex','sphere','SPD.LogEuclid','SPD.AffInv')){
    basis = basis.manifold(Ymu,Ydim,Yspace)
    dim = nrow(basis)
    
    error.cov = matrix(error.rho,dim,dim) + diag(1-error.rho,dim)
    error.cov = error.cov * error.std^2 / dim
    error = MASS::mvrnorm(n,rep(0,dim),error.cov)
    error = error %*% basis
  } else{
    m = length(Ymu)
    const = 5/(m-1)
    error.cov = outer(1:m,1:m,function(x,y){error.rho^(const*abs(x-y))})
    error.cov = error.cov * error.std^2
    error = MASS::mvrnorm(n,rep(0,m),error.cov)
    
    if (Yspace=='BayesHilbert'){
      rho2 = error.rho^(const)
      error.std.correction = sqrt(1-1/(m^2 * (1-rho2))*(m*(1+rho2) - 2*(rho2-rho2^(m+1))/(1-rho2)))
      error = error/error.std.correction
      error = error - rowMeans(error)
    }
  }
  
  return(error)
}


############################################################
############################################################
### Generate simulation data


#' @title Generate linear regression data
#' 
#' @description
#' Generate linear regression data.
#' The seed for generating Hilbert-Schmidt operators \eqn{\mathfrak{B}_j} is fixed as zero, ensuring that each \eqn{\mathfrak{B}_j} is the same across simulations.
#' 
#' @param n the number of data points
#' @param Xspaces a \eqn{p} vector of the names of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.
#' @param Yspace the name of the underlying space \eqn{\mathcal{M}_Y} of \eqn{Y}.
#' @param Xdims a \eqn{p} vector of the intrinsic dimension of \eqn{\mathcal{M}_j}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param proper.indices an index set \eqn{\mathcal{S}=\{1\le j\le p : \mathfrak{B}_j\neq0\}}.
#' @param beta.norm a \eqn{p} vector (or a single integer) of the \eqn{\|\mathfrak{B}_j\|_{\mathcal{HS}}}, used only for \eqn{j\in\mathcal{S}}. If an integer is provided, the same value is applied to all \eqn{\|\mathfrak{B}_j\|_{\mathcal{HS}}}.
#' @param Xrho a correlation parameter for \eqn{X_j} ranging from -1 to 1, with a default value of 0.5.
#' @param Xsigma a common standard deviation parameter for \eqn{X_j}, with a default value of 1.
#' @param error.rho a correlation parameter for \eqn{\varepsilon}.
#' @param error.std a standard deviation parameter for \eqn{\varepsilon}.
#' @param ngrid the number of grids. This is only used for infinite-dimensional spaces, see \code{\link{Check.manifold}}: "functional", "BayesHilbert" and "Wasserstein" spaces.
#' @param seed a random seed, which must be greater than 1.
#' 
#' @return a list of data containing:
#'    \describe{
#'       \item{X}{a list of manifold-valued covariates, see \code{\link{covariates.generate}}.}
#'       \item{Y}{an \eqn{n\times m} matrix of manifold-valued responses.}
#'       \item{p}{the number of \eqn{X_j}.}
#'       \item{beta}{a \eqn{p} list of Hilbert-Schmidt operators, see \code{\link{tensor.beta.generate}}.}
#'       \item{Ymu}{the Frechet mean \eqn{\mu_Y} of \eqn{Y}.}
#'       \item{Xmu}{a \eqn{p} list of the Frechet means \eqn{\mu_j} of \eqn{X_j}.}
#'       \item{LogY}{the Riemannian logarithmic transformation \eqn{Log_{\mu_Y}Y}.}
#'       \item{LogX}{a \eqn{p} list of \eqn{Log_{\mu_j}X_j}, the Riemannian logarithmic transformations of \eqn{X_j}.}
#'       \item{Xbeta.each}{a \eqn{p} list of \eqn{\mathfrak{B}_j(Log_{\mu_j}X_j)}.}
#'       \item{Xbeta}{the sum of Xbeta.each, i.e., \eqn{\sum_{j=1}^p\mathfrak{B}_j(Log_{\mu_j}X_j)}.}
#'       \item{ExpXbeta}{the Riemannian exponential of Xbeta, i.e., \eqn{Exp_{\mu_Y}\Big(\sum_{j=1}^p\mathfrak{B}_j(Log_{\mu_j}X_j)\Big)}.}
#'       \item{...}{other input parameters.}
#' }
#' @export
LM.data.generate = function(n,Xspaces,Yspace,Xdims,Ydim,proper.indices=NULL,beta.norm=1,Xrho=0.5,Xsigma=1,
                            error.rho=0.5,error.std=1,ngrid=100,seed=1){
  
  # generate Xmu and Ymu
  Xmu.list = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j],ngrid)})
  Ymu = Ymu.generate(Ydim,Yspace,ngrid)
  
  # generate the same beta for each simulation
  set.seed(0)
  beta = tensor.beta.generate(Xspaces,Yspace,Xmu.list,Ymu,Xdims,Ydim,proper.indices,beta.norm)
  beta.norm = beta[['beta.norm']]
  
  # generate X and error
  set.seed(seed)
  X = covariates.generate(n,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)$X
  error = error.generate(n,Yspace,Ymu,Ydim,error.rho,error.std)
  
  # compute Xbeta
  LogX = lapply(1:length(Xdims),function(j){RieLog.manifold(Xmu.list[[j]],X[[j]],Xspaces[j])})
  Xbeta.each = lapply(1:length(Xdims),function(j){operator.tensor(beta[[j]],LogX[[j]])})
  Xbeta = Reduce('+',Xbeta.each)
  
  # make Y
  LogY = Xbeta + error
  ExpXbeta = RieExp.manifold(Ymu,Xbeta,Yspace)
  Y = RieExp.manifold(Ymu,LogY,Yspace)
  
  data = list(X=X,Y=Y,beta=beta,error=error,Ymu=Ymu,Xmu.list=Xmu.list,
              LogY=LogY,LogX=LogX,Xbeta.each=Xbeta.each,Xbeta=Xbeta,ExpXbeta=ExpXbeta,
              n=n,p=length(Xdims),Xspaces=Xspaces,Yspace=Yspace,Xdims=Xdims,Ydim=Ydim,
              proper.indices=proper.indices,beta.norm=beta.norm,Xrho=Xrho,Xsigma=Xsigma,
              error.rho=error.rho,error.std=error.std,seed=seed)
  return(data)
}



