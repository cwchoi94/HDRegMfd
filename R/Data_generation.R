### data generation codes for simulations.


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
#' 
#' @return the Frechet mean \eqn{\mu_j} of \eqn{X_j}, an \eqn{m} vector.
Xmu.generate = function(m,space){
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
    m = 100
    mu = rep(0,m)
  } else if (space=='BayesHilbert'){
    m = 100
    mu = inv.clr.BayesHilbert(rep(0,m))[1,]
  } else if (space=='Wasserstein'){
    m = 100
    t = seq(0,1,length.out=m+2)[2:(m+1)]
    mu = 10*t
  }
  return(mu)
}


#' @title Generate the Frechet mean for \eqn{Y}
#' 
#' @description
#' This code is the same as \code{\link{Xmu.generate}}.
#' 
#' @param m The length of each observation. Typically, this equals the dimension, but it can differ in some cases.
#' @param space The name of space, see \code{\link{Check.manifold}}.
#' 
#' @return The Frechet mean \eqn{\mu_Y} of \eqn{Y}, an \eqn{m} vector.
Ymu.generate = function(m,space){
  mu = Xmu.generate(m,space)
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


#' @title Generate a manifold-valued covariate
#' 
#' @description
#' Given \eqn{j}, generate manifold-valued covariates \eqn{X_{ij}} for \eqn{1\le i\le n}.
#' 
#' @param Xi an \eqn{n\times m'} score matrix, where \eqn{m'} is the intrinsic dimension (the number of scores).
#' @param dim the dimension of \eqn{X}.
#' @param space the underlying spaces, see \code{\link{Check.manifold}}.
#' 
#' @return an \eqn{n\times m} matrix. Each row represents a manifold-valued covariate.
#' @export
covariates.generate.each = function(Xi,dim,space='Euclid'){
  Check.manifold(space)
  mu = Xmu.generate(dim,space)
  basis = basis.manifold(mu,ncol(Xi),space)
  
  if (space=='Wasserstein'){
    # restrict the range of Xi to ensure that LogXj is the quantile function.
    Xi = 2*pnorm(Xi) - 1 
    Xi = Xi %*% diag(sapply(1:ncol(Xi),function(k){w.ftn(k)/sqrt(2)}))
  } else{
    Xi = Xi %*% diag(sapply(1:ncol(Xi),function(k){w.ftn(k)}))
  }
  
  LogX = Xi %*% basis
  X = RieExp.manifold(mu,LogX,space)
  return(X)
}

# weight bound for w_jk
w.ftn = function(k){
  if (k<=4){
    w = k^(-1)
  } else if (k>4){
    w = k^(-1)
  }
  return(w)
}


#' @title Generate a list of manifold-valued covariates.
#' 
#' @description
#' Generate a list of manifold-valued covariates.
#' If the underlying space of \eqn{X_j} is finite-dimensional, the dimension of \eqn{X_j} is set to \eqn{D_j}.
#' If the underlying space of \eqn{X_j} is inifite-dimensional, the dimension of \eqn{X_j} is set to 100.
#' 
#' @param n the number of data points.
#' @param Xspaces a \eqn{p} vector specifying the underlying spaces of \eqn{X_j}.
#' @param dims a \eqn{p} vector specifying the dimensions of \eqn{X_j}.
#' @param Xrho a correlation parameter ranging from -1 to 1, with a default value of 0.5.
#' @param Xsigma a common standard deviation parameter, with a default value of 1.
#' 
#' @return a list of data containing:
#'    \describe{
#'       \item{j}{a \eqn{p}-list of generated data, where each \eqn{j}th element is an \eqn{n\times D_j} matrix.}
#'       \item{Xspaces}{a \eqn{p} vector specifying the underlying spaces of \eqn{X_j}.}
#'       \item{p}{the number of \eqn{X_j}.}
#' }
covariates.generate = function(n,Xspaces,dims,Xrho=0.5,Xsigma=1){
  # compute intrinsic dimension
  p = length(dims)
  dims_ = sapply(1:p,function(i){
    if(Xspaces[i] %in% c('simplex','sphere')){
      dims[i]-1
    }else if (Xspaces[i]=='SPD.LogEuclid'){
      (dims[i]+sqrt(dims[i]))/2
    }else{
      dims[i]
    }})
  
  # generate scores
  Zeta = lapply(1:p,function(j){covariates.generate.real(n,dims_[j],0,Xsigma)})
  Xi = lapply(1:p,function(j){
    if (j==1){
      I2 = make.nonsym.Idmat(dims_[j+1],dims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j+1]] %*% I2)
    } else if (j==p){
      I1 = make.nonsym.Idmat(dims_[j-1],dims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j-1]] %*% I1)
    } else{
      I1 = make.nonsym.Idmat(dims_[j-1],dims_[j])
      I2 = make.nonsym.Idmat(dims_[j+1],dims_[j])
      Xi.each = Zeta[[j]] + Xrho * (Zeta[[j-1]] %*% I1 + Zeta[[j+1]] %*% I2)
    }
    return(Xi.each)
  })
  
  # generate X_j
  Xdata = lapply(1:p,function(j){covariates.generate.each(Xi[[j]],dims[j],Xspaces[j])})
  Xdata[['spaces']] = Xspaces
  Xdata[['p']] = p
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
#' @param Xdim the intrinsic dimension of \eqn{\mathcal{M}_X}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param beta.norm the tensor norm, with a default value of 1.
#' 
#' @return a Hilbert-Schmidt operators (an tensor element) in \eqn{T_{\mu_X}\mathcal{M}_X \otimes T_{\mu_Y}\mathcal{M}_Y} generated by \code{\link{make.tensor}}.
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
tensor.beta.generate = function(Xspaces,Yspace,Xdims,Ydim,proper.indices=NULL,beta.norm=1){
  p = length(Xspaces)
  if(is.null(proper.indices)){proper.indices = seq(p)}
  if(length(beta.norm)==1){beta.norm = rep(beta.norm,p)}
  beta.norm[-proper.indices] = 0
  beta.norm = beta.norm[1:p]
  
  beta.list = list()
  for (j in 1:p){
    beta.list[[j]] = tensor.beta.generate.each(Xspaces[[j]],Yspace,Xdims[j],Ydim,beta.norm[j])
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
#' @param space the underlying space of \eqn{\mathcal{M}_Y}.
#' @param dim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' @param error.rho a correlation parameter for \eqn{\varepsilon}.
#' @param error.std a standard deviation parameter for \eqn{\varepsilon}.
#' 
#' @return an \eqn{n\times m} matrix of random errors.
error.generate = function(n,space,dim,error.rho=0.5,error.std=1){
  mu = Ymu.generate(dim,space)
  
  if (space %in% c('Euclid','simplex','sphere','SPD.LogEuclid','SPD.AffInv')){
    basis = basis.manifold(mu,dim,space)
    dim = nrow(basis)
    
    error.cov = matrix(error.rho,dim,dim) + diag(1-error.rho,dim)
    error.cov = error.cov * error.std^2 / dim
    error = MASS::mvrnorm(n,rep(0,dim),error.cov)
    error = error %*% basis
  } else{
    m = length(mu)
    const = 5/(m-1)
    error.cov = outer(1:m,1:m,function(x,y){error.rho^(const*abs(x-y))})
    error.cov = error.cov * error.std^2
    error = MASS::mvrnorm(n,rep(0,m),error.cov)
    
    if (space=='BayesHilbert'){
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
#' @param beta.norm a \eqn{p} vector (or a single integer) of the tensor norms, used only for \eqn{j\notin\mathcal{S}}. If an integer is provided, the same value is applied for all \eqn{X_j}.
#' @param Xrho a correlation parameter for \eqn{X_j} ranging from -1 to 1, with a default value of 0.5.
#' @param Xsigma a common standard deviation parameter for \eqn{X_j}, with a default value of 1.
#' @param error.rho a correlation parameter for \eqn{\varepsilon}.
#' @param error.std a standard deviation parameter for \eqn{\varepsilon}.
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
                            error.rho=0.5,error.std=1,seed=1){
  
  # generate the same beta for each simulation
  set.seed(0)
  beta = tensor.beta.generate(Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm)
  beta.norm = beta[['beta.norm']]
  
  # generate X and error
  set.seed(seed)
  X = covariates.generate(n,Xspaces,Xdims,Xrho,Xsigma)
  Xmu = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j])})
  error = error.generate(n,Yspace,Ydim,error.rho,error.std)
  
  # compute Xbeta
  LogX = lapply(1:length(Xdims),function(j){RieLog.manifold(Xmu[[j]],X[[j]],Xspaces[j])})
  Xbeta.each = lapply(1:length(Xdims),function(j){operator.tensor(beta[[j]],LogX[[j]])})
  Xbeta = Reduce('+',Xbeta.each)
  
  # make Y
  Ymu = Ymu.generate(Ydim,Yspace)
  LogY = Xbeta + error
  ExpXbeta = RieExp.manifold(Ymu,Xbeta,Yspace)
  Y = RieExp.manifold(Ymu,LogY,Yspace)
  
  data = list(X=X,Y=Y,beta=beta,error=error,Ymu=Ymu,Xmu=Xmu,
              LogY=LogY,LogX=LogX,Xbeta.each=Xbeta.each,Xbeta=Xbeta,ExpXbeta=ExpXbeta,
              n=n,p=length(Xdims),Xspaces=Xspaces,Yspace=Yspace,Xdims=Xdims,Ydim=Ydim,
              proper.indices=proper.indices,beta.norm=beta.norm,Xrho=Xrho,Xsigma=Xsigma,
              error.rho=error.rho,error.std=error.std,seed=seed)
  return(data)
}




#' @title Generate generalized linear regression data
#' 
#' @description
#' Generate generalized linear regression data.
#' The response \eqn{Y} given the covariates \eqn{X} follows an exponential family distribution, and thus \eqn{Y} is real-valued.
#' The seed for generating Hilbert-Schmidt operators \eqn{\mathfrak{B}_j} is fixed as zero, ensuring that each \eqn{\mathfrak{B}_j} is the same across simulations.
#' 
#' @param n the number of data points
#' @param Xspaces a \eqn{p} vector of the names of the underlying spaces \eqn{\mathcal{M}_j} of \eqn{X_j}.
#' @param Xdims a \eqn{p} vector of the intrinsic dimension of \eqn{\mathcal{M}_j}.
#' @param link a link function, see \code{\link{Check.link}}.
#' @param Ydim the intrinsic dimension of \eqn{Y}, typically 1.
#' @param proper.indices an index set \eqn{\mathcal{S}=\{1\le j\le p : \mathfrak{B}_j\neq0\}}.
#' @param beta.norm a \eqn{p} vector (or a single integer) of the \eqn{L^2} norms of \eqn{\| \mathfrak{B_j}(Log_{\mu_j}X_j) \|}, used only for \eqn{j\notin\mathcal{S}}. If an integer is provided, the same value is applied for all \eqn{X_j}.
#' @param beta0.norm an integer specifying the norm of \eqn{\beta_0^*}, with a default value of 1.
#' @param Xrho a correlation parameter for \eqn{X_j} ranging from -1 to 1, with a default value of 0.5.
#' @param Xsigma a common standard deviation parameter for \eqn{X_j}, with a default value of 1.
#' @param seed a random seed, which must be greater than 1.
#' @param c.beta a \eqn{p} vector of pre-computed parameters to ensure that \eqn{\mathbb{E}\| \mathfrak{B_j}(Log_{\mu_j}X_j) \|^2 = \text{beta.norm[j]}^2}, only computed if not provided.
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
GLM.data.generate = function(n,Xspaces,Xdims,link='binomial',Ydim=1,proper.indices=NULL,beta.norm=1,beta0.norm=1,Xrho=0.5,Xsigma=1,seed=1,c.beta=NULL){
  
  Yspace = 'Euclid'
  Check.link(link)
  if (link!='multinomial'){Ydim = 1}
  if(is.null(proper.indices)){proper.indices = seq(length(Xdims))}
  if(length(beta.norm)==1){beta.norm = rep(beta.norm,length(Xdims))}
  p = length(Xdims)
  s = length(proper.indices)
  
  # generate same beta for each simulation
  set.seed(0)
  beta = tensor.beta.generate(Xspaces,Yspace,Xdims,Ydim,proper.indices,1)
  beta.oracle = lapply(proper.indices,function(j){beta[[j]]})
  
  # generate nuisance X to set |B_j(LogX_j)|=beta.norm[j]
  if (is.null(c.beta)){
    Xmu = lapply(1:p,function(j){Xmu.generate(Xdims[j],Xspaces[j])})
    X.base = covariates.generate(10000,Xspaces,Xdims,Xrho,Xsigma)
    LogX.base = lapply(proper.indices,function(j){RieLog.manifold(Xmu[[j]],X.base[[j]],Xspaces[j])})
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
  X = covariates.generate(n,Xspaces,Xdims,Xrho,Xsigma)
  Xmu = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j])})
  
  # compute Xbeta
  LogX = lapply(1:length(Xdims),function(j){RieLog.manifold(Xmu[[j]],X[[j]],Xspaces[j])})
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
              Xmu=Xmu,LogX=LogX,Xbeta=Xbeta,Xbeta.each=Xbeta.each,
              n=n,p=length(Xdims),Xspaces=Xspaces,Xdims=Xdims,Ydim=Ydim,
              proper.indices=proper.indices,beta.norm=beta.norm,beta0.norm=beta0.norm,
              Xrho=Xrho,Xsigma=Xsigma,seed=seed)
  return(data)
}



