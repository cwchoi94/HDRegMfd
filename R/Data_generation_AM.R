



mX.functional = function(j,x,t){
  i = j %% 10
  if (i==1){
    z = log(2+t**2)*x**2
  }else if (i==2){
    z = 0.5*exp(x*t)
  }else if (i==3){
    z = sin(2*pi*(x**2+1)*t)
  }else if (i==4){
    z = 2*cos(2*pi*x*t)/(2 - sin(2*pi*x*t))
  }else if (i==5){
    z = sqrt(1+x**3-exp(-2*t))
  }else if (i==6){
    z = 2*(x-t**2)/(3-exp(x*t))
  }else if (i==7){
    z = exp(-2*x-t)
  }else if (i==8){
    z = log(2*x**3+1) * exp(-sqrt(t))
  }else if (i==9){
    z = (x+t**3-2)**2 / 3
  }else if (i==0){
    z = exp(-x)*cos(2*pi*(x+t))
  }
  
  return(z)
}

# xall = seq(0,1,length.out=20)
# tmp = t(sapply(1:10,function(j){
#   y = sapply(xall,function(x){sapply(xall,function(t){mX.functional(j,x,t)})})
#   a = quantile(y,c(0,0.2,0.4,0.6,0.8,1))
#   b = quantile(y**2,c(0,0.2,0.4,0.6,0.8,1))
#   cbind(a,b)
# }))
# round(tmp[,1:6],3)
# round(tmp[,7:12],3)



mX.basis = function(j,k,x){
  k = min(k,6)
  i = j %% 10
  if (i==1){
    z = exp(-k*x**2)
  }else if (i==2){
    z = exp(-k/3) * sin(4*pi*(x-k/2))
  }else if (i==3){
    z = (2-cos(pi*x)**3) / (1+sqrt(k))
  }else if (i==4){
    z = abs(x-1/2)/(1+k**2) 
  }else if (i==5){
    z = sin(2*pi*(k-3/2)*(x**2+1))
  }else if (i==6){
    z = (k-1.5)**2/((1+k)**2) * log(1+ k +3*x)
  }else if (i==7){
    z = 2*k*x**3 / (2 + x**2 + k**2)
  }else if (i==8){
    z = (cosh(x+k)-1)/(1+sinh(x+k)) + x**2/(2+k)
  }else if (i==9){
    z = (1+(x+k/2)**2) / (3+(x+k)**2)
  }else if (i==0){
    z = x**2 * cos(2*pi*(k**2+1)*x) / (k**2 + exp(k*x))
  }
  
  return(z)
}


# xall = seq(0,1,length.out=20)
# tmp = Reduce(rbind,lapply(1:10,function(j){t(sapply(1:6,function(k){
#   y = sapply(xall,function(x){mX.basis(j,k,x)})
#   a = quantile(y,c(0,0.2,0.4,0.6,0.8,1))
#   b = quantile(y**2,c(0,0.2,0.4,0.6,0.8,1))
#   cbind(a,b)
# }))}))
# j=8
# round(tmp[6*(j-1)+1:6,1:6],3)
# round(tmp[6*(j-1)+1:6,7:12],3)



######################################################
######################################################
### Generate additive mean functions


#' @title Generate an additive mean function 
#' 
#' @description
#' Generate an additive mean function \eqn{\mathfrak{m}: [0,1] \to T_{\mu_Y}\mathcal{M}_Y}, where \eqn{\mu_Y} is the Frechet mean of \eqn{Y}.
#' 
#' @param Xi.each an \eqn{n} vector of scores supported on \eqn{[0,1]}.
#' @param Yspace the name of the underlying space \eqn{\mathcal{M}_Y} of \eqn{Y}.
#' @param Ymu an \eqn{m} vector of the Frechet mean of \eqn{Y}.
#' @param Ydim the intrinsic dimension of \eqn{\mathcal{M}_Y}.
#' 
#' @return an \eqn{n\times m} matrix of an additive mean \eqn{\mathfrak{m}: [0,1] \to T_{\mu_Y}\mathcal{M}_Y}.
add.mean.generate.each = function(j,Xi.each,Yspace,Ymu,Ydim=1){
  
  if (Yspace %in% c('functional','BayesHilbert','Wasserstein')){
    t.all = seq(0,1,length.out=length(Ymu))
    add.mean = sapply(t.all,function(t){mX.functional(j,Xi.each,t)})
    if (Yspace=='BayesHilbert'){
      add.mean = clr(inv.clr.BayesHilbert(add.mean))
    }
  }else{
    Ybasis = basis.manifold(Ymu,Ydim,Yspace)
    add.mean.tmp = sapply(1:nrow(Ybasis),function(l){mX.basis(j,l,Xi.each) * b.ftn(1,l)}) # see 'b.ftn' in Data_generation_LM.R
    add.mean = add.mean.tmp %*% Ybasis
  }
  return(add.mean)
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
add.mean.generate = function(Xi,Yspace,Ymu,Ydim,all.indices){
  
  s = length(all.indices)
  
  add.mean.list = lapply(1:s,function(j){
    ind = all.indices[j]
    add.mean.generate.each(j,Xi[,ind],Yspace,Ymu,Ydim)
  })
  
  return(add.mean.list)
}




############################################################
############################################################
### Generate simulation data


#' @title Generate additive regression data
#' 
#' @description
#' Generate additive regression data.
#' The seed for generating additive mean functiions \eqn{m_{jk}^*} is fixed as zero, ensuring that each \eqn{m_{jk}^*} is the same across simulations.
#'
#' @inheritParams LM.data.generate
#' @inheritParams Transform.Score
#'
#' @param proper.ind.mat a \eqn{s\times 2} index matrix such that \eqn{\mathcal{S}=\{(j,k) : m_{jk}^*\neq0\}}.
#' @param add.mean.norm a \eqn{s} vector (or a single integer) of \eqn{\|m_{jk}^*\|_f}, used only for \eqn{(j,k)\in\mathcal{S}}. If an integer is provided, the same value is applied to all \eqn{(j,k)\in\mathcal{S}}.
#' @param seed a random seed, which must be greater than 1.
#' @param c.add.mean a list of pre-computed parameters to ensujre that \eqn{\E\|m_{jk}^*(\xi_{jk})\|^2 = (\text{add.mean.norm[j]})^2}, only computed if not provided.
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
AM.data.generate = function(n,Xspaces,Yspace,Xdims,Ydim,proper.ind.mat,add.mean.norm=1,Xrho=0.5,Xsigma=1,
                            error.rho=0.5,error.std=1,ngrid=50,transform='Gaussian',seed=1,c.add.mean=NULL){
  
  p = length(Xdims)
  s = nrow(proper.ind.mat)
  if(length(add.mean.norm)==1){add.mean.norm = rep(add.mean.norm,s)}
  
  # generate Xmu and Ymu
  Xmu.list = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j],ngrid)})
  Ymu = Ymu.generate(Ydim,Yspace,ngrid)
  
  # generate c.add.mean for ensuring \|m_{jk}^*\| = add.mean.norm[j']
  if (is.null(c.add.mean)){
    set.seed(0)
    n0 = 100000
    
    ## compute var(xi_{jk})
    Xdata.base0 = covariates.generate(n0,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)
    Xi.base0 = Xdata.base0$Xi.scaled
    sd.list = lapply(1:p,function(j){apply(Xi.base0[[j]],2,function(x){mean(x^2)})})
    
    ## generate new Xi and normalize
    Xdata.base = covariates.generate(n0,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)
    Xi.base = Xdata.base$Xi.scaled
    Xi.base = lapply(1:p,function(j){
      if(length(sd.list[[j]])>1){
        Xi.base[[j]] %*% diag(sd.list[[j]]^(-1/2))
      }else{
        Xi.base0[[j]]*(sd.list[[j]]^(-1/2))
      }})
    
    ## transform the normalized score
    object.transform = Transform.Score(Xi.base,transform,FALSE)
    index.mat = object.transform$index.mat
    Xi.base = predict(object.transform,Xi.base)
    
    ## compute matrices of additive mean functions
    all.indices = apply(proper.ind.mat,1,function(x){
      which(apply(index.mat[,-1],1,function(row){all(x==row)}))
    })
    add.mean.each = add.mean.generate(Xi.base,Yspace,Ymu,Ydim,all.indices)
    add.mean.mean = sapply(add.mean.each,colMeans)
    add.mean.std = sapply(1:s,function(j){
      tmp = add.mean.each[[j]] - matrix(add.mean.mean[,j],n0,length(Ymu),byrow=TRUE)
      sqrt(mean(norm.manifold(tmp,Ymu,Yspace)**2))
    })
    c.add.mean = list(sd.list=sd.list,mean=add.mean.mean,std=add.mean.std)
  }
  
  # generate X and error
  set.seed(seed)
  Xdata = covariates.generate(n,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)
  X = Xdata$X
  Xi = Xdata$Xi.scaled
  error = error.generate(n,Yspace,Ymu,Ydim,error.rho,error.std)
  
  ## transform the normalized score
  Xi = lapply(1:p,function(j){
    if(length(c.add.mean$sd.list[[j]])>1){
      Xi[[j]] %*% diag(c.add.mean$sd.list[[j]]^(-1/2))
    }else{
      Xi[[j]]*(c.add.mean$sd.list[[j]]^(-1/2))
    }})
  
  object.transform = Transform.Score(Xi,transform,FALSE)
  index.mat = object.transform$index.mat
  Xi.transform = predict(object.transform,Xi)
  
  ## compute matrices of additive mean functions
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  add.mean.each = add.mean.generate(Xi.transform,Yspace,Ymu,Ydim,all.indices)
  
  ## centering and constant adjustment
  for (j in 1:s){
    add.mean.each[[j]] = (add.mean.each[[j]] - matrix(c.add.mean$mean[,j],n,length(Ymu),byrow=TRUE)) / c.add.mean$std[j] * add.mean.norm[j]
  }
  add.mean = Reduce('+',add.mean.each)
  
  # make Y
  LogY = add.mean + error
  Exp.add.mean = RieExp.manifold(Ymu,add.mean,Yspace)
  Y = RieExp.manifold(Ymu,LogY,Yspace)
  
  proper.ind.mat.all = cbind(all.indices,proper.ind.mat,add.mean.norm)
  colnames(proper.ind.mat.all) = c('index','j','k','add.mean.norm')
  
  data = list(X=X,Y=Y,Xi=Xi,error=error,Ymu=Ymu,Xmu.list=Xmu.list,c.add.mean=c.add.mean,
              LogY=LogY,error=error,add.mean=add.mean,add.mean.each=add.mean.each,Exp.add.mean=Exp.add.mean,
              n=n,p=length(Xdims),Xspaces=Xspaces,Yspace=Yspace,Xdims=Xdims,Ydim=Ydim,
              proper.ind.mat=proper.ind.mat,proper.ind.mat.all=proper.ind.mat.all,
              Xrho=Xrho,Xsigma=Xsigma,error.rho=error.rho,error.std=error.std,seed=seed)
  return(data)
}




AM.data.generate2 = function(n,Xspaces,Yspace,Xdims,Ydim,proper.ind.mat,add.mean.norm=1,Xrho=0.5,Xsigma=1,
                            error.rho=0.5,error.std=1,ngrid=50,transform='Gaussian',seed=1,c.add.mean=NULL){
  
  p = length(Xdims)
  s = nrow(proper.ind.mat)
  if(length(add.mean.norm)==1){add.mean.norm = rep(add.mean.norm,s)}
  
  # generate Xmu and Ymu
  Xmu.list = lapply(1:length(Xdims),function(j){Xmu.generate(Xdims[j],Xspaces[j],ngrid)})
  Ymu = Ymu.generate(Ydim,Yspace,ngrid)
  
  # generate c.add.mean for ensuring \|m_{jk}^*\| = add.mean.norm[j']
  if (is.null(c.add.mean)){
    set.seed(0)
    n0 = 100000
    
    Xdata.base = covariates.generate(n0,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)
    Xi.base = Xdata.base$Xi.scaled
    
    ## transform the normalized score
    object.transform = Transform.Score(Xi.base,transform,FALSE)
    index.mat = object.transform$index.mat
    Xi.base = predict(object.transform,Xi.base)
    
    ## compute matrices of additive mean functions
    all.indices = apply(proper.ind.mat,1,function(x){
      which(apply(index.mat[,-1],1,function(row){all(x==row)}))
    })
    add.mean.each = add.mean.generate(Xi.base,Yspace,Ymu,Ydim,all.indices)
    add.mean.mean = sapply(add.mean.each,colMeans)
    add.mean.std = sapply(1:s,function(j){
      tmp = add.mean.each[[j]] - matrix(add.mean.mean[,j],n0,length(Ymu),byrow=TRUE)
      sqrt(mean(norm.manifold(tmp,Ymu,Yspace)**2))
    })
    c.add.mean = list(mean=add.mean.mean,std=add.mean.std)
  }
  
  # generate X and error
  set.seed(seed)
  Xdata = covariates.generate(n,Xspaces,Xmu.list,Xdims,Xrho,Xsigma)
  X = Xdata$X
  Xi = Xdata$Xi.scaled
  error = error.generate(n,Yspace,Ymu,Ydim,error.rho,error.std)
  
  ## transform the normalized score
  object.transform = Transform.Score(Xi,transform,FALSE)
  index.mat = object.transform$index.mat
  Xi.transform = predict(object.transform,Xi)
  
  ## compute matrices of additive mean functions
  all.indices = apply(proper.ind.mat,1,function(x){
    which(apply(index.mat[,-1],1,function(row){all(x==row)}))
  })
  add.mean.each = add.mean.generate(Xi.transform,Yspace,Ymu,Ydim,all.indices)
  
  ## centering and constant adjustment
  for (j in 1:s){
    add.mean.each[[j]] = (add.mean.each[[j]] - matrix(c.add.mean$mean[,j],n,length(Ymu),byrow=TRUE)) / c.add.mean$std[j] * add.mean.norm[j]
  }
  add.mean = Reduce('+',add.mean.each)
  
  # make Y
  LogY = add.mean + error
  Exp.add.mean = RieExp.manifold(Ymu,add.mean,Yspace)
  Y = RieExp.manifold(Ymu,LogY,Yspace)
  
  proper.ind.mat.all = cbind(all.indices,proper.ind.mat,add.mean.norm)
  colnames(proper.ind.mat.all) = c('index','j','k','add.mean.norm')
  
  data = list(X=X,Y=Y,Xi=Xi,error=error,Ymu=Ymu,Xmu.list=Xmu.list,c.add.mean=c.add.mean,
              LogY=LogY,error=error,add.mean=add.mean,add.mean.each=add.mean.each,Exp.add.mean=Exp.add.mean,
              n=n,p=length(Xdims),Xspaces=Xspaces,Yspace=Yspace,Xdims=Xdims,Ydim=Ydim,
              proper.ind.mat=proper.ind.mat,proper.ind.mat.all=proper.ind.mat.all,
              Xrho=Xrho,Xsigma=Xsigma,error.rho=error.rho,error.std=error.std,seed=seed)
  return(data)
}

