

library(HDRegMfd)
# source('initialize.R')

n = 400
n2 = 100
n3 = 1000


tau = 0.5
error.var.type = 'homoscedastic'
# error.var.type = 'heteroscedastic'

error.type = 'normal'
error.type = 't'
error.type = 'Cauchy'


df = 3


get.dim.each = function(x){
  if (x=='Euclid'){
    dim = 1
  } else if (x=='simplex'){
    dim = 3
  } else if (x=='sphere'){
    dim = 3
  } else if (x=='SPD.LogEuclid'){
    dim = 4
  } else if (x=='functional'){
    dim = 50
  } else if (x=='BayesHilbert'){
    dim = 40
  } else if (x=='Wasserstein'){
    dim = 5
  }
  return(dim)
}

Xrho = 0.5
Xsigma = 1
beta0.norm = 1


Yspace = 'Euclid'
Ydim = 1
error.rho = 0.2
error.sigma = 1


Xspaces = c(rep('Euclid',10),rep('simplex',10),rep('sphere',10),rep('SPD.LogEuclid',0),rep('functional',2),rep('BayesHilbert',0),rep('Wasserstein',2))
Xdims = sapply(1:length(Xspaces),function(j){get.dim.each(Xspaces[j])})
proper.indices = c(10,20,30,32,34)
beta.norm = rep(0,length(Xspaces))
beta.norm[proper.indices] = c(0.5,0.5,0.5,2,2)


p = 100
q = 4
s = 3

Xspaces = c(rep('functional',q/2),rep('Wasserstein',q/2),rep('Euclid',(p-q)/2),rep('simplex',(p-q)/4),rep('sphere',(p-q)/4))
Xdims = sapply(Xspaces,get.dim.each,USE.NAMES=FALSE)
proper.indices = c(2,if(s>1){sapply(1:(s-1),function(j){q+floor((p-q)/(s-1)*j)})}else{})
beta.norm = 1.0

proper.indices
Xspaces[proper.indices]


iter = 6

data = QM.data.generate(n,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter)
datanew = QM.data.generate(n2,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+1000,c.beta=data$c.beta)
datatest = QM.data.generate(n3,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+2000,c.beta=data$c.beta)



# norm of LogX
sapply(proper.indices,function(j){mean(norm.manifold(data$LogX[[j]],data$Xmu[[j]],data$Xspaces[j])^2)})

# norm of beta
sapply(proper.indices,function(j){norm.tensor(data$beta[[j]])})

# norm of Xbeta
sapply(proper.indices,function(j){sqrt(mean(norm.manifold(data$Xbeta.each[[j]])^2))})




