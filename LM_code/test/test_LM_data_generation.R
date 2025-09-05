

library(HDRegMfd)

n = 200
n2 = 100
n3 = 1000

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
Xsigma = 0.5
Yspace = 'BayesHilbert'
Ydim = 50
# Yspace = 'SPD.LogEuclid'
# Ydim = 9
beta.norm = 1
error.rho = 0.6
error.sigma = 1

Xspaces = c(rep('Euclid',100),rep('simplex',50),rep('sphere',46),rep('SPD.LogEuclid',0),rep('functional',2),rep('BayesHilbert',0),rep('Wasserstein',2))
Xdims = sapply(1:length(Xspaces),function(j){get.dim.each(Xspaces[j])})
proper.indices = c(20,40,60,80,120,140,160,180,198,200)

Xspaces = c(rep('Euclid',10),rep('simplex',10),rep('sphere',10),rep('SPD.LogEuclid',0),rep('functional',2),rep('BayesHilbert',0),rep('Wasserstein',2))
Xdims = sapply(1:length(Xspaces),function(j){get.dim.each(Xspaces[j])})
proper.indices = c(10,20,30,32,34)


n = 200
p = 10
q = 2
s = 3

Xspaces = c(rep('functional',q/2),rep('Wasserstein',q/2),rep('Euclid',(p-q)/2),rep('simplex',(p-q)/4),rep('sphere',(p-q)/4))
Xdims = sapply(Xspaces,get.dim.each,USE.NAMES=FALSE)
proper.indices = c(c(1,q/2+1),sapply(1:(s-2),function(j){q+floor((p-q)/(s-2)*j)}))
beta.norm = 1

proper.indices
Xspaces[proper.indices]

Yspace = 'Euclid'
Ydim = 1


data = LM.data.generate(n,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=1)
datanew = LM.data.generate(n2,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=1+1000)
datatest = LM.data.generate(n3,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=1+1000)


# norm of LogX
sapply(proper.indices,function(j){mean(norm.manifold(data$LogX[[j]],data$Xmu[[j]],data$Xspaces[j])^2)})

# norm of beta
sapply(proper.indices,function(j){norm.tensor(data$beta[[j]])})

# norm of Xbeta
sapply(proper.indices,function(j){sqrt(mean(norm.manifold(data$Xbeta.each[[j]])^2))})




