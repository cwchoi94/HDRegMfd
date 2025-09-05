
library(HDRegMfd)
library(hms)

# must include normalizer=(cpm.max-cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './LM_code/NHANES/data/'
result.path0 = './LM_code/NHANES/result/'
load(paste0(data.path0,'PA_density.Rdata'))
load(paste0(data.path0,'covariate_final.Rdata'))
Covariate = Covariate[,-1] # remove subject

dim(Y)
dim(Covariate)
p = ncol(Covariate)

X.simplex = as.matrix(Covariate[,15:17])
X.simplex = X.simplex/rowSums(X.simplex)
rowSums(X.simplex)

Covariate = Covariate[,-(15:17)]
p = ncol(Covariate)

Yall = Y
Xall = lapply(1:p,function(j){matrix(Covariate[,j],ncol=1)})
Xall[[p+1]] = X.simplex
Xall[['p']] = p+1
Xall[['spaces']] = c(rep('Euclid',p),'simplex')

# basic parameters
Yspace = 'BayesHilbert'

kfold = 5
seed = 10000
gamma = 0

lambda.list = seq(0.001,0.1,length.out=100)
Xdim.max.list = 1:2
R.list = seq(5,10,length.out=2)
penalty.list = c('LASSO','SCAD','MCP')


# model fitting
start.time=Sys.time()


for (penalty in penalty.list){
  model = LM.kfold(Xdata,Y,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list)
  save(model,file=paste0(result.path0,'NHANES_',penalty,'.Rdata'))
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('penalty: ',penalty,', runtime: ',runtime))
}


