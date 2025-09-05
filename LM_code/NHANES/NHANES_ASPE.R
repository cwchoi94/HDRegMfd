
library(HDRegMfd)
library(hms)

# must include normalizer=(cpm.max-cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './LM_code/NHANES/data/'
result.path0 = './LM_code/NHANES/result/'
load(paste0(data.path0,'PA_density.Rdata'))
load(paste0(data.path0,'covariate_final.Rdata'))
Covariate = Covariate[,-1]

dim(Y)
dim(Covariate)
n = nrow(Y)
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
Xspaces = c(rep('Euclid',p),'simplex')
Xall[['spaces']] = Xspaces

# basic parameters
Yspace = 'BayesHilbert'

kfold = 5
seed = 10000
gamma = 0

lambda.list = seq(0.001,0.1,length.out=100)
Xdim.max.list = 1:2
R.list = seq(5,10,length.out=2)
penalty.list = c('LASSO','SCAD','MCP')


# Compute ASPE
start.time=Sys.time()

error = matrix(0,5,3)
colnames(error) = penalty.list

set.seed(240101)
all.index = sample(1:n)

for (penalty in penalty.list){
  beta.norm.each = matrix(0,5,length(Xspaces))
  for(i in 1:5){
    idx1 = as.integer(n/5*(i-1))+1
    idx2 = as.integer(n/5*i)
    test.index = all.index[idx1:idx2]
    
    # test data
    Ytest = Y[test.index,]
    Xtest = lapply(1:p,function(j){matrix(Covariate[test.index,j],ncol=1)})
    Xtest[[p+1]] = X.simplex[test.index,]
    Xtest[['p']] = p+1
    Xtest[['spaces']] = Xspaces
    
    # train data
    Ytrain = Y[-test.index,]
    Xtrain = lapply(1:p,function(j){matrix(Covariate[-test.index,j],ncol=1)})
    Xtrain[[p+1]] = X.simplex[-test.index,]
    Xtrain[['p']] = p+1
    Xtrain[['spaces']] = Xspaces
    
    # compute prediction error
    model = LM.kfold(Xtrain,Ytrain,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list)
    Ypred = predict(model,Xtest)
    error[i,penalty] = mean(dist.manifold(Ypred,Ytest,Yspace)^2)
    # error[i,penalty] = vector.norm(Ypred-Ytest,model$Ymu,Yspace,'L2')/sqrt(length(test.index))
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('penalty: ',penalty,', iter: ',i,'/',5,', runtime: ',runtime))
  }
}

ASPE = colMeans(error)

write.csv(ASPE,paste0(result.path0,'ASPE.csv'))


# load ASPE 
ASPE = read.csv(paste0(result.path0,'ASPE.csv'))
ASPE


