
library(HDRegMfd)
library(hms)


# must include normalizer=(Y.max-Y.min)^(-1) in "inv.clr.density" function.
data.path0 = './LM_code/KLIPS/data/'
result.path0 = './LM_code/KLIPS/result/'
load(paste0(data.path0,'KLIPS.Rdata'))
load(paste0(data.path0,'Covariate.Rdata'))
colnames(Covariate)

Y = Y.density
Yspace = 'BayesHilbert'

n = nrow(Y)
p = Xdata[['p']]
Xspaces = Xdata[['spaces']]
Xdims = sapply(1:Xdata[['p']],function(j){ncol(Xdata[[j]])})


dim(Y)
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})
Xdims
Xspaces

# basic parameters
kfold = 5
seed = 10000
gamma = 0

lambda.list = seq(0.001,0.05,length.out=50)
Xdim.max.list = 1:10
R.list = seq(5,10,length.out=2)
penalty.list = c('LASSO','SCAD','MCP')


# Compute ASPE
start.time=Sys.time()

error = matrix(0,5,3)
colnames(error) = penalty.list

set.seed(240101)
all.index = sample(1:n)

for (penalty in penalty.list){
  for(i in 1:5){
    idx1 = as.integer(n/5*(i-1))+1
    idx2 = as.integer(n/5*i)
    test.index = all.index[idx1:idx2]
    
    # test data
    Ytest = Y[test.index,]
    Xtest = lapply(1:p,function(j){matrix(Xdata[[j]][test.index,],ncol=Xdims[j])})
    Xtest[['p']] = p
    Xtest[['spaces']] = Xspaces
    
    # train data
    Ytrain = Y[-test.index,]
    Xtrain = lapply(1:p,function(j){matrix(Xdata[[j]][-test.index,],ncol=Xdims[j])})
    Xtrain[['p']] = p
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

