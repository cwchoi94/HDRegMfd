
library(HDRegMfd)
library(hms)
library(caret)
library(pROC)


data.path0 = './QM_code/PM/data/'
result.path0 = './QM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))


cv.type = 'BIC'
gamma = 0
tau.list = seq(0.1,0.9,0.1)
tau.list = c(0.25,0.5,0.75)

kernel = 'Gaussian'
cv.const = 2
h = NULL


p = Xdata[['p']]
Xspaces = Xdata[['spaces']]

dim(Y)
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# basic parameters
kfold = 5
proper.indices = c(12) # wind
lambda.list = seq(0.0005,0.01,length.out=20)
Xdim.max.list = 1:5
penalty.list = c('LASSO','SCAD','MCP')


quantile.loss = function(y,yhat,tau=0.5){
  u = y - yhat
  loss = ifelse(u>=0, tau*u,(tau-1)*u)
  return(mean(loss))
}


# Compute ASPE under quantile loss
start.time=Sys.time()

set.seed(250501)
groups = createFolds(Y,kfold)

result.all = array(0,dim=c(kfold,length(tau.list),4))

for(k in 1:kfold){
  
  if ((k<=9) & (kfold>=10)){
    test.index = groups[[paste0('Fold0',k)]]
  }else{
    test.index = groups[[paste0('Fold',k)]]
  }
  
  # test data
  Ytest = Y[test.index,,drop=FALSE]
  Xtest = lapply(1:p,function(j){Xdata[[j]][test.index,,drop=FALSE]})
  Xtest[['p']] = p
  Xtest[['spaces']] = Xspaces
  
  # train data
  Ytrain = Y[-test.index,,drop=FALSE]
  Xtrain = lapply(1:p,function(j){Xdata[[j]][-test.index,,drop=FALSE]})
  Xtrain[['p']] = p
  Xtrain[['spaces']] = Xspaces
  
  for (idx in 1:length(tau.list)){
    tau = tau.list[idx]
    
    result = matrix(0,1,4)
    rownames(result) = NULL
    colnames(result) = c(penalty.list,'only wind')
    
    for (penalty in penalty.list){
      
      # compute prediction error
      model = QM.CV(Xtrain,Ytrain,tau,h,kernel,cv.type,penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const)
      Ypred = predict(model,Xtest)
      
      ## loss
      loss = quantile.loss(Ytest,Ypred,tau)
      result[,penalty] = loss
      
      runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
      print(paste0('iter: ',k,'/',kfold,', tau: ',tau,', penalty: ',penalty,', runtime: ',runtime))
    }
    
    ####################################################################################
    # only wind
    model = QM.oracle.CV(Xdata,Y,tau,h,kernel,proper.indices,cv.type,Xdim.max.list,cv.const=cv.const)
    Ypred = predict(model,Xtest)
    
    ## loss
    loss = quantile.loss(Ytest,Ypred,tau)
    result[,'only wind'] = loss
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('iter: ',k,'/',kfold,', tau: ',tau,', penalty: ','only wind',', runtime: ',runtime))
    
    ####################################################################################
    
    result.all[k,idx,] = result
  }
}


result = colMeans(result.all)
result = round(result,4)
rownames(result) = tau.list
colnames(result) = c(penalty.list,'only wind')
write.csv(result,paste0(result.path0,'result_avg.csv'))


##################################################################

result = read.csv(paste0(result.path0,'result_avg.csv'))
result




