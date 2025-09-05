
library(HDRegMfd)
library(hms)
library(caret)
library(pROC)


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))


link = 'binomial'
cv.type = 'AIC'
gamma = 0

p = Xdata[['p']]
Xspaces = Xdata[['spaces']]

dim(Y)
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# basic parameters
kfold = 5
proper.indices = c(12) # wind
lambda.list = seq(0.001,0.15,length.out=100)
Xdim.max.list = 2:5
R.list = c(100,200)
penalty.list = c('LASSO','SCAD','MCP')


# Compute ASPE
start.time=Sys.time()

set.seed(241101)
groups = createFolds(Y,kfold)

result.all = array(0,dim=c(kfold,6,4))

for(k in 1:kfold){
  result = matrix(0,6,4)
  rownames(result) = c('accuracy','auc','specificity','sensitivity','ratio.train','ratio.test')
  colnames(result) = c(penalty.list,'only wind')
  
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
  
  for (penalty in penalty.list){
    # compute prediction error
    model = GLM.CV(Xtrain,Ytrain,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list)
    Ypred = predict(model,Xtest)
    
    ## accuracy
    Yclass = 1.0 * (Ypred>=0.5)
    accuracy = mean(Ytest == Yclass)
    
    ## auc, sensitivity, specificity
    roc.result = roc(Ytest[,1] ~ Ypred[,1],levels=c(0,1),direction='<')
    tmp.result = coords(roc.result,'best',ret=c('specificity','sensitivity'))[1,]
    
    result.each = c(accuracy,roc.result$auc,tmp.result,mean(Ytrain),mean(Ytest))
    result.each = unlist(result.each)
    names(result.each) = NULL
    result[,penalty] = result.each
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('iter: ',k,'/',kfold,', penalty: ',penalty,', runtime: ',runtime))
  }
  
  ####################################################################################
  # only wind
  model = GLM.oracle.CV(Xtrain,Ytrain,link,proper.indices,cv.type,Xdim.max.list)
  Ypred = predict(model,Xtest)
  
  ## accuracy
  Yclass = 1.0 * (Ypred>=0.5)
  accuracy = mean(Ytest == Yclass)
  
  ## auc, sensitivity, specificity
  roc.result = roc(Ytest[,1] ~ Ypred[,1],levels=c(0,1),direction='<')
  tmp.result = coords(roc.result,'best',ret=c('specificity','sensitivity'))[1,]
  
  result.each = c(accuracy,roc.result$auc,tmp.result,mean(Ytrain),mean(Ytest))
  result.each = unlist(result.each)
  names(result.each) = NULL
  result[,'only wind'] = result.each
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('iter: ',k,'/',kfold,', penalty: ','only wind',', runtime: ',runtime))
  ####################################################################################
  
  
  result.all[k,,] = result
}


result = colMeans(result.all)
rownames(result) = c('accuracy','auc','specificity','sensitivity','ratio.train','ratio.test')
colnames(result) = c(penalty.list,'only wind')
write.csv(result,paste0(result.path0,'class_result_avg.csv'))


##################################################################

result = read.csv(paste0(result.path0,'class_result_avg.csv'))
result





