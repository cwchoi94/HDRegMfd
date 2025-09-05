
library(HDRegMfd)
library(hms)
library(caret)
library(pROC)

# must include normalizer=(log10.cpm.max-log10.cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './GLM_code/NHANES/data/'
result.path0 = './GLM_code/NHANES/result/'
result.path = paste0(result.path0,'density/')
dir.create(result.path,showWarnings=FALSE)

load(paste0(data.path0,'PA_density_final.Rdata'))


link = 'binomial'
cv.type = 'AIC'
gamma = 0

# basic parameters
kfold = 5
proper.indices = c(40) # PA
lambda.list = seq(0.001,0.15,length.out=100)
Xdim.max.list = 2:4
R.list = c(100,200)
penalty.list = c('LASSO','SCAD','MCP')


split.Y.Xdata = function(PA,Covariate,X.simplex,Ycol){
  Y = Covariate[,Ycol,drop=FALSE]
  Covariate = Covariate[,!(colnames(Covariate) %in% c(Ycol))]
  p2 = ncol(Covariate)
  
  Xdata = lapply(1:p2,function(j){matrix(Covariate[,j],ncol=1)})
  Xdata[[p2+1]] = X.simplex
  Xdata[[p2+2]] = PA
  Xdata[['p']] = p2+2
  Xdata[['spaces']] = c(rep('Euclid',p2),'simplex','BayesHilbert')
  
  Y = as.matrix(Y)
  colnames(Y) = NULL
  
  return(list(Y=Y,Xdata=Xdata))
}


row.list = c('accuracy','auc','specificity','sensitivity','ratio.train','ratio.test')
col.list = c(penalty.list,'only PA')


# Compute classification error
start.time=Sys.time()

print(Ycols)
for (Ycol in Ycols){
  tmp = split.Y.Xdata(PA,Covariate,X.simplex,Ycol)
  Y = tmp$Y
  Xdata = tmp$Xdata
  p = Xdata[['p']]
  Xspaces = Xdata[['spaces']]
  
  set.seed(241101)
  groups = createFolds(Y,kfold)
  
  result.all = array(0,dim=c(kfold,6,4))
  
  for(k in 1:kfold){
    result = matrix(0,6,4)
    rownames(result) = row.list
    colnames(result) = col.list
    
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
      print(paste0('Ycol: ',Ycol,', iter: ',k,'/',kfold,', penalty: ',penalty,', runtime: ',runtime))
    }
    
    ####################################################################################
    # only PM
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
    result[,'only PA'] = result.each
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('Ycol: ',Ycol,', iter: ',k,'/',kfold,', penalty: ','only wind',', runtime: ',runtime))
    ####################################################################################
    
    result.all[k,,] = result
  }
  
  result = colMeans(result.all)
  rownames(result) = row.list
  colnames(result) = col.list
  write.csv(result,paste0(result.path,'class_result_avg_density_',Ycol,'.csv'))
}


# accuracy and AUC
accuracy = matrix(0,5,4)
rownames(accuracy) = Ycols
colnames(accuracy) = col.list
AUC = matrix(0,5,4)
rownames(AUC) = Ycols
colnames(AUC) = col.list
for (Ycol in Ycols){
  result = read.csv(paste0(result.path,'class_result_avg_density_',Ycol,'.csv'))
  rownames(result) = result[,1]
  result = result[,-1]
  
  accuracy[Ycol,] = unlist(result['accuracy',])
  AUC[Ycol,] = unlist(result['auc',])
}

accuracy = round(accuracy,4)
AUC = round(AUC,4)

write.csv(accuracy,paste0(result.path0,'Accuracy_density_all.csv'))
write.csv(AUC,paste0(result.path0,'AUC_density_all.csv'))

##################################################################


result.path0 = './GLM_code/NHANES/result/'
result.path = paste0(result.path0,'density/')

for (Ycol in Ycols){
  result = read.csv(paste0(result.path,'class_result_avg_density_',Ycol,'.csv'))
  
  print(Ycol)
  print(result)
}


read.csv(paste0(result.path0,'Accuracy_density_all.csv'))
read.csv(paste0(result.path0,'AUC_density_all.csv'))






