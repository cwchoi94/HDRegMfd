
library(HDRegMfd)
library(hms)

# must include normalizer=(log10.cpm.max-log10.cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './GLM_code/NHANES/data/'
result.path0 = './GLM_code/NHANES/result/'
load(paste0(data.path0,'PA_quantile_final.Rdata'))


link = 'binomial'
cv.type = 'AIC'
gamma = 0


# basic parameters
proper.indices = c(40)
lambda.list = seq(0.001,0.1,length.out=100)
Xdim.max.list = 2:5
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
  Xdata[['spaces']] = c(rep('Euclid',p2),'simplex','Wasserstein')
  
  Y = as.matrix(Y)
  colnames(Y) = NULL
  
  return(list(Y=Y,Xdata=Xdata))
}



# model fitting
start.time=Sys.time()

print(Ycols)
for (Ycol in Ycols){
  tmp = split.Y.Xdata(PA,Covariate,X.simplex,Ycol)
  Y = tmp$Y
  Xdata = tmp$Xdata
  
  for (penalty in penalty.list){
    model = GLM.CV(Xdata,Y,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list)
    save(model,file=paste0(result.path0,'NHANES_quantile_',Ycol,'_',penalty,'.Rdata'))
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('Ycol: ',Ycol,', penalty: ',penalty,', runtime: ',runtime))
  }
  
  # only PM
  model = GLM.oracle.CV(Xdata,Y,link,proper.indices,cv.type,Xdim.max.list)
  save(model,file=paste0(result.path0,'NHANES_density_',Ycol,'_only PM.Rdata'))
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('Ycol: ',Ycol,'penalty: only wind, runtime: ',runtime))
}






