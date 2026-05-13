
library(HDRegMfd)
library(hms)


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

dim(Y)
Xdata[['spaces']]
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# basic parameters
proper.indices = c(12) # wind
lambda.list = seq(0.0005,0.01,length.out=20)
Xdim.max.list = 1:5
penalty.list = c('LASSO','SCAD','MCP')


# model fitting
start.time=Sys.time()

for (tau in tau.list){
  for (penalty in penalty.list){
    model = QM.CV(Xdata,Y,tau,h,kernel,cv.type,penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const)
    save(model,file=paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))
    
    runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
    print(paste0('tau: ',tau,', penalty: ',penalty,', runtime: ',runtime))
  }
  
  # only wind
  model = QM.oracle.CV(Xdata,Y,tau,h,kernel,proper.indices,cv.type,Xdim.max.list,cv.const=cv.const)
  save(model,file=paste0(result.path0,'PM2.5_',tau,'_only wind.Rdata'))
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('tau: ',tau,', penalty: only wind, runtime: ',runtime))
}










