
library(HDRegMfd)
library(hms)


# must include normalizer=(Y.max-Y.min)^(-1) in "inv.clr.density" function.
data.path0 = './LM_code/BikeRental/data/'
result.path0 = './LM_code/BikeRental/result/'
load(paste0(data.path0,'BikeRental.Rdata'))

Yspace = 'functional'

dim(Y)
Xdata[['spaces']]
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# basic parameters
kfold = 5
seed = 10000
gamma = 0

lambda.list = seq(0.0005,0.03,length.out=60)
Xdim.max.list = 1:10
R.list = seq(5,10,length.out=2)
penalty.list = c('LASSO','SCAD','MCP')


# model fitting
start.time=Sys.time()

for (penalty in penalty.list){
  model = LM.kfold(Xdata,Y,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list)
  save(model,file=paste0(result.path0,'BikeRental_',penalty,'.Rdata'))
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('penalty: ',penalty,', runtime: ',runtime))
}


