
library(HDRegMfd)
library(hms)


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))


link = 'binomial'
cv.type = 'AIC'
gamma = 0

dim(Y)
Xdata[['spaces']]
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# basic parameters
proper.indices = c(12) # wind
lambda.list = seq(0.002,0.04,length.out=20)
Xdim.max.list = 2:5
R.list = c(100,200)
penalty.list = c('LASSO','SCAD','MCP')


# model fitting
start.time=Sys.time()

for (penalty in penalty.list){
  model = GLM.CV(Xdata,Y,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list)
  save(model,file=paste0(result.path0,'PM2.5_',penalty,'.Rdata'))
  
  runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
  print(paste0('penalty: ',penalty,', runtime: ',runtime))
}

# only wind
model = GLM.oracle.CV(Xdata,Y,link,proper.indices,cv.type,Xdim.max.list)
save(model,file=paste0(result.path0,'PM2.5_only wind.Rdata'))

runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
print(paste0('penalty: only wind, runtime: ',runtime))





