### simulation source


library(HDRegMfd)
library(hms)


# set basic simulation settings
n.iteration = 100


if (base.model=='QM'){
  tau.list = c(0.1,0.5,0.9)
}else if (base.model=='LM'){
  tau.list = c(0.5)
}

# basic parameters
# error.type.idx: index of error.type, 1=='normal-homoscedastic', 2=='normal-heteroscedastic', 3=='t-homoscedastic' and 4=='t-heterogeneous'
error.type.idx.list = c(1,2,3,4) 
n.list = c(200,400)
beta.norm.range = c(1.0,1.5)


# define common parameters
if (base.model=='QM'){
  lambda.list = seq(0.01,0.1,length.out=19)
  Xdim.max.list = 1:6
}else if (base.model=='LM'){
  lambda.list = seq(0.06,0.5,length.out=45)
  Xdim.max.list = 1:6
  R.list = c(10000,20000)
}


cv.type = 'AIC'
cv.type = 'BIC'
h = NULL

kernel = 'Gaussian'
# kernel = 'Uniform'
# kernel = 'Epanechnikov'
# kernel = 'Triangular'


Yspace = 'Euclid'

beta0.norm = 1
Xrho = 0.5
Xsigma = 1
error.sigma = 1
df = 2

cv.const = 2
gamma = 0
c.beta = NULL





# define additional parameters

get.dim.each = function(x){
  if (x=='Euclid'){
    dim = 1
  } else if (x=='simplex'){
    dim = 3
  } else if (x=='sphere'){
    dim = 3
  } else if (x=='functional'){
    dim = 50
  } else if (x=='BayesHilbert'){
    dim = 40
  } else if (x=='Wasserstein'){
    dim = 5
  }
  return(dim)
}

Xspaces = c(rep('functional',q/2),rep('Wasserstein',q/2),rep('Euclid',(p-q)/2),rep('simplex',(p-q)/4),rep('sphere',(p-q)/4))
Xdims = sapply(Xspaces,get.dim.each,USE.NAMES=FALSE)
proper.indices = c(2,if(s>1){sapply(1:(s-1),function(j){q+floor((p-q)/(s-1)*j)})}else{})
s = length(proper.indices)


# create save folder
save_path = paste0('./sim_result/',kernel,'_',cv.type,'_',p,'_',s,'/')
save_rdata_path = paste0(save_path,'Rdata/')
dir.create('./sim_result/',showWarnings=FALSE)
dir.create(save_path,showWarnings=FALSE)
dir.create(paste0(save_path,'result_each_',base.model,'/'),showWarnings=FALSE)


model.names = c('LASSO','SCAD','MCP','ORACLE')
scores = c('TP','FP','RMSE','MSE','MAD','RMSPE','MSPE','runtime.second','runtime.opt.second')
cols.each = c('iter',scores,'opt.lambda','opt.Xdim.max')
cols = c('tau','n',scores)


# save parameters
# sim.info: 1 dim parameters
# parameters: list parameters
sim.info = list(p=p,q=q,s=s,Xrho=Xrho,Xsigma=Xsigma)
sim.info = data.frame(sim.info)
parameters = list(Xspaces=Xspaces,Xdims=Xdims,proper.indices=proper.indices,
                  lambda.list=lambda.list,Xdim.max.list=Xdim.max.list,
                  n.list=n.list)

save(sim.info,file=paste0(save_path,'sim_info.RData'))
save(parameters,file=paste0(save_path,'parameters.RData'))


# simulation
start.time = Sys.time()

message()

loops = expand.grid(n=n.list,tau=tau.list)
n.count = nrow(loops)
n.count.all = n.count * length(error.type.idx.list)

count.all = 0
for (error.type.idx in error.type.idx.list){
  
  # set error distribution
  if (error.type.idx==1){
    error.type = 'normal'
    error.var.type = 'homoscedastic'
  }else if (error.type.idx==2){
    error.type = 'normal'
    error.var.type = 'heteroscedastic'
  }else if (error.type.idx==3){
    error.type = 't'
    error.var.type = 'homoscedastic'
  }else if (error.type.idx==4){
    error.type = 't'
    error.var.type = 'heteroscedastic'
  }
  error.type.org = paste0(error.type,'-',error.var.type)
  
  # make folders
  result.name.base.each = paste0(save_path,'result_each_',base.model,'/result_',error.type.org,'_')
  result.name.base = paste0(save_path,'result_',base.model,'_',error.type.org,'_')
  dir.create(save_rdata_path,showWarnings=FALSE)
  
  
  # make all result csv file
  for (name in model.names){
    assign(paste0('result.',name),t(data.frame(rep(0,length(cols)),row.names=cols)))
  }
  
  for (count in 1:n.count){
    count.all = count.all + 1
    if (count.all<start.count){next}
    
    tau = loops[count,'tau']
    n = loops[count,'n']
    
    
    if ((base.model=='LM') & (error.type=='t')){
      lambda.list = seq(0.2,0.5,length.out=31)
    }else if ((base.model=='LM') & (error.type=='normal')){
      lambda.list = seq(0.06,0.5,length.out=45)
    }
    
    
    # make each result csv file if it doesn't exist
    for (name in model.names){
      if (!file.exists(paste0(result.name.base.each,name,'_',tau,'_',n,'.csv'))){
        result.each = t(data.frame(rep(0,length(cols.each)),row.names=cols.each))
        rownames(result.each) = NULL
        write.csv(result.each,file=paste0(result.name.base.each,name,'_',tau,'_',n,'.csv'))
      }
    }
    
    # generate beta.norm
    set.seed(13579+count.all)
    beta.norm = rep(0,p)
    beta.norm[proper.indices] = runif(s,beta.norm.range[1],beta.norm.range[2])
    
    
    for (iter in 1:n.iteration){
      if ((count.all==start.count) & (iter<start.iteration)){next}
      
      # check if the simulation for 'iter'th iteration was implemented
      model.run.check = sapply(model.names,function(name){
        result.each = read.csv(paste0(result.name.base.each,name,'_',tau,'_',n,'.csv'))[,-1]
        return(iter %in% result.each[,'iter'])
      })
      if (all(model.run.check)){next}
      
      message('model: ',base.model,', cv.type: ',cv.type, ', p: ',p,', s: ',s)
      message('count: ',count.all,'/',n.count.all,', error.type: ',error.type.org,', tau: ',tau,', n: ',n,', iter: ',iter,'/',n.iteration)
      
      # data generate
      if (is.null(c.beta)){
        data = QM.data.generate(n,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+1000)
        datatest = QM.data.generate(n2,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+2000,c.beta=data$c.beta)
        c.beta = data$c.beta
      }else{
        data = QM.data.generate(n,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+1000,c.beta=c.beta)
        datatest = QM.data.generate(n2,Xspaces,Xdims,tau,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,error.sigma,error.var.type,error.type,df,seed=iter+2000,c.beta=data$c.beta)
      }
      
      # LASSO
      penalty = 'LASSO'
      result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        if (base.model=='QM'){
          model = QM.CV(data$X,data$Y,tau,h,kernel,cv.type,penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const)
        }else if (base.model=='LM'){
          model = LM.CV(data$X,data$Y,Yspace,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
        }
        
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        
        # variable selection
        TP = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        FP = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(datatest$Y,Yhat)^2)
        MSE = mean(dist.manifold(datatest$Xbeta,Yhat)^2)
        MAD = mean(dist.manifold(datatest$Xbeta,Yhat))
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        
        tmp.scores = list(iter,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,model$runtime.second,model$runtime.opt.second,model$lambda,model$Xdim.max)
        result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # SCAD
      penalty = 'SCAD'
      result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        if (base.model=='QM'){
          model = QM.CV(data$X,data$Y,tau,h,kernel,cv.type,penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const)
        }else if (base.model=='LM'){
          model = LM.CV(data$X,data$Y,Yspace,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
        }
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        
        # variable selection
        TP = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        FP = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(datatest$Y,Yhat)^2)
        MSE = mean(dist.manifold(datatest$Xbeta,Yhat)^2)
        MAD = mean(dist.manifold(datatest$Xbeta,Yhat))
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        
        tmp.scores = list(iter,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,model$runtime.second,model$runtime.opt.second,model$lambda,model$Xdim.max)
        result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # MCP
      penalty = 'MCP'
      result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        if (base.model=='QM'){
          model = QM.CV(data$X,data$Y,tau,h,kernel,cv.type,penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const)
        }else if (base.model=='LM'){
          model = LM.CV(data$X,data$Y,Yspace,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
        }
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        
        # variable selection
        TP = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        FP = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(datatest$Y,Yhat)^2)
        MSE = mean(dist.manifold(datatest$Xbeta,Yhat)^2)
        MAD = mean(dist.manifold(datatest$Xbeta,Yhat))
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        
        tmp.scores = list(iter,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,model$runtime.second,model$runtime.opt.second,model$lambda,model$Xdim.max)
        result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # oracle
      penalty = 'ORACLE'
      result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        if (base.model=='QM'){
          model = QM.oracle.CV(data$X,data$Y,tau,h,kernel,proper.indices,cv.type,Xdim.max.list,cv.const=cv.const)
        }else if (base.model=='LM'){
          model = LM.oracle.CV(data$X,data$Y,Yspace,proper.indices,cv.type,Xdim.max.list,cv.const=cv.const)
        }
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        
        # variable selection
        TP = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        FP = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(datatest$Y,Yhat)^2)
        MSE = mean(dist.manifold(datatest$Xbeta,Yhat)^2)
        MAD = mean(dist.manifold(datatest$Xbeta,Yhat))
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        
        tmp.scores = list(iter,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,model$runtime.second,model$runtime.opt.second,0,model$Xdim.max)
        result.each = read.csv(paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base.each,penalty,'_',tau,'_',n,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      message()
    }
    
    
    # Update simulation result
    # LASSO
    result.each = read.csv(paste0(result.name.base.each,'LASSO','_',tau,'_',n,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('TP','FP')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(tau,n,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,runtime.second,runtime.opt.second)
    result.LASSO = rbind(result.LASSO,result.tmp)
    rownames(result.LASSO) = NULL
    write.csv(result.LASSO,paste0(result.name.base,'LASSO.csv'))
    
    # SCAD
    result.each = read.csv(paste0(result.name.base.each,'SCAD','_',tau,'_',n,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('TP','FP')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(tau,n,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,runtime.second,runtime.opt.second)
    result.SCAD = rbind(result.SCAD,result.tmp)
    rownames(result.SCAD) = NULL
    write.csv(result.SCAD,paste0(result.name.base,'SCAD.csv'))
    
    # MCP
    result.each = read.csv(paste0(result.name.base.each,'MCP','_',tau,'_',n,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('TP','FP')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(tau,n,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,runtime.second,runtime.opt.second)
    result.MCP = rbind(result.MCP,result.tmp)
    rownames(result.MCP) = NULL
    write.csv(result.MCP,paste0(result.name.base,'MCP.csv'))
    
    # oracle
    result.each = read.csv(paste0(result.name.base.each,'ORACLE','_',tau,'_',n,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('TP','FP')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(tau,n,TP,FP,RMSE,MSE,MAD,RMSPE,MSPE,runtime.second,runtime.opt.second)
    result.ORACLE = rbind(result.ORACLE,result.tmp)
    rownames(result.ORACLE) = NULL
    write.csv(result.ORACLE,paste0(result.name.base,'ORACLE.csv'))
    
    message()
    
  }
  
  
}








