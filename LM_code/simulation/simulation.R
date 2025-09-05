### simulation source
### Implemented for both GCV and Kfold


library(HDRegMfd)
library(hms)


# define common parameters
lambda.list = seq(0.005,0.2,length.out=40)
Xdim.max.list = 1:8
R.list = seq(100,200,length.out=2)

n.iteration = 100

Xrho = 0.5
Xsigma = 1
error.rho = 0.5
error.sigma = 1

gamma = 0
phi = 1
cv.type = 'GCV'


# define basic parameters

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
proper.indices = c(c(2,q/2+2),sapply(1:(s-2),function(j){q+floor((p-q)/(s-2)*j)}))
#proper.indices = proper.indices[proper.indices<=p]


# create save folder
save_path = paste0('./sim_result/',sim,'/')
save_rdata_path = paste0(save_path,'Rdata/')
dir.create('./sim_result/',showWarnings=FALSE)
dir.create(save_path,showWarnings=FALSE)
dir.create(save_rdata_path,showWarnings=FALSE)


# save parameters
# sim.info: 1 dim parameters
# parameters: list parameters
sim.info = list(p=p,q=q,s=s,Xrho=Xrho,Xsigma=Xsigma,error.rho=error.rho,error.sigma=error.sigma,Yspace=Yspace)
sim.info = data.frame(sim.info)
parameters = list(Xspaces=Xspaces,Xdims=Xdims,Ydim=Ydim,proper.indices=proper.indices,
                  lambda.list=lambda.list,Xdim.max.list=Xdim.max.list,R.list=R.list,
                  n.list=n.list,beta.norm.list=beta.norm.list)

save(sim.info,file=paste0(save_path,'sim_info.RData'))
save(parameters,file=paste0(save_path,'parameters.RData'))


# make or read all result csv file
model.names = c('lasso','scad','mcp','oracle')
scores = c('CS','ICS','RMSPE','MSPE','RMSE','MSE','RMSPE2','MSPE2','RMSE2','MSE2')
cols.each = c('iter',scores,'opt.lambda','opt.Xdim.max','opt.R')
cols = c('n','beta.norm',scores)

for (name in model.names){
  assign(paste0('result.',name),t(data.frame(rep(0,length(cols)),row.names=cols)))
}


# simulation
start.time = Sys.time()

count = 0
n.count = length(n.list)*length(beta.norm.list)
for (n in n.list){
  for (beta.norm in beta.norm.list){
    count = count+1
    # if (count<=2){next}
    # if (count<=4){next}
    
    # make or read each result csv file
    for (name in model.names){
      if (file.exists(paste0(save_path,'result_',name,'_',count,'.csv'))){
        assign(paste0('result.',name,'.each'),read.csv(paste0(save_path,'result_',name,'_',count,'.csv'))[,-1])
      } else{
        assign(paste0('result.',name,'.each'),t(data.frame(rep(0,length(cols.each)),row.names=cols.each)))
        eval(parse(text=paste0('rownames(result.',name,'.each) = NULL')))
        eval(parse(text=paste0('write.csv(result.',name,'.each,file="',save_path,'result_',name,'_',count,'.csv")')))
      }
    }
    
    for (iter in 1:n.iteration){
      if (iter+1<=min(nrow(result.lasso.each),nrow(result.scad.each),nrow(result.mcp.each),nrow(result.oracle.each))){next}
      
      message('count: ',count,'/',n.count,', n: ',n,', beta.norm: ',beta.norm,', iter: ',iter,'/',n.iteration)
      
      # data generate
      data = LM.data.generate(n,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=iter+1000)
      datanew = LM.data.generate(n2,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=iter+2000)
      datatest = LM.data.generate(n3,Xspaces,Yspace,Xdims,Ydim,proper.indices,beta.norm,Xrho,Xsigma,error.rho,error.sigma,seed=iter+3000)
      
      # LASSO
      if (iter+1>nrow(result.lasso.each)){
        penalty = 'LASSO'
        model = LM.GCV(data$X,data$Y,datanew$X,datanew$Y,Yspace,penalty,gamma,lambda.list,Xdim.max.list,R.list)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        LogYhat = RieLog.manifold(model$Ymu,Yhat,Yspace)
        
        LogYtest = RieLog.manifold(model$Ymu,datatest$Y,Yspace)
        LogExpXbeta = RieLog.manifold(model$Ymu,datatest$ExpXbeta,Yspace)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(Yhat,datatest$Y,Yspace)^2)
        MSE = mean(dist.manifold(Yhat,datatest$ExpXbeta,Yspace)^2)
        MSPE2 = mean(norm.manifold(LogYhat-LogYtest,model$Ymu,Yspace)^2)
        MSE2 = mean(norm.manifold(LogYhat-LogExpXbeta,model$Ymu,Yspace)^2)
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        RMSPE2 = sqrt(MSPE2)
        RMSE2 = sqrt(MSE2)
        
        tmp.scores = list(iter,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2,model$lambda,model$Xdim.max,model$R)
        result.lasso.each = rbind(result.lasso.each,tmp.scores)
        write.csv(result.lasso.each,paste0(save_path,'result_lasso_',count,'.csv'))
        save(model,file=paste0(save_rdata_path,'lasso_',count,'_',iter,'.RData'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # SCAD
      if (iter+1>nrow(result.scad.each)){
        penalty = 'SCAD'
        model = LM.GCV(data$X,data$Y,datanew$X,datanew$Y,Yspace,penalty,gamma,lambda.list,Xdim.max.list,R.list)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        LogYhat = RieLog.manifold(model$Ymu,Yhat,Yspace)
        
        LogYtest = RieLog.manifold(model$Ymu,datatest$Y,Yspace)
        LogExpXbeta = RieLog.manifold(model$Ymu,datatest$ExpXbeta,Yspace)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(Yhat,datatest$Y,Yspace)^2)
        MSE = mean(dist.manifold(Yhat,datatest$ExpXbeta,Yspace)^2)
        MSPE2 = mean(norm.manifold(LogYhat-LogYtest,model$Ymu,Yspace)^2)
        MSE2 = mean(norm.manifold(LogYhat-LogExpXbeta,model$Ymu,Yspace)^2)
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        RMSPE2 = sqrt(MSPE2)
        RMSE2 = sqrt(MSE2)
        
        tmp.scores = list(iter,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2,model$lambda,model$Xdim.max,model$R)
        result.scad.each = rbind(result.scad.each,tmp.scores)
        write.csv(result.scad.each,paste0(save_path,'result_scad_',count,'.csv'))
        save(model,file=paste0(save_rdata_path,'scad_',count,'_',iter,'.RData'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # MCP
      if (iter+1>nrow(result.mcp.each)){
        penalty = 'MCP'
        model = LM.GCV(data$X,data$Y,datanew$X,datanew$Y,Yspace,penalty,gamma,lambda.list,Xdim.max.list,R.list)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        LogYhat = RieLog.manifold(model$Ymu,Yhat,Yspace)
        
        LogYtest = RieLog.manifold(model$Ymu,datatest$Y,Yspace)
        LogExpXbeta = RieLog.manifold(model$Ymu,datatest$ExpXbeta,Yspace)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(Yhat,datatest$Y,Yspace)^2)
        MSE = mean(dist.manifold(Yhat,datatest$ExpXbeta,Yspace)^2)
        MSPE2 = mean(norm.manifold(LogYhat-LogYtest,model$Ymu,Yspace)^2)
        MSE2 = mean(norm.manifold(LogYhat-LogExpXbeta,model$Ymu,Yspace)^2)
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        RMSPE2 = sqrt(MSPE2)
        RMSE2 = sqrt(MSE2)
        
        tmp.scores = list(iter,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2,model$lambda,model$Xdim.max,model$R)
        result.mcp.each = rbind(result.mcp.each,tmp.scores)
        write.csv(result.mcp.each,paste0(save_path,'result_mcp_',count,'.csv'))
        save(model,file=paste0(save_rdata_path,'mcp_',count,'_',iter,'.RData'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # oracle
      if (iter+1>nrow(result.oracle.each)){
        penalty = 'oracle'
        model = LM.oracle.GCV(data$X,data$Y,datanew$X,datanew$Y,Yspace,proper.indices,Xdim.max.list)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        LogYhat = RieLog.manifold(model$Ymu,Yhat,Yspace)
        
        LogYtest = RieLog.manifold(model$Ymu,datatest$Y,Yspace)
        LogExpXbeta = RieLog.manifold(model$Ymu,datatest$ExpXbeta,Yspace)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSPE = mean(dist.manifold(Yhat,datatest$Y,Yspace)^2)
        MSE = mean(dist.manifold(Yhat,datatest$ExpXbeta,Yspace)^2)
        MSPE2 = mean(norm.manifold(LogYhat-LogYtest,model$Ymu,Yspace)^2)
        MSE2 = mean(norm.manifold(LogYhat-LogExpXbeta,model$Ymu,Yspace)^2)
        RMSPE = sqrt(MSPE)
        RMSE = sqrt(MSE)
        RMSPE2 = sqrt(MSPE2)
        RMSE2 = sqrt(MSE2)
        
        tmp.scores = list(iter,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2,0,model$Xdim.max,0)
        result.oracle.each = rbind(result.oracle.each,tmp.scores)
        write.csv(result.oracle.each,paste0(save_path,'result_oracle_',count,'.csv'))
        save(model,file=paste0(save_rdata_path,'oracle_',count,'_',iter,'.RData'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      message()
    }
    
    
    # Update simulation result
    
    # LASSO
    result.each = read.csv(paste0(save_path,'result_lasso_',count,'.csv'))[-1,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(n.iteration)
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2)
    result.lasso = rbind(result.lasso,result.tmp)
    write.csv(result.lasso,paste0(save_path,'result_lasso.csv'))
    
    # SCAD
    result.each = read.csv(paste0(save_path,'result_scad_',count,'.csv'))[-1,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(n.iteration)
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2)
    result.scad = rbind(result.scad,result.tmp)
    write.csv(result.scad,paste0(save_path,'result_scad.csv'))
    
    # MCP
    result.each = read.csv(paste0(save_path,'result_mcp_',count,'.csv'))[-1,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(n.iteration)
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2)
    result.mcp = rbind(result.mcp,result.tmp)
    write.csv(result.mcp,paste0(save_path,'result_mcp.csv'))
    
    # oracle
    result.each = read.csv(paste0(save_path,'result_oracle_',count,'.csv'))[-1,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(n.iteration)
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSPE,MSPE,RMSE,MSE,RMSPE2,MSPE2,RMSE2,MSE2)
    result.oracle = rbind(result.oracle,result.tmp)
    write.csv(result.oracle,paste0(save_path,'result_oracle.csv'))
    
    message()
  }
}







