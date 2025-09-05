### simulation source
### Implemented for each link ('binomial' or 'poisson') and cv.type ('AIC' or 'BIC') 


library(HDRegMfd)
library(hms)


# define common parameters

lambda.list = seq(0.01,0.4,length.out=40)
Xdim.max.list = 1:8
R.list = seq(200,400,length.out=2)

beta0.norm = 1
Xrho = 0.5
Xsigma = 1

gamma = 0
phi = 1
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
proper.indices = c(2,q/2+2,if(s>2){sapply(1:(s-2),function(j){q+floor((p-q)/(s-2)*j)})}else{})
#proper.indices = proper.indices[proper.indices<=p]


# create save folder
save_path = paste0('./sim_result/',link,'_',cv.type,'_',p,'_',s,'/')
save_rdata_path = paste0(save_path,'Rdata/')
result.name.base = paste0(save_path,'result_')
dir.create('./sim_result/',showWarnings=FALSE)
dir.create(save_path,showWarnings=FALSE)
dir.create(save_rdata_path,showWarnings=FALSE)


# save parameters
# sim.info: 1 dim parameters
# parameters: list parameters
sim.info = list(p=p,q=q,s=s,link=link,Xrho=Xrho,Xsigma=Xsigma)
sim.info = data.frame(sim.info)
parameters = list(Xspaces=Xspaces,Xdims=Xdims,proper.indices=proper.indices,
                  lambda.list=lambda.list,Xdim.max.list=Xdim.max.list,R.list=R.list,
                  n.list=n.list,beta.norm.list=beta.norm.list)

save(sim.info,file=paste0(save_path,'sim_info.RData'))
save(parameters,file=paste0(save_path,'parameters.RData'))


# make or read all result csv file
model.names = c('lasso','scad','mcp','oracle')
scores = c('CS','ICS','RMSE','MSE','RMSE2','MSE2','CV.iteration','runtime.second','runtime.opt.second')
cols.each = c('iter',scores,'opt.lambda','opt.Xdim.max','opt.R')
cols = c('n','beta.norm',scores)

for (name in model.names){
  assign(paste0('result.',name),t(data.frame(rep(0,length(cols)),row.names=cols)))
}


# simulation
start.time = Sys.time()

message()

count = 0
n.count = length(n.list)*length(beta.norm.list)
for (n in n.list){
  for (beta.norm in beta.norm.list){
    count = count+1
    if (count<start.count){next}
    
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
      if ((count==start.count) & (iter<start.iteration)){next}
      
      # check if the simulation for 'iter'th iteration was implemented
      model.run.check = sapply(model.names,function(name){
        result.each = read.csv(paste0(result.name.base,name,'_',count,'.csv'))[,-1]
        return(iter %in% result.each[,'iter'])
      })
      if (all(model.run.check)){next}
      
      message('link: ',link,', cv.type: ',cv.type, ', p: ',p,', s: ',s)
      message('count: ',count,'/',n.count,', n: ',n,', beta.norm: ',beta.norm,', iter: ',iter,'/',n.iteration)
      
      # data generate
      if (is.null(c.beta)){
        data = GLM.data.generate(n,Xspaces,Xdims,link,Ydim,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,seed=iter+1000)
        datatest = GLM.data.generate(n3,Xspaces,Xdims,link,Ydim,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,seed=iter+2000,c.beta=data$c.beta)
        c.beta = data$c.beta
      }else{
        data = GLM.data.generate(n,Xspaces,Xdims,link,Ydim,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,seed=iter+1000,c.beta=c.beta)
        datatest = GLM.data.generate(n3,Xspaces,Xdims,link,Ydim,proper.indices,beta.norm,beta0.norm,Xrho,Xsigma,seed=iter+2000,c.beta=c.beta)
      }
      
      # LASSO
      result.each = read.csv(paste0(result.name.base,'lasso','_',count,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        penalty = 'LASSO'
        model = GLM.CV(data$X,data$Y,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        theta = predict(model,datatest$X,is.inv.link=FALSE)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSE = mean(dist.manifold(datatest$Ymu,Yhat)^2)
        MSE2 = mean(dist.manifold(datatest$theta,theta)^2)
        RMSE = sqrt(MSE)
        RMSE2 = sqrt(MSE2)
        
        # runtime
        CV.iteration = (nrow(model$parameter.list)+2)%/%3
        runtime.second = model$runtime.second
        runtime.opt.second = model$runtime.opt.second
        
        tmp.scores = list(iter,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second,model$lambda,model$Xdim.max,model$R)
        result.each = read.csv(paste0(result.name.base,'lasso','_',count,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base,'lasso','_',count,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # SCAD
      result.each = read.csv(paste0(result.name.base,'scad','_',count,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        penalty = 'SCAD'
        model = GLM.CV(data$X,data$Y,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        theta = predict(model,datatest$X,is.inv.link=FALSE)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSE = mean(dist.manifold(datatest$Ymu,Yhat)^2)
        MSE2 = mean(dist.manifold(datatest$theta,theta)^2)
        RMSE = sqrt(MSE)
        RMSE2 = sqrt(MSE2)
        
        # runtime
        CV.iteration = (nrow(model$parameter.list)+2)%/%3
        runtime.second = model$runtime.second
        runtime.opt.second = model$runtime.opt.second
        
        tmp.scores = list(iter,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second,model$lambda,model$Xdim.max,model$R)
        result.each = read.csv(paste0(result.name.base,'scad','_',count,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base,'scad','_',count,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # MCP
      result.each = read.csv(paste0(result.name.base,'mcp','_',count,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        penalty = 'MCP'
        model = GLM.CV(data$X,data$Y,link,cv.type,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        theta = predict(model,datatest$X,is.inv.link=FALSE)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSE = mean(dist.manifold(datatest$Ymu,Yhat)^2)
        MSE2 = mean(dist.manifold(datatest$theta,theta)^2)
        RMSE = sqrt(MSE)
        RMSE2 = sqrt(MSE2)
        
        # runtime
        CV.iteration = (nrow(model$parameter.list)+2)%/%3
        runtime.second = model$runtime.second
        runtime.opt.second = model$runtime.opt.second
        
        tmp.scores = list(iter,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second,model$lambda,model$Xdim.max,model$R)
        result.each = read.csv(paste0(result.name.base,'mcp','_',count,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base,'mcp','_',count,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      
      # oracle
      result.each = read.csv(paste0(result.name.base,'oracle','_',count,'.csv'))[,-1]
      if (!(iter %in% result.each[,'iter'])){
        penalty = 'oracle'
        model = GLM.oracle.CV(data$X,data$Y,link,proper.indices,cv.type,Xdim.max.list)
        
        # predict Ytest
        Yhat = predict(model,datatest$X)
        theta = predict(model,datatest$X,is.inv.link=FALSE)
        
        # variable selection
        CS = length(model$proper.indices[which(model$proper.indices %in% proper.indices)])
        ICS = length(model$proper.indices[-which(model$proper.indices %in% proper.indices)])
        
        # prediction
        MSE = mean(dist.manifold(datatest$Ymu,Yhat)^2)
        MSE2 = mean(dist.manifold(datatest$theta,theta)^2)
        RMSE = sqrt(MSE)
        RMSE2 = sqrt(MSE2)
        
        # runtime
        CV.iteration = 0
        runtime.second = model$runtime.second
        runtime.opt.second = model$runtime.opt.second
        
        tmp.scores = list(iter,CS,ICS,RMSE,MSE,RMSE2,MSE2,0,runtime.second,runtime.opt.second,0,model$Xdim.max,0)
        result.each = read.csv(paste0(result.name.base,'oracle','_',count,'.csv'))[,-1]
        result.each = rbind(result.each,tmp.scores)
        result.each = result.each[order(result.each[,'iter']),]
        result.each = result.each[!duplicated(result.each[,'iter']),]
        rownames(result.each) = NULL
        write.csv(result.each,paste0(result.name.base,'oracle','_',count,'.csv'))
        
        runtime = hms(round(as.numeric(difftime(Sys.time(),start.time,units='secs'))))
        message('iter: ',iter,'/',n.iteration,', ',penalty,', runtime: ',runtime)
      }
      
      message()
    }
    
    
    # Update simulation result
    
    # LASSO
    result.each = read.csv(paste0(result.name.base,'lasso','_',count,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second)
    result.lasso = rbind(result.lasso,result.tmp)
    rownames(result.lasso) = NULL
    write.csv(result.lasso,paste0(result.name.base,'lasso.csv'))
    
    # SCAD
    result.each = read.csv(paste0(result.name.base,'scad','_',count,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second)
    result.scad = rbind(result.scad,result.tmp)
    rownames(result.scad) = NULL
    write.csv(result.scad,paste0(result.name.base,'scad.csv'))
    
    # MCP
    result.each = read.csv(paste0(result.name.base,'mcp','_',count,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second)
    result.mcp = rbind(result.mcp,result.tmp)
    rownames(result.mcp) = NULL
    write.csv(result.scad,paste0(result.name.base,'mcp.csv'))
    
    # oracle
    result.each = read.csv(paste0(result.name.base,'oracle','_',count,'.csv'))[-1,]
    result.each = result.each[rowSums(is.na(result.each))==0,]
    for (score in scores){
      score.mean = mean(result.each[[score]])
      score.se = sd(result.each[[score]])/sqrt(nrow(result.each))
      if (score %in% c('CS','ICS')){
        assign(score,paste0(sprintf('%0.2f',score.mean),'(',sprintf('%0.2f',score.se),')'))
      } else{
        assign(score,paste0(sprintf('%0.3f',score.mean),'(',sprintf('%0.3f',score.se),')'))
      }
    }
    result.tmp = list(n,beta.norm,CS,ICS,RMSE,MSE,RMSE2,MSE2,CV.iteration,runtime.second,runtime.opt.second)
    result.oracle = rbind(result.oracle,result.tmp)
    rownames(result.oracle) = NULL
    write.csv(result.scad,paste0(result.name.base,'oracle.csv'))
    
    message()
  }
}







