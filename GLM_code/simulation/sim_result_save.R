### simulation results


path0 = './GLM_code/simulation/sim_result/'

# basic parameters
link.list = c('binomial','poisson')
cv.type.list = c('AIC','BIC')
penalty.list = c('lasso','scad','mcp','oracle')

p.list = c(100,200)


# define score
score = 'RMSE'
score.list = c(score,'CV.iteration','runtime.opt.second')
vs.score = c('CS','ICS')


# save results as csv files
loops = expand.grid(link=link.list,cv.type=cv.type.list)

for (idx in 1:nrow(loops)){
  link = loops[idx,'link']
  cv.type = loops[idx,'cv.type']
  
  print(paste(link,cv.type))
  
  if (link=='binomial'){
    s.list = c(3,6)
  }else if (link=='poisson'){
    s.list = c(2,4)
  }
  loop2 = expand.grid(s=s.list,p=p.list)
  
  # score results
  for (score in score.list){
    result.list = list()
    for (i in 1:nrow(loop2)){
      p = loop2[i,'p']
      s = loop2[i,'s']
      
      path = paste0(path0,link,'_',cv.type,'_',p,'_',s,'/')
      load(paste0(path,'sim_info.RData'))
      
      result.score = list()
      for (j in 1:length(penalty.list)){
        model = penalty.list[j]
        file.name = paste0(path,'result_',model,'.csv')
        result.each = read.csv(file.name)[-1,-1]
        
        # change the order of 'n' and 'beta.norm'
        result.each = result.each[,c(2,1,3:ncol(result.each))]
        result.each = result.each[order(result.each[['beta.norm']],result.each[['n']]),]
        
        tmp = result.each[score]
        rownames(tmp) = NULL
        colnames(tmp) = c(model)
        result.score[[j+2]] = tmp
      }
      param = rep(paste0('(',sim.info$p,',',sim.info$s,')'),nrow(result.each))
      param2 = result.each[1:2]
      rownames(param2) = NULL
      
      result.score[[1]] = data.frame(param)
      result.score[[2]] = param2
      
      n.row = min(sapply(result.score,nrow))
      result.score = lapply(result.score,function(x){x[1:n.row,,drop=FALSE]})
      
      result.list[[i]] = do.call(cbind,result.score)
    }
    result.list1 = do.call(rbind,result.list)
    write.csv(result.list1,paste0(path0,link,'_',cv.type,'_',score,'.csv'))
  }
  
  
  
  # variable selection results (CS, ICS)
  result.list = list()
  for (i in 1:nrow(loop2)){
    p = loop2[i,'p']
    s = loop2[i,'s']
    
    path = paste0(path0,link,'_',cv.type,'_',p,'_',s,'/')
    load(paste0(path,'sim_info.RData'))
    
    result.score = list()
    for (j in 1:length(penalty.list)){
      model = penalty.list[j]
      file.name = paste0(path,'result_',model,'.csv')
      result.each = read.csv(file.name)[-1,-1]
      
      # change the order of 'n' and 'beta.norm'
      result.each = result.each[,c(2,1,3:ncol(result.each))]
      result.each = result.each[order(result.each[['beta.norm']],result.each[['n']]),]
      
      tmp = result.each[vs.score]
      rownames(tmp) = NULL
      colnames(tmp) = as.vector(sapply(vs.score,FUN=function(x){paste(model,x)}))
      result.score[[j+2]] = tmp
    }
    param = rep(paste0('(',sim.info$p,',',sim.info$s,')'),nrow(result.each))
    param2 = result.each[1:2]
    rownames(param2) = NULL
    
    result.score[[1]] = data.frame(param)
    result.score[[2]] = param2
    
    n.row = min(sapply(result.score,nrow))
    result.score = lapply(result.score,function(x){x[1:n.row,,drop=FALSE]})
    
    result.list[[i]] = do.call(cbind,result.score)
  }
  result.list2 = do.call(rbind,result.list)
  write.csv(result.list2,paste0(path0,link,'_',cv.type,'_CS_ICS.csv'))
  
  
  if (link=='binomial'){
    score.class = 'accuracy'
    
    result.list = list()
    for (i in 1:nrow(loop2)){
      p = loop2[i,'p']
      s = loop2[i,'s']
      
      path = paste0(path0,link,'_',cv.type,'_',p,'_',s,'/')
      load(paste0(path,'sim_info.RData'))
      
      result.score = list()
      for (j in 1:length(penalty.list)){
        model = penalty.list[j]
        file.name = paste0(path,'result_',model,'.csv')
        result.each = read.csv(file.name)[-1,-1]
        
        # change the order of 'n' and 'beta.norm'
        result.each = result.each[,c(2,1,3:ncol(result.each))]
        result.each = result.each[order(result.each[['beta.norm']],result.each[['n']]),]
        
        tmp = result.each[score.class]
        rownames(tmp) = NULL
        colnames(tmp) = c(model)
        result.score[[j+2]] = tmp
      }
      tmp = result.each['accuracy.true']
      rownames(tmp) = NULL
      colnames(tmp) = c('true')
      result.score[[length(penalty.list)+3]] = tmp
      param = rep(paste0('(',sim.info$p,',',sim.info$s,')'),nrow(result.each))
      param2 = result.each[1:2]
      rownames(param2) = NULL
      
      result.score[[1]] = data.frame(param)
      result.score[[2]] = param2
      
      n.row = min(sapply(result.score,nrow))
      result.score = lapply(result.score,function(x){x[1:n.row,,drop=FALSE]})
      
      result.list[[i]] = do.call(cbind,result.score)
    }
    result.list3 = do.call(rbind,result.list)
    write.csv(result.list3,paste0(path0,link,'_',cv.type,'_','accuracy','.csv'))
    
    
    
    
  }
}





