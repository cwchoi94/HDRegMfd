### simulation results


model.list = c('lasso','scad','mcp','oracle')

name = 'sim_BayesHilbert'
sim.list = c('sim1','sim2','sim3','sim4')

# name = 'sim_SPD_Affine'
# sim.list = c('sim5','sim6','sim7','sim8')
 
# name = 'sim_SPD_LogEuclid'
# sim.list = c('sim9','sim10','sim11','sim12')


path = './LM_code/simulation/sim_result/'


score = 'MSPE'

result.list = list()
for (i in 1:length(sim.list)){
  sim = sim.list[i]
  save_path = paste0(path,sim,'/')
  load(paste0(save_path,'sim_info.RData'))
  
  result.score = list()
  for (j in 1:length(model.list)){
    model = model.list[j]
    file.name = paste0(save_path,'result_',model,'.csv')
    result.each = read.csv(file.name)[-1,-1]
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
  
  result.list[[i]] = do.call(cbind,result.score)
}

result.list1 = do.call(rbind,result.list)
write.csv(result.list1,paste0(path,name,'_',score,'.csv'))


score = c('CS','ICS')

result.list = list()
for (i in 1:length(sim.list)){
  sim = sim.list[i]
  save_path = paste0(path,sim,'/')
  load(paste0(save_path,'sim_info.RData'))
  
  result.score = list()
  for (j in 1:(length(model.list)-1)){
    model = model.list[j]
    file.name = paste0(save_path,'result_',model,'.csv')
    result.each = read.csv(file.name)[-1,-1]
    tmp = result.each[score]
    rownames(tmp) = NULL
    colnames(tmp) = as.vector(sapply(score,FUN=function(x){paste(model,x)}))
    result.score[[j+2]] = tmp
  }
  param = rep(paste0('(',sim.info$p,',',sim.info$s,')'),nrow(result.each))
  param2 = result.each[1:2]
  rownames(param2) = NULL
  
  result.score[[1]] = data.frame(param)
  result.score[[2]] = param2
  
  result.list[[i]] = do.call(cbind,result.score)
}

result.list2 = do.call(rbind,result.list)
write.csv(result.list2,paste0(path,name,'_','CS_ICS.csv'))

result.list1
result.list2

