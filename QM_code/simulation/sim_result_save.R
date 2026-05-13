### simulation results


path0 = './QM_code/simulation/sim_result/'

# basic parameters
model.list = c('QM','LM')

tau.list = c(0.1,0.5,0.9)
n.list = c(200,400)
beta.norm.list = c(0.5,1.0)
error.type.list = c('normal-homoscedastic','t-homoscedastic','normal-heteroscedastic','t-heteroscedastic')
penalty.list = c('LASSO','SCAD','MCP','ORACLE')

kernel = 'Gaussian'
cv.type = 'BIC'

p.list = c(200,400)
s.list = c(3,6)


# define score
score = 'RMSE'
score = 'MAD'
score.list = c('TP','FP','MSE','MAD','runtime.opt.second')
result.cols = c('(p,n)','MODEL',score.list)


# save results as csv files
loop1 = expand.grid(tau=tau.list,error.type=error.type.list,s=s.list)
loop2 = expand.grid(n=n.list,p=p.list)
loop3 = expand.grid(penalty=penalty.list,model=model.list)[-8,]

for (idx in 1:nrow(loop1)){
  tau = loop1[idx,'tau']
  error.type = loop1[idx,'error.type']
  s = loop1[idx,'s']
  
  print(paste(tau,error.type))
  
  # score results
  result.list.all = list()
  result.list = list()
  for (i in 1:nrow(loop2)){
    n = loop2[i,'n']
    p = loop2[i,'p']
    
    path = paste0(path0,kernel,'_',cv.type,'_',p,'_',s,'/')
    
    result.score = list()
    for (j in 1:nrow(loop3)){
      model = loop3[j,'model']
      penalty = loop3[j,'penalty']
      
      file.name = paste0(path,'result_',model,'_',error.type,'_',penalty,'.csv')
      result.each = read.csv(file.name)[-1,-1]
      rownames(result.each) = NULL
      
      result.each = result.each[(result.each$tau==tau) & (result.each$n==n),]
      tmp = result.each[score.list] 
      
      if ((model=='QM') && (penalty=='ORACLE')){
        tmp['TP'] = '-'
        tmp['FP'] = '-'
      }
      
      if (!(model=='LM') || (tau==0.5)){
        tmp['(p,n)'] = paste0('(',p,',',n,')')
        tmp['MODEL'] = paste0(model,'-',penalty)
        result.score[[j]] = tmp
      }
    }
    
    result.score = do.call(rbind,result.score)
    result.score = result.score[result.cols]
    rownames(result.score) = NULL
    
    result.list[[i]] = result.score
  }
  result.list.all = do.call(rbind,result.list)
  write.csv(result.list.all,paste0(path0,kernel,'_',cv.type,'_',tau,'_',s,'_',error.type,'.csv'))
}





