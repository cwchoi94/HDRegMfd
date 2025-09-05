

path0 = './GLM_code/simulation/sim_result/'

link = 'binomial'
link = 'poisson'

cv.type = 'AIC'
# cv.type = 'BIC'


score = 'RMSE'
# score = 'CS_ICS'
score = 'class_error'



# AIC vs BIC
print(paste(link,score))
read.csv(paste0(path0,link,'_','AIC','_',score,'.csv'))
read.csv(paste0(path0,link,'_','BIC','_',score,'.csv'))


