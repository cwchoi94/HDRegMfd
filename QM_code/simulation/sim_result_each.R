### simulation results


kernel = 'Gaussian'

# cv.type = 'AIC'
cv.type = 'BIC'

p = 100
s = 3

save_path = paste0('./QM_code/simulation/sim_result/',kernel,'_',cv.type,'_',p,'_',s,'/')
save_rdata_path = paste0(save_path,'Rdata/')

load(paste0(save_path,'sim_info.RData'))
load(paste0(save_path,'parameters.RData'))


parameters$Xspaces
parameters$Xdims

# basic information

sim.info
parameters$proper.indices


# total result
base.model.name = 'QM'
# base.model.name = 'LM'

error.type.list = c('normal-homoscedastic','t-homoscedastic','normal-heteroscedastic','t-heteroscedastic')
error.type = error.type.list[1]

read.csv(paste0(save_path,'result_',base.model.name,'_',error.type,'_lasso.csv'))[-1,-1]
read.csv(paste0(save_path,'result_',base.model.name,'_',error.type,'_scad.csv'))[-1,-1]
read.csv(paste0(save_path,'result_',base.model.name,'_',error.type,'_mcp.csv'))[-1,-1]
read.csv(paste0(save_path,'result_',base.model.name,'_',error.type,'_oracle.csv'))[-1,-1]


# each result
tau.list = c(0.1,0.5,0.9)
penalty.list = c('LASSO','SCAD','MCP','ORACLE')
n.list = c(200,400)

tau = tau.list[1]
penalty = penalty.list[1]
n = n.list[1]
error.type = error.type.list[1]

save_path_each = paste0(save_path,'result_each_',base.model.name,'/')
read.csv(paste0(save_path_each,'result_',error.type,'_LASSO_',tau,'_',n,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path_each,'result_',error.type,'_SCAD_',tau,'_',n,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path_each,'result_',error.type,'_MCP_',tau,'_',n,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path_each,'result_',error.type,'_ORACLE_',tau,'_',n,'.csv'),row.names=2)[-1,-1]





