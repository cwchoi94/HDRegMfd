### simulation results



cv.type = 'AIC'

p = 100
s = 2

link = 'binomial'
link = 'poisson'

save_path = paste0('./GLM_code/simulation/sim_result/',link,'_',cv.type,'_',p,'_',s,'/')
save_rdata_path = paste0(save_path,'Rdata/')

load(paste0(save_path,'sim_info.RData'))
load(paste0(save_path,'parameters.RData'))

parameters$Xspaces
parameters$Xdims

# basic information

sim.info
parameters$proper.ind.mat
parameters$proper.indices


# total result
base.model.name = 'NW'
# base.model.name = 'LL'
# base.model.name = 'LM'
read.csv(paste0(save_path,'result_lasso.csv'))[-1,-1]
read.csv(paste0(save_path,'result_scad.csv'))[-1,-1]
read.csv(paste0(save_path,'result_mcp.csv'))[-1,-1]
read.csv(paste0(save_path,'result_oracle.csv'))[-1,-1]


# each result
count = 1
save_path_each = paste0(save_path,'result_each_',base.model.name,'/')
read.csv(paste0(save_path,'result_lasso_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_scad_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_mcp_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_oracle_',count,'.csv'),row.names=2)[-1,-1]





