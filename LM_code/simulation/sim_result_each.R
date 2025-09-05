### simulation results




sim = 'sim5'

save_path = paste0('./LM_code/simulation/sim_result/',sim,'/')
save_rdata_path = paste0(save_path,'Rdata/')

load(paste0(save_path,'sim_info.RData'))
load(paste0(save_path,'parameters.RData'))

sim.info$Yspace
parameters$Xspaces
parameters$Xdims

# basic information

sim.info
parameters$proper.indices


# total result

read.csv(paste0(save_path,'result_lasso.csv'))[-1,-1]
read.csv(paste0(save_path,'result_scad.csv'))[-1,-1]
read.csv(paste0(save_path,'result_mcp.csv'))[-1,-1]
read.csv(paste0(save_path,'result_oracle.csv'))[-1,-1]


# each result

count = 1
read.csv(paste0(save_path,'result_lasso_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_scad_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_mcp_',count,'.csv'),row.names=2)[-1,-1]
read.csv(paste0(save_path,'result_oracle_',count,'.csv'),row.names=2)[-1,-1]



# Compare two simulation result

read.csv(paste0('./LM_code/simulation/sim_result/sim1/result_','mcp','.csv'))[-1,-1]
read.csv(paste0('./LM_code/sim_result/sim1/result_','oracle','.csv'))[-1,-1]

count = 3
read.csv(paste0('./LM_code/simulation/sim_result/sim1/result_','lasso','_',count,'.csv'))[-1,-1]#[1:5,]
read.csv(paste0('./LM_code/simulation/sim_result/sim12/result_','lasso','_',count,'.csv'))[-1,-1]#[1:5,]


# Compare two models

count = 1
iter = 1
modelname = 'scad'
load(paste0('./LM_code/simulation/sim_result/sim2/Rdata/',modelname,'_',count,'_',iter,'.RData'))
model1 = model
load(paste0('./LM_code/simulation/sim_result/sim2/Rdata/oracle_',count,'_',iter,'.RData'))
model2 = model

model1$proper.indices
model2$proper.indices
c(model2$lambda,model2$K.max)
c(model3$lambda,model3$K.max)






