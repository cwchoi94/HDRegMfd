

data.path0 = './GLM_code/NHANES/data/'
result.path0 = './GLM_code/NHANES/result/'
load(paste0(data.path0,'PA_final.Rdata'))

Ycols
Ycol = Ycols[5]
result = read.csv(paste0(result.path0,'class_result_avg_density_',Ycol,'.csv'))
result.Euclidean = read.csv(paste0(result.path0,'ncvreg_class_result_avg_',Ycol,'.csv'))

Ycol
result
result.Euclidean

