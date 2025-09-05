

data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))

result = read.csv(paste0(result.path0,'class_result_avg.csv'))
result.Euclidean = read.csv(paste0(result.path0,'ncvreg_class_result_avg.csv'))

result
result.Euclidean

