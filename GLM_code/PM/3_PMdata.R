# Seoul PM data in 2023
# data obtained from air Korea
# https://www.airkorea.or.kr/web/

library(readxl)

data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'
pm.path0 = paste0(data.path0,'PMdata/')

list.files(pm.path0)

year = 2023

dates.all = c()
dates.removed = c()
Y = c()
for (month in 1:12){
  if (month<10){
    month = paste0('0',month)
  }
  
  pm.data = read_excel(paste0(pm.path0,year,month,'.xls'))[,-c(2,3,4)]
  colnames(pm.data) = c('date','SO2','CO','O3','NO2','PM10','PM2.5')
  
  # split date and time
  dates.month.all = unique(pm.data[['date']])
  pm.data[,c('date','time')] = t(sapply(dates.month.all,function(x){strsplit(x,' ')[[1]]}))
  pm.data[,'time'] = sapply(pm.data[,'time'],as.integer)
  pm.data = pm.data[,c('date','time','SO2','CO','O3','NO2','PM10','PM2.5')]
  
  
  # compute binary response
  dates.month.all = unique(pm.data[['date']])
  dates.included = c()
  Y.month = c()
  for (date in dates.month.all){
    indices = which(pm.data==date)
    
    pm.each = pm.data[indices,][['PM2.5']]
    pm.each = as.integer(pm.each)
    # if there are NA values more than 6 hours, we omit this daily data
    if (sum(is.na(pm.each))>6){
      dates.removed = c(dates.removed,date)
      next
    }
    # # check if pm2.5 > 35 more than the consecutive 2 hours
    # rle_result = rle(pm.each>35)
    # Y.each = 1.0 * any(rle_result$values & rle_result$lengths>=2)
    
    # check if (the daily average of pm2.5) > 35
    Y.each = 1.0 * (mean(pm.each,na.rm=TRUE) > 35)
    
    dates.included = c(dates.included,date)
    Y.month = c(Y.month,Y.each)
  }
  
  dates.all = c(dates.all,dates.included)
  Y = c(Y,Y.month)
  
  print(paste0('month: ',month))
}

dates.removed
dates.all

dates.PM = dates.all
Y = matrix(Y)

# check data
length(dates.PM)
summary(as.factor(Y))


save(dates.PM,Y,file=paste0(data.path0,'PM_binary.RData'))





