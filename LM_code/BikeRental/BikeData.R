# Seoul hourly bike rental data in 2020~2021
# data obtained from Korea public data portal
# https://data.seoul.go.kr/dataList/OA-21229/F/1/datasetView.do


data.path0 = './LM_code/BikeRental/data/'
result.path0 = './LM_code/BikeRental/result/'
bike.path0 = paste0(data.path0,'Bike/')

list.files(bike.path0)

dates.bike = c()
BikeRental = c()
for (month in 1:12){
  if (month<10){
    month = paste0('0',month)
  }
  zipname = paste0(bike.path0,'tpss_bcycl_od_statnhm_2020',month,'.zip')
  filename = paste0(bike.path0,'tpss_bcycl_od_statnhm_2020',month,'.csv')
  
  bike = read.csv(unzip(zipname,filename),header=TRUE,fileEncoding='CP949')
  colnames(bike) = c('date','code','time','id.start','name.start','id.end','name.end','count','use.time','use.distance')
  bike = bike[,c('date','time','count')]
  
  dates_each = unique(bike[['date']])
  
  rental = c()
  for (date in dates_each){
    rental.daily = c()
    for (hour in 1:24){
      h1 = 100*(hour-1)
      h2 = 100*hour
      indices = which(bike$date==date & h1<=bike$time & bike$time<h2)
      
      rental.daily = c(rental.daily,sum(bike[indices,'count']))
    }
    rental = rbind(rental,rental.daily)
  }
  rownames(rental) = NULL
  
  dates.bike = c(dates.bike,dates_each)
  BikeRental = rbind(BikeRental,rental)
  
  print(paste0('month: ',month))
}

dates.bike = as.Date(as.character(dates.bike),'%Y%m%d')

# check data
length(dates)
dim(BikeRental)


save(dates.bike,BikeRental,file=paste0(data.path0,'Bike.RData'))





