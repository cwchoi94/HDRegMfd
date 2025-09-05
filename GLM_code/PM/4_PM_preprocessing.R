


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'


# weather data
source(paste0('./GLM_code/PM/2_weather.R'))

# PM2.5 data
load(paste0(data.path0,'PM_binary.Rdata'))

# compute the intersection of dates
dates.all = intersect(dates,dates.PM)
Y = Y[which(dates.PM %in% dates.all),,drop=FALSE]
dim(Y)

indices = which(dates %in% dates.all)
weather = weather[indices,,drop=FALSE]
days = days[indices,,drop=FALSE]
direction = direction[indices,,drop=FALSE]
temperature = temperature[indices,,drop=FALSE]
wind = wind[indices,,drop=FALSE]
humidity = humidity[indices,,drop=FALSE]
pressure = pressure[indices,,drop=FALSE]


# holiday
holiday = read.csv(paste0(data.path0,'holiday.csv'),header=TRUE,row.names=1)
holiday = as.vector(holiday)$x
holiday = sapply(dates.all,function(date){1*(date %in% holiday)},USE.NAMES=FALSE)

# season
month = matrix(unlist(strsplit(dates.all,split='-')),ncol=3,byrow=TRUE)[,2]
compute.season = function(x){
  x = as.integer(x)
  if (x %in% c(3,4,5)){
    y = 'spring'
  } else if (x %in% c(6,7,8)){
    y = 'summer'
  } else if (x %in% c(9,10,11)){
    y = 'fall'
  } else if (x %in% c(12,1,2)){
    y = 'winter'
  }
  return(y)
}
season = sapply(month,compute.season,USE.NAMES=FALSE)


# data dimension

## PM2.5 data
dim(Y)

## data length
length(dates.all)

## 4 Euclidean covariate
dim(weather) 
length(holiday)
length(season)

## 2 sphere covariate
dim(days)
dim(direction)

## 4 functional covariates
dim(temperature)
dim(wind)
dim(humidity)
dim(pressure)



# Make Xdata
Xdata = list()
Xdata[[1]] = matrix(1*(season=='summer'),ncol=1)
Xdata[[2]] = matrix(1*(season=='fall'),ncol=1)
Xdata[[3]] = matrix(1*(season=='winter'),ncol=1)
Xdata[[4]] = matrix(holiday,ncol=1)
Xdata[[5]] = matrix(weather[,1],ncol=1)
Xdata[[6]] = matrix(weather[,2],ncol=1)
Xdata[[7]] = matrix(weather[,3],ncol=1)
Xdata[[8]] = matrix(weather[,4],ncol=1)
Xdata[[9]] = days
Xdata[[10]] = direction
Xdata[[11]] = temperature
Xdata[[12]] = wind
Xdata[[13]] = humidity
Xdata[[14]] = pressure

Xdata[['p']] = 14
Xdata[['spaces']] = c(rep('Euclid',8),rep('sphere',2),rep('functional',4))

Xcols = c('summer','fall','winter','holiday','cloudy or misty','snowy or rainy','precipitation','sunshine',
          'days','wind direction','temperature','wind','humidity','pressure')


save(Y,Xdata,Xcols,file=paste0(data.path0,'PM2.5.Rdata'))





