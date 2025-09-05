


data.path0 = './LM_code/BikeRental/data/'
result.path0 = './LM_code/BikeRental/result/'


# weather data
source(paste0('BikeRental/weather.R'))

# data length
length(dates)

# 4 Euclidean covariate
dim(weather) 

# 2 sphere covariate
dim(days)
dim(direction)

# 4 functional covariates
dim(temperature)
dim(wind)
dim(humidity)
dim(pressure)


# holiday
holiday = read.csv(paste0(data.path0,'holiday.csv'),header=TRUE,row.names=1)
holiday = as.vector(holiday)$x
holiday = sapply(dates,function(date){1*(date %in% holiday)},USE.NAMES=FALSE)

# season
month = matrix(unlist(strsplit(dates,split='-')),ncol=3,byrow=TRUE)[,2]
f = function(x){
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
season = sapply(month,f,USE.NAMES=FALSE)


# Bike data
load(paste0(data.path0,'Bike.Rdata'))
BikeRental = BikeRental[which(dates.bike %in% dates),]
dim(BikeRental)




# Make Y and Xdata

Y = log10(1+BikeRental)

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


save(Y,Xdata,file=paste0(data.path0,'BikeRental.Rdata'))





