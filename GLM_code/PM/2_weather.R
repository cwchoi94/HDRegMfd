# Seoul hourly weather data in 2020~2021
# data obtained from Korea meteorological administration
# https://data.kma.go.kr/resources/html/en/aowdp.html


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'


# weather data
weather.org = read.csv(paste0(data.path0,'weather_hourly.csv'),fileEncoding='CP949')[,c(3,4,6,8,12,14,20)]
colnames(weather.org) = c('date','temperature','precipitation','wind','humidity','pressure','sunshine')
colnames(weather.org)

weather.org[c('date','time')] = matrix(unlist(strsplit(weather.org[,'date'],split=" ")),ncol=2,byrow=TRUE)
weather.org = weather.org[c('date','time','temperature','precipitation','wind','humidity','pressure','sunshine')]

weather.org[is.na(weather.org)] = 0
head(weather.org,20)


# weather and window direction data
weather.daily.org = read.csv(paste0(data.path0,'weather_daily.csv'),fileEncoding='CP949')[,c(3,12,13)]
colnames(weather.daily.org) = c('date','direction','situation')

dates_tmp = weather.daily.org$date

## window direction
direction_tmp = weather.daily.org['direction']
direction_tmp = as.matrix(cbind(cos(direction_tmp),sin(direction_tmp)))
colnames(direction_tmp) = NULL

## weather situation -> written in Korean
situation = weather.daily.org['situation']
situation = apply(situation,1,function(x){gsub("[\\{\\}]", "", regmatches(x, gregexpr("\\{.*?\\}",x))[[1]])})
# table(unlist(situation))

# weather condition
# sit1: misty or cloudy
# sit2: snowy or rainy
sit1 = 1*sapply(situation,function(x){('박무' %in% x) | ('안개' %in% x) | ('연무' %in% x)})
sit2 = 1*sapply(situation,function(x){('눈' %in% x) | ('진눈깨비' %in% x) | ('비' %in% x) | ('소나기' %in% x) | ('안개비' %in% x)})
sit1[which((sit1==1) & (sit2==1))] = 0


# days
days_tmp = 1:length(dates_tmp)
days_tmp = 2*pi*days_tmp/(length(days_tmp)+1)
days_tmp = as.matrix(cbind(cos(days_tmp),sin(days_tmp)))[-length(dates_tmp),]


# split data for each day
ind = c()
temperature = c()
wind = c()
humidity = c()
pressure = c()
precipitation = c()
sunshine = c()
for (i in 1:length(dates_tmp)){
  date = dates_tmp[i]
  indices = which(weather.org$date==date)
  tmp = weather.org[indices,]
  
  if (nrow(tmp)==23){
    tmp = rbind(tmp[1,],tmp)
  } else if (nrow(tmp)<23){
    ind = c(ind,i) # 12/31 data
    next
  }
  
  temperature = rbind(temperature,unlist(tmp['temperature']))
  precipitation = rbind(precipitation,unlist(tmp['precipitation']))
  wind = rbind(wind,unlist(tmp['wind']))
  humidity = rbind(humidity,unlist(tmp['humidity']))
  pressure = rbind(pressure,unlist(tmp['pressure']))
  sunshine = rbind(sunshine,unlist(tmp['sunshine']))
}

weather_tmp = cbind(sit1[-ind],sit2[-ind],rowMeans(precipitation),rowSums(sunshine))
colnames(weather_tmp) = NULL

dates = dates_tmp[-ind]
days = days_tmp[-ind,]
direction = direction_tmp[-ind,]
weather = weather_tmp



# remove NA values
NA.indices = c(which(rowSums(is.na(weather))>0),
               which(rowSums(is.na(days))>0),
               which(rowSums(is.na(direction))>0),
               which(rowSums(is.na(temperature))>0),
               which(rowSums(is.na(wind))>0),
               which(rowSums(is.na(humidity))>0),
               which(rowSums(is.na(pressure))>0))

NA.indices = unique(NA.indices)

dates = dates[-NA.indices]
weather = weather[-NA.indices,,drop=FALSE]
days = days[-NA.indices,,drop=FALSE]
direction = direction[-NA.indices,,drop=FALSE]
temperature = temperature[-NA.indices,,drop=FALSE]
wind = wind[-NA.indices,,drop=FALSE]
humidity = humidity[-NA.indices,,drop=FALSE]
pressure = pressure[-NA.indices,,drop=FALSE]


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










