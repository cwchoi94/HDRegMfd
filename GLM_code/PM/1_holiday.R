

# Korea holiday data
# crawling from Korea public data portal
# https://www.data.go.kr/index.do


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'


## public holiday

library(jsonlite)
library(httr)

service.key = NULL
service.key = 'ao5aQcZGJQ3K7yoDk%2BSGRIZ5FFF4EsNh%2BJGnL%2Bwz8X2CHrBzvLt8wKgD3ENDZqGClWaaXYRUT1dk0limUOzjqQ%3D%3D'

year.list = c(2023)

holiday = c()
for (year in year.list){
  for (month in 1:12){
    if (month<10){
      month = paste0('0',month)
    }
    
    url.path = 'http://apis.data.go.kr/B090041/openapi/service/SpcdeInfoService/getRestDeInfo'
    url.path = paste0(url.path,'?serviceKey=',service.key,'&solYear=',year,'&solMonth=',month)
    
    data = GET(url.path)
    
    data = content(data,as='text',encoding='UTF-8')
    data = fromJSON(data)
    
    if (!(data$response$body$items=="")){
      data = data$response$body$items$item
      holiday = rbind(holiday,data)
    }
    
    print(paste0('year: ',year,', month: ',month))
  } 
}


holiday = as.Date(as.character(holiday$locdate),format='%Y%m%d')

## Saturday and Sunday

df = data.frame(seq(as.Date(paste0(year.list[1],'-01-01')),as.Date(paste0(year.list[length(year.list)],'-12-31')),1))
names(df) = 'date'
df$weekday = weekdays(as.Date(df$date))
df = df[(df$weekday=='토요일') | (df$weekday=='일요일'),]

## append

holiday = unique(sort(c(holiday,df$date)))
write.csv(holiday,paste0(data.path0,'holiday.csv'))


# test

holiday = read.csv(paste0(data.path0,'holiday.csv'),header=TRUE,row.names=1)
holiday = as.vector(holiday)$x
holiday




