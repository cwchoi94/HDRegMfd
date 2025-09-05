

library(tidyr)
library(readxl)

data.path0 = './LM_code/KLIPS/data/'
result.path0 = './LM_code/KLIPS/result/'


Y.min = 2.5
Y.max = 4



get_KLIPS_data_each = function(year){
  wave = year - 1997
  data = read_excel(paste0(data.path0,'eklips',wave,'h.xlsx'),col_types='numeric')
  
  # remove na household income
  col.number = '2102'
  col = paste0('h',wave,col.number)
  indices = !is.na(data[,col])
  data = data[indices,]
  
  n = nrow(data)
  p = ncol(data)
  
  #' brief description of coavriates
  #' 
  #' 0141 ~ 0142: address
  #' 0150: total number of household members
  #' 0151: Household attrition -> ???
  #' 
  #' Household member information
  #' 0241 ~ 0255: gender
  #' 0361 ~ 0375: age
  #' 0661 ~ 0675: education
  #' 1401: moved to another place or not
  #' 1402: reason for moving to current place
  #' 1406: ownership of home
  #' 1407: type of residence
  #' 1501: existence of children in high school or younger
  #' 2001: existence of children in college or higher
  #' 
  #' Income
  #' 2101: Whether any family members with earned income
  #' 2102: average total earnings (unit: KRW 10,000)
  #' 
  #' Expenditures, savings, assets, debt
  #' 2301: average monthly living expenses
  #' 2401: existence of savings
  #' 2402: average monthly savings
  #' 2501: own real estate other than current dwelling
  #' 2632: existence of debts
  #' 2705: current financial condition
  
  
  # find owners among household members and their gender, age and education
  # gender: 1(male), 2(female)
  # age.group: 1(young), 2(middle), 3(old)
  # education: 1(middle school), 2(high school), 3(college), 4(graduate)
  
  ## region
  col.number = '0141'
  col = paste0('h',wave,col.number)
  table(data[,col])
  
  region.col1 = c(1) # Seoul
  region.col2 = c(5,8) # Incheon, Gyeonggi
  region.col3 = c(4,10,11,19) # Daejeon, Sejong, Chungcheong
  region.col4 = c(2,3,7,14,15) # Busan, Daegu, Ulsan, Gyeongsang
  region.col5 = c(6,12,13) # Gwanju, Jeonra
  region.col6 = c(9,16) # Jeju, Gangwon
  
  region.cols = c('Seoul','Metropolitan','Center','South East','South West','Other')
  
  region = as.data.frame(data[,col])
  region[[col]] = as.factor(region[[col]])
  levels(region[[col]]) = list('1'=region.col1,'2'=region.col2,'3'=region.col3,
                               '4'=region.col4,'5'=region.col5,'6'=region.col6)
  table(region[[col]])
  data['region'] = region
  
  ## household owner's gender, age and education
  data.owner = data.frame(matrix(NA,n,3))
  colnames(data.owner) = c('gender','age','education')
  for (col.number in 1:15){
    col = paste0('h',wave,'02',60+col.number)
    indices = which(data[,col]==10)
    
    col.gender = paste0('h',wave,'02',40+col.number)
    col.age = paste0('h',wave,'03',60+col.number)
    col.education = paste0('h',wave,'06',60+col.number)
    
    data.owner[indices,] = as.data.frame(data[indices,c(col.gender,col.age,col.education)])
  }
  data[,c('gender','age','education')] = data.owner
  data[['gender']] = cut(data[['gender']],breaks=c(-Inf,1.5,Inf),labels=c('1','2'))
  data[['age.group']] = cut(data[['age']],breaks=c(-Inf,39.5,59.5,Inf),labels=c('1','2','3'))
  data[['education']] = cut(data[['education']],breaks=c(-Inf,4.5,5.5,7.5,Inf),labels=c('1','2','3','4'))
  # data[['education']] = apply(data[sapply(1:15,function(j){paste0('h',wave,'06',60+j)})],1,function(x){max(x[!is.na(x)])})
  
  
  # compute household density for each group of (region,age.group)
  # n=15 (=5*3) observation for each year 
  levels.region = levels(data[['region']])
  levels.gender = levels(data[['gender']])
  levels.age.group = levels(data[['age.group']])
  
  Y.density = c()
  Y.quantile = c()
  Covariate = c()
  Yorg = list()
  iter = 0
  for (region in levels.region){
    for (age.group in levels.age.group){
      iter = iter + 1
      indices = which((data[['region']]==region) & (data[['age.group']]==age.group))
      
      # household income density
      col.number = '2102'
      col = paste0('h',wave,col.number)
      
      q.seq = seq(0,1,length.out=200)
      
      y = data[indices,col][[col]]
      y[y<10^(Y.min)] = 10^(Y.min)
      y[y>10^(Y.max)] = 10^(Y.max)
      y = log10(y)
      y.density = density(y,from=Y.min,to=Y.max)$y
      y.quantile = unname(quantile(y,q.seq))
      
      Y.density = rbind(Y.density,y.density)
      Y.quantile = rbind(Y.quantile,y.quantile)
      Yorg[[iter]] = y
      
      # year, region, gender, age.group, covid
      gender = data[indices,'gender']
      gender = mean(gender==1)
      covid = 1*(year %in% c(2020,2021))
      info = c(year,region,age.group,covid)
      info = as.integer(info)
      names(info) = c('year','region','age.group','covid')
      
      # gender
      gender = data[indices,'gender']
      gender = mean(gender==1)
      names(gender) = c('gender')
      
      # total number of household members
      col.number = '0150'
      col = paste0('h',wave,col.number)
      num.member = data[indices,col]
      num.member[num.member>4] = 4
      num.member = table(factor(num.member[[col]],levels=1:4))
      names(num.member) = sapply(1:4,function(j){paste0('tot.member',j)})
      
      # education
      education = data[indices,'education']
      education = table(education)
      names(education) = sapply(1:4,function(j){paste0('edu',j)})
      
      # ratio of moved to another place
      col.number = '1401'
      col = paste0('h',wave,col.number)
      move = data[indices,col]
      move = mean(move==1)
      names(move) = 'move'
      
      # number of child
      col.number = '1502'
      col = paste0('h',wave,col.number)
      child = data[indices,col][[col]]
      child[is.na(child)] = 0
      child[child>2] = 2
      child = table(factor(child,levels=0:2))
      names(child) = sapply(0:2,function(j){paste0('child',j)})
      
      # total earning
      col.number = '2102'
      col = paste0('h',wave,col.number)
      earning = data[indices,col][[col]]
      earning = cut(earning,breaks=c(-Inf,2000,4000,6000,Inf),labels=c(1,2,3,4))
      earning = table(earning)
      names(earning) = sapply(1:4,function(j){paste0('earning',j)})
      
      # monthly living expense
      col.number = '2301'
      col = paste0('h',wave,col.number)
      expense = data[indices,col][[col]]
      expense = cut(expense,breaks=c(-Inf,150,300,Inf),labels=c(1,2,3))
      expense = table(expense)
      names(expense) = sapply(1:3,function(j){paste0('expense',j)})
      
      # monthly savings
      col.number = '2402'
      col = paste0('h',wave,col.number)
      saving = data[indices,col][[col]]
      saving[is.na(saving)] = 0
      saving = cut(saving,breaks=c(-Inf,50,100,Inf),labels=c(1,2,3))
      saving = table(saving)
      names(saving) = sapply(1:3,function(j){paste0('saving',j)})
      
      # own real state
      col.number = '2501'
      col = paste0('h',wave,col.number)
      estate = data[indices,col][[col]]
      estate = mean(estate==1)
      names(estate) = 'estate'
      
      # debt
      col.number = '2632'
      col = paste0('h',wave,col.number)
      debt = data[indices,col][[col]]
      debt = mean(debt==1)
      names(debt) = 'debt'
      
      # current financial condition
      col.number = '2705'
      col = paste0('h',wave,col.number)
      financial = data[indices,col][[col]]
      financial[financial<=2] = 1
      financial[financial==3] = 2
      financial[financial>=4] = 3
      financial = table(factor(financial,levels=1:3))
      names(financial) = sapply(1:3,function(j){paste0('financial',j)})
      
      # age density
      age = data[indices,'age'][['age']]
      if (age.group==1){
        age[age<20] = 20
        age = density(age,n=20,from=20,to=40)$y
      } else if (age.group==2){
        age = density(age,n=20,from=40,to=60)$y
      } else if (age.group==3){
        age[age>80] = 79
        age = density(age,n=20,from=60,to=80)$y
      }
      names(age) = sapply(1:20,function(j){paste0('age',j)})
      
      # age density2
      q.seq = seq(0,1,length.out=20)
      age2 = data[indices,'age'][['age']]
      age2 = unname(quantile(age2,q.seq))
      if (age.group==1){
        age2[age2<20] = 20
        age2 = unname(quantile(age2,q.seq)) - 20
      } else if (age.group==2){
        age2 = unname(quantile(age2,q.seq)) - 40
      } else if (age.group==3){
        age2[age2>80] = 80
        age2 = unname(quantile(age2,q.seq)) - 60
      }
      names(age2) = sapply(1:20,function(j){paste0('age.quantile',j)})
      
      
      # append
      each.covariate = c(info,gender,num.member,education,move,child,earning,expense,saving,estate,debt,financial,age,age2)
      Covariate = rbind(Covariate,each.covariate)
      
    }
  }
  rownames(Y.density) = NULL
  rownames(Y.quantile) = NULL
  rownames(Yorg) = NULL
  rownames(Covariate) = NULL
  
  # normalize to density
  cols.age = sapply(1:20,function(j){paste0('age',j)})
  cols.age2 = sapply(1:20,function(j){paste0('age.quantile',j)})
  Y.density = Y.density/(rowMeans(Y.density)*(Y.max-Y.min))
  Covariate[,cols.age] = Covariate[,cols.age]/(rowMeans(Covariate[,cols.age])*20)
  
  result = list(Y.density=Y.density,Y.quantile=Y.quantile,X=Covariate,Yorg=Yorg)
  return(result)
}









