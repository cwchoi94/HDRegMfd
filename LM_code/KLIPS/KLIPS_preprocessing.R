
library(HDRegMfd)
library(tidyr)
library(readxl)
library(fastDummies)

data.path0 = './LM_code/KLIPS/data/'
result.path0 = './LM_code/KLIPS/result/'


source('./LM_code/KLIPS/KLIPS_preprocessing_ftn.R')



Y.density.org = c()
Y.quantile.org = c()
Covariate.org = c()

for (year in 2017:2021){
  result = get_KLIPS_data_each(year)
  Y.density.org = rbind(Y.density.org,result$Y.density)
  Y.quantile.org = rbind(Y.quantile.org,result$Y.quantile)
  Covariate.org = rbind(Covariate.org,result$X)
}

# change to simplex value
cols.tot.member = sapply(1:4,function(j){paste0('tot.member',j)})
cols.edu = sapply(1:4,function(j){paste0('edu',j)})
cols.child = sapply(0:2,function(j){paste0('child',j)})
cols.earning = sapply(1:4,function(j){paste0('earning',j)})
cols.expense = sapply(1:3,function(j){paste0('expense',j)})
cols.saving = sapply(1:3,function(j){paste0('saving',j)})
cols.financial = sapply(1:3,function(j){paste0('financial',j)})
cols.age = sapply(1:20,function(j){paste0('age',j)})
cols.age2 = sapply(1:20,function(j){paste0('age.quantile',j)})

Covariate.org[,cols.tot.member] = Covariate.org[,cols.tot.member]/rowSums(Covariate.org[,cols.tot.member])
Covariate.org[,cols.edu] = Covariate.org[,cols.edu]/rowSums(Covariate.org[,cols.edu])
Covariate.org[,cols.child] = Covariate.org[,cols.child]/rowSums(Covariate.org[,cols.child])
Covariate.org[,cols.earning] = Covariate.org[,cols.earning]/rowSums(Covariate.org[,cols.earning])
Covariate.org[,cols.expense] = Covariate.org[,cols.expense]/rowSums(Covariate.org[,cols.expense])
Covariate.org[,cols.saving] = Covariate.org[,cols.saving]/rowSums(Covariate.org[,cols.saving])
Covariate.org[,cols.financial] = Covariate.org[,cols.financial]/rowSums(Covariate.org[,cols.financial])

row.na1 = rowSums(is.na(clr(Covariate.org[,cols.tot.member])))
row.na2 = rowSums(is.na(clr(Covariate.org[,cols.child])))
row.na3 = rowSums(is.na(clr(Covariate.org[,cols.earning])))
row.na4 = rowSums(is.na(clr(Covariate.org[,cols.expense])))
row.na5 = rowSums(is.na(clr(Covariate.org[,cols.saving])))
row.na6 = rowSums(is.na(clr(Covariate.org[,cols.financial])))
row.naY = rowSums(is.na(clr(Y.density.org)))
# row.naY2 = (apply(abs(clr(Y.density.org)),1,max)>10)+0
row.na = row.na1+row.na2+row.na3+row.na4+row.na5+row.na6+row.naY#+row.naY2
row.na[is.na(row.na)] = 1
idx = which(row.na>0)

print(idx) # idx 70: NA in Y.density. idx 88: zero in expense
Y.density = Y.density.org[-idx,]
Y.quantile = Y.quantile.org[-idx,]
Covariate.org = Covariate.org[-idx,]

apply(abs(clr(Y.density)),1,max)


# make dummy variables
head(Covariate.org)
Covariate = dummy_cols(Covariate.org,select_columns=c('region','age.group'),remove_selected_columns=TRUE)
p1 = ncol(Covariate)
colnames(Covariate)
colnames(Covariate)[(p1-8):p1] = c('Seoul','Metropolitan','Center','South East','South West','Other region','young','middle age','old')

## remove "Seoul", "young" and "edu2" variabes
cols.real = c('year','covid','Metropolitan','Center','South East','South West','Other region',
              'middle age','old','gender','edu1','edu3','edu4','move','estate','debt')
cols.simplex = c(cols.tot.member,cols.child,cols.earning,cols.expense,cols.saving,cols.financial)
cols.density = c(cols.age,cols.age2)
cols.all = c(cols.real,cols.simplex,cols.density)
Covariate = Covariate[,cols.all]
colnames(Covariate)


# make Xdata

p.real = length(cols.real)

Xdata = list()
# real covariates
for (j in 1:p.real){
  col = cols.real[j]
  Xdata[[j]] = matrix(Covariate[,col],ncol=1)
}

# simplex covariates
## tot.member
cols = cols.tot.member
Xdata[[p.real+1]] = as.matrix(Covariate[,cols])

## child
cols = cols.child
Xdata[[p.real+2]] = as.matrix(Covariate[,cols])

## earning
cols = cols.earning
Xdata[[p.real+3]] = as.matrix(Covariate[,cols])

## expense
cols = cols.expense
Xdata[[p.real+4]] = as.matrix(Covariate[,cols])

## saving
cols = cols.saving
Xdata[[p.real+5]] = as.matrix(Covariate[,cols])

## financial
cols = cols.financial
Xdata[[p.real+6]] = as.matrix(Covariate[,cols])


# density
## age
cols = cols.age
Xdata[[p.real+7]] = as.matrix(Covariate[,cols])

Xdata[['p']] = p.real+7
Xdata[['spaces']] = c(rep('Euclid',p.real),rep('simplex',6),'BayesHilbert')


save(Y.density,Y.quantile,Xdata,Y.min,Y.max,file=paste0(data.path0,'KLIPS.Rdata'))
save(Covariate,Covariate.org,file=paste0(data.path0,'Covariate.Rdata'))



# density
## age2
cols = cols.age2
Xdata[[p.real+7]] = as.matrix(Covariate[,cols])

Xdata[['p']] = p.real+7
Xdata[['spaces']] = c(rep('Euclid',p.real),rep('simplex',6),'Wasserstein')


save(Y.density,Y.quantile,Xdata,Y.min,Y.max,file=paste0(data.path0,'KLIPS2.Rdata'))


