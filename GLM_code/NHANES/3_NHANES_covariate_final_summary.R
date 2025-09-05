
library(HDRegMfd)

library(knitr)
library(kableExtra)
library(gridExtra)
library(corrplot)
library(reshape2)
library(magrittr)
library(plyr)
library(dplyr)
library(survey)
library(mgcv)
library(refund)
library(rnhanesdata)
RNGkind(sample.kind="Rounding")
library(tidyr)
library(haven)
library(XML)
library(fastDummies)


data.path0 = './GLM_code/NHANES/data/'
load(paste0(data.path0,'PA.Rdata'))
load(paste0(data.path0,'covariate_categorized_PA.Rdata'))

dim(Covariate)

### category to dummy


Covariate = dummy_cols(Covariate,remove_selected_columns=TRUE,remove_most_frequent_dummy=TRUE)
cols = colnames(Covariate)
cols[19:23] = sapply(cols[19:23],function(x){gsub('.*_','',x)}) # marital(19-21), income(22-23)
cols[24:25] = c('BMI.Normal','BMI.Obese') # BMI
cols[26:29] = sapply(cols[26:29],function(x){gsub('.*_','',x)}) # race(26-28), gender(29)
cols[30:34] = sapply(cols[30:34],function(x){gsub('_.*','',x)}) # diabetes(30), CHF(31), CHD(32), Cancer(33), Stroke(34)
cols[35:38] = c('9-11th','Above.college','High.school','Less.9th') # education
cols[39] = 'MobilityProblem'
cols[40:41] = sapply(cols[40:41],function(x){gsub('.*_','',x)}) # DrinkStatus
cols[42:43] = c('Smoking','Smoking.Former') # smoking

colnames(Covariate) = cols


X.simplex = as.matrix(Covariate[,16:18]) # 16-18: carbohydrate, protein, fat
X.simplex = X.simplex/rowSums(X.simplex)
rowSums(X.simplex)

Covariate = Covariate[,-c(1,16:18)] # remove id and simplex covariate

Ycols = c('Diabetes','CHF','CHD','Cancer','Stroke','MobilityProblem')
Ycols = c('Diabetes','CHF','CHD','Cancer','Stroke')

PA = PA.density
save(PA,Covariate,X.simplex,Ycols,log10.cpm.min,log10.cpm.max,file=paste0(data.path0,'PA_density_final.Rdata'))

PA = PA.quantile
save(PA,Covariate,X.simplex,Ycols,log10.cpm.min,log10.cpm.max,file=paste0(data.path0,'PA_quantile_final.Rdata'))



# summary of variables

load(paste0(data.path0,'Covariate.Rdata'))

indices = Covariate$SEQN

load(paste0(data.path0,'Covariate_categorized.Rdata'))

Covariate = Covariate[Covariate$SEQN %in% indices,]


# summary of categorical variables
cols.category = c('Gender','Marital','Income','BMI_cat','Race','EducationAdult','DrinkStatus','SmokeCigs',
                  'Diabetes','CHF','CHD','Cancer','Stroke','MobilityProblem','Hypertensive')

n = nrow(Covariate)
print(n)
for (col in cols.category){
  print(col)
  print(rbind(table(Covariate[,col]),round(100*table(Covariate[,col])/n,2)))
}


# summary of continuous variables
cols.continuous = c('Age','BMI','Systolic.BP','Diastolic.BP','Height','Weight',
                    'Total.chol','HDL','LDL','Triglyceride','Glucose','Kcal','carbohydrate','protein','fat')
round(rbind(apply(Covariate[cols.continuous],2,function(x){c(mean(x),sd(x))})),2)

