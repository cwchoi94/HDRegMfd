# require "rnhanesdata" library
# library(devtools)
# devtools::install_github("andrew-leroux/rnhanesdata")

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


# demographic data
# use data made by "rnhanesdata" library: https://github.com/andrew-leroux/rnhanesdata
# https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Dietary&CycleBeginYear=2005

data.path0 = './GLM_code/NHANES/data/'
data("Covariate_D")



# 2005~2006 NHANES data

data.path = 'https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/'

Covariate = Covariate_D
Covariate = Covariate[Covariate['RIDAGEYR']>=20,][,-c(2:8)]
colnames(Covariate) = c('SEQN','Age',colnames(Covariate)[-(1:2)])

## add more demographic data: marital, income
## description: https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/DEMO_D.htm
# demographic = read_xpt(paste0(data.path,'DEMO_D.XPT'))
demographic = read_xpt(paste0(data.path,'DEMO_D.XPT'))
demographic = data.frame(demographic[c('SEQN','DMDMARTL','INDHHINC')])
colnames(demographic) = c('SEQN','Marital','Income')


Covariate = merge(demographic,Covariate,by='SEQN')
rm(demographic)


## blood presure: systolic, diastolic
## https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/BPX_D.htm
BP = read_xpt(paste0(data.path,'BPX_D.XPT'))
BP$Systolic.BP = (BP$BPXSY1+BP$BPXSY2+BP$BPXSY3)/3
BP$Diastolic.BP = (BP$BPXDI1+BP$BPXDI2+BP$BPXDI3)/3
BP = data.frame(BP[c('SEQN','Systolic.BP','Diastolic.BP')])
BP$Hypertensive = 1*((BP$Systolic.BP>=140) | (BP$Diastolic.BP>=90))

Covariate = merge(Covariate,BP,by='SEQN')
rm(BP)

## body index: weight, tall
## https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/BMX_D.htm
BM = read_xpt(paste0(data.path,'BMX_D.XPT'))
BM = data.frame(BM[c('SEQN','BMXHT','BMXWT','BMXBMI')])
colnames(BM) = c('SEQN','Height','Weight')
BM = BM[c('SEQN','Height','Weight')]

Covariate = merge(Covariate,BM,by='SEQN')
rm(BM)


## cholesterol: total, HDL, LDL, triglyceride (mg/dL)
## https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/HDL_D.htm, etc
## total
total = read_xpt(paste0(data.path,'TCHOL_D.XPT'))
total = data.frame(total[c('SEQN','LBXTC')])
colnames(total) = c('SEQN','Total.chol')
### HDL
HDL = read_xpt(paste0(data.path,'HDL_D.XPT'))
HDL = data.frame(HDL[c('SEQN','LBDHDD')])
colnames(HDL) = c('SEQN','HDL')
### LDL
LDL = read_xpt(paste0(data.path,'TRIGLY_D.XPT'))
LDL = data.frame(LDL[c('SEQN','LBDLDL','LBXTR')])
colnames(LDL) = c('SEQN','LDL','Triglyceride')
### merge
cholesterol = merge(total,HDL,by='SEQN')
cholesterol = merge(cholesterol,LDL,by='SEQN')

Covariate = merge(Covariate,cholesterol,by='SEQN')
rm(total,HDL,LDL,cholesterol)


## plasma glucose
## https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/GLU_D.htm
glucose = read_xpt(paste0(data.path,'GLU_D.XPT'))
glucose = data.frame(glucose[c('SEQN','LBXGLU')])
colnames(glucose) = c('SEQN','Glucose')

Covariate = merge(Covariate,glucose,by='SEQN')
rm(glucose)


## intake: total energy, carbohydrate, protein, fat
## https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/DR1TOT_D.htm, etc
intake = read_xpt(paste0(data.path,'DR1TOT_D.XPT'))
intake2 = read_xpt(paste0(data.path,'DR2TOT_D.XPT'))
intake$Kcal = (intake$DR1TKCAL + intake2$DR2TKCAL)/2
intake$carbohydrate = (intake$DR1TCARB + intake2$DR2TCARB)/2
intake$protein = (intake$DR1TPROT + intake2$DR2TPROT)/2
intake$fat = (intake$DR1TTFAT + intake2$DR2TTFAT)/2
intake = data.frame(intake[c('SEQN','Kcal','carbohydrate','protein','fat')])

Covariate = merge(Covariate,intake,by='SEQN')
rm(intake,intake2)


save(Covariate,file=paste0(data.path0,'covariate_org.Rdata'))



## remove poor observations

load(paste0(data.path0,'covariate_org.Rdata'))

Covariate$Marital = as.factor(Covariate$Marital)
Covariate$Income = as.factor(Covariate$Income)

Covariate = Covariate[apply(Covariate,1,function(x){sum(is.na(x))})==0,]
Covariate = Covariate[Covariate$Income %in% 1:12,]
Covariate = Covariate[Covariate$Diabetes %in% c('No','Yes'),]
Covariate = Covariate[Covariate$CHF %in% c('No','Yes'),]
Covariate = Covariate[Covariate$CHD %in% c('No','Yes'),]
Covariate = Covariate[Covariate$Cancer %in% c('No','Yes'),]
Covariate = Covariate[Covariate$Stroke %in% c('No','Yes'),]
Covariate = Covariate[!(Covariate$EducationAdult %in% c('Refused',"Don't know")),]


### re-categorize categorical covariates

levels(Covariate$Marital) = list('Married'=c(1,6),'Widowed'=c(2),'Separated'=c(3,4),'Never.married'=c(5))
levels(Covariate$Income) = list('low.income'=1:4,'middle.income'=5:8,'high.income'=9:12)
levels(Covariate$BMI_cat) = list('Normal'=c('Normal','Underweight'),'Overweight'=c('Overweight'),'Obese'=c('Obese'))
levels(Covariate$Race) = list('White'=c('White'),'Mexican'=c('Mexican American'),'Black'=c('Black'),'Other'=c('Other Hispanic','Other'))
levels(Covariate$Diabetes) = list('No'=c('No'),'Yes'=c('Yes'))
levels(Covariate$CHF) = list('No'=c('No'),'Yes'=c('Yes'))
levels(Covariate$CHD) = list('No'=c('No'),'Yes'=c('Yes'))
levels(Covariate$Cancer) = list('No'=c('No'),'Yes'=c('Yes'))
levels(Covariate$Stroke) = list('No'=c('No'),'Yes'=c('Yes'))
levels(Covariate$EducationAdult) = list('Less than 9th grade'=c('Less than 9th grade'),
                                        '9-11th grade'=c('9-11th grade'),
                                        'High school grad/GED or equivalent'=c('HIgh school grad/GED or equivalent'),
                                        'Some College or AA degree'=c('Some College or AA degree'),
                                        'College graduate or above'=c('College graduate or above'))


save(Covariate,file=paste0(data.path0,'covariate_categorized.Rdata'))









