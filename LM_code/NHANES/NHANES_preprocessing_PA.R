# require "rnhanesdata" library
# library(devtools)
# devtools::install_github("andrew-leroux/rnhanesdata")

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


data.path0 = './LM_code/NHANES/data/'

# density data
load(paste0(data.path0,'covariate_categorized.Rdata'))



data("PAXINTEN_C");data("PAXINTEN_D")
data("Flags_C");data("Flags_D")

col_vars = paste0("MIN",1:1440)
cov.subjects = Covariate$SEQN
PAXINTEN_C = PAXINTEN_C[PAXINTEN_C$SEQN %in% cov.subjects,]
Flags_C = Flags_C[Flags_C$SEQN %in% cov.subjects,]
PAXINTEN_D = PAXINTEN_D[PAXINTEN_D$SEQN %in% cov.subjects,]
Flags_D = Flags_D[Flags_D$SEQN %in% cov.subjects,]

PAXINTEN_C[,col_vars] = PAXINTEN_C[,col_vars]*Flags_C[,col_vars]
PAXINTEN_D[,col_vars] = PAXINTEN_D[,col_vars]*Flags_D[,col_vars]


cpm.min = 0
cpm.max = 1000
cpm.count.min = 100

q.seq = seq(0,1,length.out=200)

# 2003~2004 density data
subjects = unique(PAXINTEN_C$SEQN)
subjects_C = c()
density_C = c()
quantile_C = c()
for (i in 1:length(subjects)){
  idx = subjects[i]
  y = PAXINTEN_C[which(PAXINTEN_C$SEQN==idx),-(1:5)]
  #if (nrow(y)<4){next}
  y = unlist(y)
  y = y[!is.na(y)]
  y = y[(y>cpm.min) & (y<=cpm.max)]
  if (length(y)<=cpm.count.min){next}
  
  y.density = density(y,from=cpm.min,to=cpm.max)$y
  y.quantile = unname(quantile(y,q.seq))
  
  subjects_C = c(subjects_C,idx)
  density_C = rbind(density_C,y.density)
  quantile_C = rbind(quantile_C,y.quantile)
  print(c(i,length(subjects)))
}

# 2005~2006 density data
subjects = unique(PAXINTEN_D$SEQN)
subjects_D = c()
density_D = c()
quantile_D = c()
for (i in 1:length(subjects)){
  idx = subjects[i]
  y = PAXINTEN_D[which(PAXINTEN_D$SEQN==idx),-(1:5)]
  #if (nrow(y)<4){next}
  y = unlist(y)
  y = y[!is.na(y)]
  y = y[(y>cpm.min) & (y<=cpm.max)]
  if (length(y)<=cpm.count.min){next}
  
  y.density = density(y,from=cpm.min,to=cpm.max)$y
  y.quantile = unname(quantile(y,q.seq))
  
  subjects_D = c(subjects_D,idx)
  density_D = rbind(density_D,y.density)
  quantile_D = rbind(quantile_D,y.quantile)
  print(c(i,length(subjects)))
}

subjects = c(subjects_C,subjects_D)
Y.density = rbind(density_C,density_D)
Y.density = Y.density/(rowMeans(Y.density)*(cpm.max-cpm.min))
Y.quantile = rbind(quantile_C,quantile_D)

# check if Y.density is a density on [cpm.min,cpm.max]
# check if Y.quantile is a quantile 
rowMeans(Y.density)*(cpm.max-cpm.min)
mean(apply(Y.quantile,1,function(y){all(diff(y)>=0)}))

dim(Covariate)
dim(Y.density)
dim(Y.quantile)



# remove outlier Y
Y2 = clr(Y.density)
Y.min = Y2 %>% apply(1,min)
indices = (!is.na(Y.min)) & (Y.min>-10)
Y.density = Y.density[indices,]
Y.quantile = Y.quantile[indices,]
subjects = subjects[indices]

## Age in c(50,85)
Covariate = Covariate[Covariate$SEQN %in% subjects,]
indices = which((Covariate$Age>=50) & (Covariate$Age<=85))
Covariate = Covariate[indices,]
Y = Y.density[indices,]

dim(Y)
dim(Covariate)

save(Y,cpm.min,cpm.max,file=paste0(data.path0,'PA_density.Rdata'))
save(Covariate,file=paste0(data.path0,'covariate_categorized_PA.Rdata'))

Y = Y.quantile[indices,]

dim(Y)
dim(Covariate)

save(Y,cpm.min,cpm.max,file=paste0(data.path0,'PA_quantile.Rdata'))







