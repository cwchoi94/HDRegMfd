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


# covariate data
data.path0 = './GLM_code/NHANES/data/'
load(paste0(data.path0,'covariate_categorized.Rdata'))


data("PAXINTEN_D")
data("Flags_D")

col_vars = paste0("MIN",1:1440)
cov.subjects = Covariate$SEQN
PAXINTEN_D = PAXINTEN_D[PAXINTEN_D$SEQN %in% cov.subjects,]
Flags_D = Flags_D[Flags_D$SEQN %in% cov.subjects,]

PAXINTEN_D[,col_vars] = PAXINTEN_D[,col_vars]*Flags_D[,col_vars]


log10.cpm.min = 0
log10.cpm.max = 3
cpm.count.min = 100

q.seq = seq(0,1,length.out=200)

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
  y = log10(y)
  y = y[!is.na(y)]
  y = y[(y>=log10.cpm.min) & (y<=log10.cpm.max)]
  if (length(y)<=cpm.count.min){next}
  
  y.density = density(y,from=log10.cpm.min,to=log10.cpm.max)$y
  y.quantile = unname(quantile(y,q.seq))
  
  subjects_D = c(subjects_D,idx)
  density_D = rbind(density_D,y.density)
  quantile_D = rbind(quantile_D,y.quantile)
  print(c(i,length(subjects)))
}

subjects = subjects_D
Y.density = density_D
Y.quantile = quantile_D
Y.density = Y.density/(rowMeans(Y.density)*(log10.cpm.max-log10.cpm.min))

# check if Y.density is a density on [log10.cpm.min,log10.cpm.max]
rowMeans(Y.density)*(log10.cpm.max-log10.cpm.min)
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


PA.density = Y.density[indices,]
PA.quantile = Y.quantile[indices,]
rownames(PA.density) = NULL
rownames(PA.quantile) = NULL

dim(PA.density)
dim(PA.quantile)
dim(Covariate)

save(PA.density,PA.quantile,log10.cpm.min,log10.cpm.max,file=paste0(data.path0,'PA.Rdata'))
save(Covariate,file=paste0(data.path0,'covariate_categorized_PA.Rdata'))



