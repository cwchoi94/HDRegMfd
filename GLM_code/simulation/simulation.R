### simulation source
### Implemented for each link ('binomial' or 'poisson') and cv.type ('AIC' or 'BIC') 


library(HDRegMfd)
library(hms)


# set basic simulation settings
cv.type = 'AIC'
cv.type = 'BIC'

n.iteration = 100


n.list = c(200,400)
if (link=='binomial'){
  beta.norm.list = c(1.0,2.0)
}else if(link=='poisson'){
  beta.norm.list = c(0.2,0.4)
}


if(link=='binomial'){
  source('simulation_binomial.R') # compute accuracy
}else if (link=='poisson'){
  source('simulation_poisson.R')
}


