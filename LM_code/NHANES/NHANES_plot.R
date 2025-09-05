
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)

# must include normalizer=(cpm.max-cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './LM_code/NHANES/data/'
result.path0 = './LM_code/NHANES/result/'
load(paste0(data.path0,'PA_density.Rdata'))
load(paste0(data.path0,'Covariate.Rdata'))
Covariate = Covariate[,-1]

normalizer = (cpm.max-cpm.min)^{-1}

dim(Y)
dim(Covariate)
p = ncol(Covariate)

X.simplex = as.matrix(Covariate[,15:17])
X.simplex = X.simplex/rowSums(X.simplex)
rowSums(X.simplex)

Covariate = Covariate[,-(15:17)]
p = ncol(Covariate)

Yall = Y
Xall = lapply(1:p,function(j){matrix(Covariate[,j],ncol=1)})
Xall[[p+1]] = X.simplex
Xall[['p']] = p+1
Xall[['spaces']] = c(rep('Euclid',p),'simplex')
cols.all = c(colnames(Covariate),'nutrient')

f = function(model){
  proper.indices = model$proper.indices
  Xvar = sapply(proper.indices,function(idx){sum(model$pca[[idx]]$values)})
  proper.beta.norm = rbind(seq(1,length(proper.indices)),proper.indices,model$beta.norm[proper.indices],Xvar)
  rownames(proper.beta.norm) = c('count','ind','beta.norm','varX')
  colnames(proper.beta.norm) = cols.all[proper.indices]
  proper.beta.norm = proper.beta.norm[,order(model$beta.norm[proper.indices],decreasing=TRUE)]
  
  return(proper.beta.norm)
}


######################################
######################################


penalty.list = c('LASSO','SCAD','MCP')
penalty = penalty.list[1]

# load model
load(paste0(result.path0,'NHANES_',penalty,'.Rdata'))

# check variable selection (nonzero norm)
proper.beta.norm = f(model)
proper.beta.norm

dim(proper.beta.norm)
Xall[['p']]

# For preliminary setting for plotting
pacman::p_load(RColorBrewer)
col.pal = RColorBrewer::brewer.pal(7, "Spectral")
new.color=colorRampPalette(col.pal)

n.grid = 2
n.t = ncol(model$beta)

# Mobility
i = 1
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xall[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
X[X==0] = 'N'
X[X==1] = 'Y'
t = rep(seq(cpm.min,cpm.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

plot1 = ggplot(data=df)+aes(x=t,y=Xbeta,color=as.factor(X),group=as.factor(X))+geom_line()+
  labs(color='')+
  scale_x_continuous('cpm')+
  scale_y_continuous('Density')+
  ggtitle('Mobility')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))



# Stroke
i = 2
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xall[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
X[X==0] = 'N'
X[X==1] = 'Y'
t = rep(seq(cpm.min,cpm.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

plot2 = ggplot(data=df)+aes(x=t,y=Xbeta,color=as.factor(X),group=as.factor(X))+geom_line()+
  labs(color='')+
  scale_x_continuous('cpm')+
  scale_y_continuous('Density')+
  ggtitle('Stroke')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))



# Never.married
i = 3
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xall[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
X[X==0] = 'N'
X[X==1] = 'Y'
t = rep(seq(cpm.min,cpm.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

plot3 = ggplot(data=df)+aes(x=t,y=Xbeta,color=as.factor(X),group=as.factor(X))+geom_line()+
  labs(color='')+
  scale_x_continuous('cpm')+
  scale_y_continuous('Density')+
  ggtitle('Never-married')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))


# CHF
i = 4
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xall[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
X[X==0] = 'N'
X[X==1] = 'Y'
t = rep(seq(cpm.min,cpm.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

plot4 = ggplot(data=df)+aes(x=t,y=Xbeta,color=as.factor(X),group=as.factor(X))+geom_line()+
  labs(color='')+
  scale_x_continuous('cpm')+
  scale_y_continuous('Density')+
  ggtitle('CHF')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))




plot_ = ggdraw()+
  draw_plot(plot1,x=0.0,y=0.5,width=0.5,height=0.5)+
  draw_plot(plot2,x=0.5,y=0.5,width=0.5,height=0.5)+
  draw_plot(plot3,x=0.0,y=0.0,width=0.5,height=0.5)+
  draw_plot(plot4,x=0.5,y=0.0,width=0.5,height=0.5)

plot_
ggsave(paste0(result.path0,'NHANES.pdf'),plot_,width=10,height=4.5)




