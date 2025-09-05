
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)


data.path0 = './GLM_code/PM/data/'
result.path0 = './GLM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))

dim(Y)
Xspaces = Xdata[['spaces']]

Xspaces
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


compute.proper.beta.norm = function(model){
  proper.indices = model$proper.indices
  Xvar = sapply(proper.indices,function(idx){sum(model$pca[[idx]]$values)})
  proper.beta.norm = rbind(seq(1,length(proper.indices)),proper.indices,model$beta.norm[proper.indices],Xvar)
  rownames(proper.beta.norm) = c('count','ind','beta.norm','var.Xi')
  colnames(proper.beta.norm) = Xcols[proper.indices]
  proper.beta.norm = proper.beta.norm[,order(model$beta.norm[proper.indices],decreasing=TRUE)]
  
  return(proper.beta.norm)
}


######################################
######################################

penalty.list = c('LASSO','SCAD','MCP')
penalty = penalty.list[2]

# load model
load(paste0(result.path0,'PM2.5_',penalty,'.Rdata'))

# check variable selection (nonzero norm)
proper.beta.norm = compute.proper.beta.norm(model)

print(penalty)
proper.beta.norm
Xspaces[proper.beta.norm['ind',]]

dim(proper.beta.norm)
Xdata[['p']]

Xdim.max = model$Xdim.max
score.list = predict(model$pca,Xdata)


######################################
######################################


# season
indices = c(2,4,6)

mu.list = list()
beta.list = list()
for (i in indices){
  idx = proper.beta.norm['ind',i]
  name = colnames(proper.beta.norm)[i]
  
  mu.list[[name]] = model$pca[[idx]]$mu
  beta.list[[name]] = model$beta.tensor[[idx]]$element2
}
mu.list = unlist(mu.list)
beta.list = unlist(beta.list)
beta.list['spring'] = - t(mu.list) %*% beta.list
season.effect = beta.list[c('spring','summer','fall','winter')]



# weather
indices = c(1,5)

mu.list = list()
beta.list = list()
for (i in indices){
  idx = proper.beta.norm['ind',i]
  name = colnames(proper.beta.norm)[i]
  
  mu.list[[name]] = model$pca[[idx]]$mu
  beta.list[[name]] = model$beta.tensor[[idx]]$element2
}
mu.list = unlist(mu.list)
beta.list = unlist(beta.list)
beta.list['sunny'] = - t(mu.list) %*% beta.list
weather.effect = beta.list[c('sunny','cloudy or misty','snowy or rainy')]




# For preliminary setting for plotting
pacman::p_load(RColorBrewer)
col.pal = RColorBrewer::brewer.pal(7, "Spectral")
new.color=colorRampPalette(col.pal)

Y.min = 1
Y.max = 24

n.grid = 40
n.t = ncol(model$beta)



# wind
i = 3
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xdata[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=2)

beta = model$beta.tensor[[idx]]
beta.new = beta$element2 %*% beta$element1
beta.new = as.vector(beta.new)

plot1 = qplot(x=1:24,y=beta.new,color=c('wind'),group='wind',geom='line')+
  labs(color='')+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('')+
  ggtitle('Wind-speed')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='')

plot1
ggsave(paste0(result.path0,'wind-speed.pdf'),plot1,width=5,height=3)




season.effect
weather.effect



