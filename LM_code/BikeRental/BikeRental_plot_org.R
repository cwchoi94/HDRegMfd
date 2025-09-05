
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)


# must include normalizer=(Y.max-Y.min)^(-1) in "inv.clr.density" function.
data.path0 = './LM_code/BikeRental/data/'
result.path0 = './LM_code/BikeRental/result/'
load(paste0(data.path0,'BikeRental.Rdata'))

dim(Y)
Xdata[['spaces']]
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})

cols.all = c('summer','fall','winter','holiday','cloudy or misty','snowy or rainy','precipitation','sunshine',
             'days','wind direction','temperature','wind','humidity','pressure')

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
penalty = penalty.list[2]

# load model
load(paste0(result.path0,'BikeRental_',penalty,'.Rdata'))

# check variable selection (nonzero norm)
proper.beta.norm = f(model)
proper.beta.norm

dim(proper.beta.norm)
Xdata[['p']]

Xdim.max = model$Xdim.max
score.list = predict(model$pca,Xdata)


# For preliminary setting for plotting
pacman::p_load(RColorBrewer)
col.pal = RColorBrewer::brewer.pal(7, "Spectral")
new.color=colorRampPalette(col.pal)

Y.min = 1
Y.max = 24

n.grid = 40
n.t = ncol(model$beta)

# winter
i = 3
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xdata[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=2)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),2)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin,Xmax,length.out=2)
group.season = as.factor(X)
levels(group.season) = list('not winter'=c(0),'winter'=c(1))

plot1 = qplot(data=df,x=t,y=Xbeta,color=group.season,group=group.season,geom='line')+
  labs(color='')+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('Bike rental')+
  ggtitle('Season')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


# precipitation
i = 2
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xdata[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin,Xmax,length.out=4)

plot2 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=new.color(n.grid),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('Bike rental')+
  ggtitle('Precipitation')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


Xdim.max = 2


# wind -> functional covariate
i = 4
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

score = score.list[[idx]][,1:Xdim.max]
Xmin = apply(score,2,min)
Xmax = apply(score,2,max)
beta = model$beta.tensor[[idx]]$element2

## k = 1
k = 1
Xseq = seq(Xmin[k],Xmax[k],length.out=n.grid)
Xbeta = Xseq %*% matrix(beta[k,],nrow=1)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot31 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=new.color(n.grid),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('Bike rental')+
  ggtitle('1st score of wind')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


## k = 2
k = 2
Xseq = seq(Xmin[k],Xmax[k],length.out=n.grid)
Xbeta = Xseq %*% matrix(beta[k,],nrow=1)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot32 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=new.color(n.grid),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('Bike rental')+
  ggtitle('2nd score of wind')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')



eigvectors.wind = model$beta.tensor[[idx]]$element1[1:Xdim.max,]
t = 1:24
category = c(rep(1,24),rep(2,24))
df = data.frame(as.vector(t(eigvectors.wind)),t,category)
names(df) = c('eig','t','category')

group = as.factor(category)
levels(group) = list('1st'=c(1),'2nd'=c(2))

plot33 = ggplot(data=df)+aes(x=t,y=eig,color=group,group=group)+geom_line()+
  labs(color='')+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('wind')+
  ggtitle('Orthonormal bases of wind')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))


#plot_ = grid.arrange(plot1,plot2,plot31,plot32,nrow=2)
plot_ = ggdraw()+
  draw_plot(plot1 ,x=0.0,y=0.63,width=0.5,height=0.37)+
  draw_plot(plot2 ,x=0.5,y=0.63,width=0.5,height=0.37)+
  draw_plot(plot31,x=0.0,y=0.26,width=0.5,height=0.37)+
  draw_plot(plot32,x=0.5,y=0.26,width=0.5,height=0.37)+
  draw_plot(plot33,x=0.25,y=0.0,width=0.5,height=0.26)

plot_
ggsave(paste0(result.path0,'BikeRental.pdf'),plot_,width=10,height=7)





