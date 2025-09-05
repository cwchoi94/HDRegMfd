
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)


# must include normalizer=(Y.max-Y.min)^(-1) in "inv.clr.density" function.
data.path0 = './LM_code/KLIPS/data/'
result.path0 = './LM_code/KLIPS/result/'
load(paste0(data.path0,'KLIPS.Rdata'))
load(paste0(data.path0,'Covariate.Rdata'))
colnames(Covariate)

normalizer = (Y.max - Y.min)^(-1)

Y = Y.density

dim(Y)
Xdata[['spaces']]
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


# defines columns
cols.tot.member = sapply(1:4,function(j){paste0('tot.member',j)})
cols.edu = sapply(1:4,function(j){paste0('edu',j)})
cols.child = sapply(0:2,function(j){paste0('child',j)})
cols.earning = sapply(1:4,function(j){paste0('earning',j)})
cols.expense = sapply(1:3,function(j){paste0('expense',j)})
cols.saving = sapply(1:3,function(j){paste0('saving',j)})
cols.financial = sapply(1:3,function(j){paste0('financial',j)})
cols.age = sapply(1:20,function(j){paste0('age',j)})

cols.real = c('year','covid','Metropolitan','Center','South East','South West','Other region',
              'middle age','old','gender','middle school','college','graduate school','move','estate','debt')
cols.simplex = c(cols.tot.member,cols.child,cols.earning,cols.expense,cols.saving,cols.financial)
cols.density = c(cols.age)

cols.all = c(cols.real,'tot.member','child','earning','expense','saving','financial','age')


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
penalty = penalty.list[3]

# load model
load(paste0(result.path0,'KLIPS_',penalty,'.Rdata'))

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

n.grid = 40
n.t = ncol(model$beta)

# estate -> Self-house
i = 3
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xdata[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin,Xmax,length.out=4)

plot1 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density',limits=c(0,1))+
  ggtitle('Self-house')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


# move -> Moved-in
i = 4
idx = proper.beta.norm['ind',i]
name = colnames(proper.beta.norm)[i]

X = Xdata[[idx]]
Xmin = min(X)
Xmax = max(X)
Xseq = seq(Xmin,Xmax,length.out=n.grid)

mu = model$pca[[idx]]$mu
beta = model$beta.tensor[[idx]]$element2
Xbeta = (Xseq-mu) %*% beta
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin,Xmax,length.out=4)

plot2 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density')+
  ggtitle('Moved-in')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


Xdim.max = 2


# expense -> simplex covariate
i = 8
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
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot31 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density',limits=c(0.3,1.1))+
  ggtitle('1st score of Expense')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


## k = 2
k = 2
Xseq = seq(Xmin[k],Xmax[k],length.out=n.grid)
Xbeta = Xseq %*% matrix(beta[k,],nrow=1)
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot32 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density',limits=c(0,1))+
  ggtitle('2nd score of Expense')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')



eigvectors.expense = model$beta.tensor[[idx]]$element1
eigvectors.expense = inv.clr.simplex(eigvectors.expense)


# age -> density covariate
i = 2
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
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot41 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density',limits=c(0.3,0.9))+
  ggtitle('1st score of Age')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


## k = 2
k = 2
Xseq = seq(Xmin[k],Xmax[k],length.out=n.grid)
Xbeta = Xseq %*% matrix(beta[k,],nrow=1)
Xbeta = inv.clr.BayesHilbert(Xbeta,normalizer=normalizer)

X = as.vector(sapply(Xseq,function(x){rep(x,n.t)}))
t = rep(seq(Y.min,Y.max,length.out=n.t),n.grid)
Xbeta = as.vector(t(Xbeta))
df = data.frame(X,t,Xbeta)

Xseq2 = seq(Xmin[k],Xmax[k],length.out=4)

plot42 = qplot(data=df,x=t,y=Xbeta,color=X,group=as.factor(X),geom='line')+
  scale_colour_gradientn(name='',colours=rev(new.color(n.grid)),breaks=Xseq2,labels=round(Xseq2,2),
                         guide=guide_colorbar(barwidth=10,barheight=0.3,))+
  scale_x_continuous('Income')+
  scale_y_continuous('Density',limits=c(0,1.6))+
  ggtitle('2nd score of Age')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='bottom')


eigvectors.age = model$beta.tensor[[idx]]$element1[1:2,]
eigvectors.age = inv.clr.BayesHilbert(eigvectors.age)
t = seq(1,20,length.out=20)
category = c(rep(1,20),rep(2,20))
df = data.frame(as.vector(t(eigvectors.age)),t,category)
names(df) = c('eig','t','category')

group = as.factor(category)
levels(group) = list('1st'=c(1),'2nd'=c(2))

plot43 = ggplot(data=df)+aes(x=t,y=eig,color=group,group=group)+geom_line()+
  labs(color='')+
  scale_x_continuous('Relative age')+
  scale_y_continuous('Density')+
  ggtitle('Orthonormal bases of Age')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5))



#plot_ = grid.arrange(plot1,plot2,plot31,plot32,nrow=2)
plot_ = ggdraw()+
  draw_plot(plot1 ,x=0.0,y=0.73,width=0.5,height=0.27)+
  draw_plot(plot2 ,x=0.5,y=0.73,width=0.5,height=0.27)+
  draw_plot(plot31,x=0.0,y=0.46,width=0.5,height=0.27)+
  draw_plot(plot32,x=0.5,y=0.46,width=0.5,height=0.27)+
  draw_plot(plot41,x=0.0,y=0.19,width=0.5,height=0.27)+
  draw_plot(plot42,x=0.5,y=0.19,width=0.5,height=0.27)+
  draw_plot(plot43,x=0.25,y=0.0,width=0.5,height=0.19)

plot_
ggsave(paste0(result.path0,'KLIPS.pdf'),plot_,width=10,height=10)

print(round(eigvectors.expense,3))





