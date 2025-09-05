
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)

# must include normalizer=(cpm.max-cpm.min)^(-1) in "inv.clr.density" function.

data.path0 = './GLM_code/NHANES/data/'
result.path0 = './GLM_code/NHANES/result/'
result.path = paste0(result.path0,'density/')
load(paste0(data.path0,'PA_final.Rdata'))

normalizer = (log10.cpm.max-log10.cpm.min)^{-1}


cols.all = c(colnames(Covariate),'nutrient','PA')

compute.proper.beta.norm = function(model,Xcols){
  proper.indices = model$proper.indices
  Xvar = sapply(proper.indices,function(idx){sum(model$pca[[idx]]$values)})
  proper.beta.norm = rbind(seq(1,length(proper.indices)),proper.indices,model$beta.norm[proper.indices],Xvar)
  rownames(proper.beta.norm) = c('count','ind','beta.norm','varX')
  colnames(proper.beta.norm) = Xcols[proper.indices]
  proper.beta.norm = proper.beta.norm[,order(model$beta.norm[proper.indices],decreasing=TRUE),drop=FALSE]
  
  return(proper.beta.norm)
}


######################################
######################################

penalty.list = c('LASSO','SCAD','MCP')
penalty = penalty.list[3]
print(Ycols)

proper.beta.norm.list = sapply(Ycols,function(Ycol){
  Xcols = cols.all[!(cols.all %in% c(Ycol))]
  # load model
  load(paste0(result.path0,'NHANES_density_',Ycol,'_',penalty,'.Rdata'))
  
  # check variable selection (nonzero norm)
  proper.beta.norm = compute.proper.beta.norm(model,Xcols)
  
  proper.beta.norm
})

for (Ycol in Ycols){
  print(paste0('Ycol: ',Ycol,', penalty: ',penalty))
  print(proper.beta.norm.list[[Ycol]])
}





######################################
######################################


# For preliminary setting for plotting
pacman::p_load(RColorBrewer)
col.pal = RColorBrewer::brewer.pal(7, "Spectral")
new.color=colorRampPalette(col.pal)

n.grid = 2
n.t = ncol(PA)
t.all = seq(log10.cpm.min,log10.cpm.max,length.out=n.t)


# Diabetes
Ycol = Ycols[1]
proper.beta.norm = proper.beta.norm.list[[Ycol]]
idx = proper.beta.norm['ind','PA']

load(paste0(result.path0,'NHANES_density_',Ycol,'_',penalty,'.Rdata'))
beta = model$beta.tensor[[idx]]
beta.new1 = beta$element2 %*% beta$element1
beta.new1 = inv.clr.BayesHilbert(beta.new1,normalizer=normalizer)
beta.new1 = as.vector(beta.new1)

plot1 = qplot(x=t.all,y=beta.new1,color='',group='',geom='line')+
  labs(color='')+
  scale_x_continuous(expression(log[10]*'(cpm)'))+
  scale_y_continuous('PA Density')+
  ggtitle('Diabetes')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='')



# CHF
Ycol = Ycols[2]
proper.beta.norm = proper.beta.norm.list[[Ycol]]
idx = proper.beta.norm['ind','PA']

load(paste0(result.path0,'NHANES_density_',Ycol,'_',penalty,'.Rdata'))
beta = model$beta.tensor[[idx]]
beta.new2 = beta$element2 %*% beta$element1
beta.new2 = inv.clr.BayesHilbert(beta.new2,normalizer=normalizer)
beta.new2 = as.vector(beta.new2)

plot2 = qplot(x=t.all,y=beta.new2,color='',group='',geom='line')+
  labs(color='')+
  scale_x_continuous(expression(log[10]*'(cpm)'))+
  scale_y_continuous('PA Density')+
  ggtitle('CHF')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='')




# CHD
Ycol = Ycols[3]
print(Ycol)

proper.beta.norm = proper.beta.norm.list[[Ycol]]
idx = proper.beta.norm['ind','PA']

load(paste0(result.path0,'NHANES_density_',Ycol,'_',penalty,'.Rdata'))
beta = model$beta.tensor[[idx]]
beta.new3 = beta$element2 %*% beta$element1
beta.new3 = inv.clr.BayesHilbert(beta.new3,normalizer=normalizer)
beta.new3 = as.vector(beta.new3)

plot3 = qplot(x=t.all,y=beta.new3,color='',group='',geom='line')+
  labs(color='')+
  scale_x_continuous(expression(log[10]*'(cpm)'))+
  scale_y_continuous('PA Density')+
  ggtitle('CHD')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='')



# Stroke
Ycol = Ycols[5]
print(Ycol)

proper.beta.norm = proper.beta.norm.list[[Ycol]]
idx = proper.beta.norm['ind','PA']

load(paste0(result.path0,'NHANES_density_',Ycol,'_',penalty,'.Rdata'))
beta = model$beta.tensor[[idx]]
beta.new4 = beta$element2 %*% beta$element1
beta.new4 = inv.clr.BayesHilbert(beta.new4,normalizer=normalizer)
beta.new4 = as.vector(beta.new4)

plot4 = qplot(x=t.all,y=beta.new4,color='',group='',geom='line')+
  labs(color='')+
  scale_x_continuous(expression(log[10]*'(cpm)'))+
  scale_y_continuous('PA Density')+
  ggtitle('Stroke')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.position='')


plot_ = ggdraw()+
  draw_plot(plot1,x=0.0,y=0.5,width=0.5,height=0.5)+
  draw_plot(plot2,x=0.5,y=0.5,width=0.5,height=0.5)+
  draw_plot(plot3,x=0.0,y=0.0,width=0.5,height=0.5)+
  draw_plot(plot4,x=0.5,y=0.0,width=0.5,height=0.5)

plot_
ggsave(paste0(result.path0,'NHANES.pdf'),plot_,width=10,height=4.5)




