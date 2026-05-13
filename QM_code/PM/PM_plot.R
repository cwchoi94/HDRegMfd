
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)


data.path0 = './QM_code/PM/data/'
result.path0 = './QM_code/PM/result/'
load(paste0(data.path0,'PM2.5.Rdata'))

dim(Y)
Xspaces = Xdata[['spaces']]

Xspaces
sapply(1:Xdata[['p']],function(j){dim(Xdata[[j]])})


compute.proper.beta.norm = function(tau,penalty){
  load(paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))
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

tau.list = seq(0.1,0.9,0.1)
tau.list = c(0.25,0.5,0.75)
penalty.list = c('LASSO','SCAD','MCP')


# check variable selection (nonzero norm)
penalty = penalty.list[2]

proper.beta.norm.list = list()
for (k in 1:length(tau.list)){
  tau = tau.list[k]
  proper.beta.norm = compute.proper.beta.norm(tau,penalty)
  
  print(c(penalty,tau,ncol(proper.beta.norm),Xdata$p))
  print(proper.beta.norm)
  
  proper.beta.norm.list[[k]] = proper.beta.norm
}


print(c(penalty,tau))
proper.beta.norm
Xspaces[proper.beta.norm['ind',]]

dim(proper.beta.norm)
Xdata[['p']]



# load model
tau = tau.list[3]
penalty = penalty.list[2]

load(paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))

Xdim.max = model$Xdim.max
score.list = predict(model$pca,Xdata)


######################################
######################################


# season
beta.list = c()
for (k in 1:length(tau.list)){
  tau = tau.list[k]
  proper.beta.norm = proper.beta.norm.list[[k]]
  
  for (season in c('summer','fall','winter')){
    if (!(season %in% colnames(proper.beta.norm))){
      beta.list = c(beta.list,0)
      next
    }
    ind = proper.beta.norm['ind',season]
    load(paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))
    beta.list = c(beta.list,model$beta.tensor[[ind]]$element2)
  }
}
season.effect = matrix(beta.list,3,3,byrow=TRUE)
rownames(season.effect) = tau.list
colnames(season.effect) = c('summer','fall','winter')
season.effect


# weather
beta.list = c()
for (k in 1:length(tau.list)){
  tau = tau.list[k]
  proper.beta.norm = proper.beta.norm.list[[k]]
  
  for (weather in c('cloudy or misty','snowy or rainy')){
    if (!(weather %in% colnames(proper.beta.norm))){
      beta.list = c(beta.list,0)
      next
    }
    ind = proper.beta.norm['ind',weather]
    load(paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))
    beta.list = c(beta.list,model$beta.tensor[[ind]]$element2)
    print(paste(weather,tau,model$beta.tensor[[ind]]$element2))
  }
}
weather.effect = matrix(beta.list,3,2,byrow=TRUE)
rownames(weather.effect) = tau.list
colnames(weather.effect) = c('cloudy or misty','snowy or rainy')
weather.effect



# For preliminary setting for plotting
pacman::p_load(RColorBrewer)
col.pal = RColorBrewer::brewer.pal(7, "Spectral")
new.color=colorRampPalette(col.pal)

Y.min = 1
Y.max = 24

n.grid = 40
n.t = ncol(model$beta)



# wind
idx = 12
Y.min = 1
Y.max = 24

beta.list = c()
for (k in 1:length(tau.list)){
  tau = tau.list[k]
  proper.beta.norm = proper.beta.norm.list[[k]]
  
  load(paste0(result.path0,'PM2.5_',tau,'_',penalty,'.Rdata'))
  beta = model$beta.tensor[[idx]]
  beta.new = beta$element2 %*% beta$element1
  beta.new = as.vector(beta.new)
  
  beta.list = rbind(beta.list,beta.new)
}
rownames(beta.list) = tau.list
df = data.frame(tau=as.factor(rep(tau.list,each=24)),
                X = rep(1:24,3),
                Xbeta=as.vector(t(beta.list)))


plot1 = ggplot(data=df)+aes(x=X,y=Xbeta,color=tau,group=tau)+geom_line()+
  labs(color=parse(text=expression(tau)))+
  scale_x_continuous('hour',breaks = ceiling(seq(Y.min,Y.max+1, by = 4)))+
  scale_y_continuous('')+
  ggtitle('Wind-speed')+
  theme_gray()+
  theme(plot.title=element_text(hjust=0.5),legend.title=element_text(size=16))

plot1
ggsave(paste0(result.path0,'wind-speed.pdf'),plot1,width=5,height=3)


season.effect = round(season.effect,4)
weather.effect = round(weather.effect,4)


write.csv(season.effect,paste0(result.path0,'season.csv'))
write.csv(weather.effect,paste0(result.path0,'weather.csv'))
season.effect
weather.effect



