
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


time = seq(Y.min,Y.max,length.out=ncol(Y))
n = nrow(Y)


matplot(x=time,t(Y), type = 'l', lty = 1, col = rainbow(n),
        xlab = "Income", ylab = "Density",
        main = "Income distribution")


