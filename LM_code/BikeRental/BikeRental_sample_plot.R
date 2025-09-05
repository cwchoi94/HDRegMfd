
library(HDRegMfd)
library(pacman)
library(ggplot2)
library(gridExtra)
library(cowplot)


# log transformation
# Xbeta = RieExp.manifold(model$Ymu,Xbeta,model$beta.tensor[[idx]]$space2)
# 
# y label
# scale_y_continuous(expression(log[10]*'(1+Bike Rental)'))+

# must include normalizer=(Y.max-Y.min)^(-1) in "inv.clr.density" function.
data.path0 = './LM_code/BikeRental/data/'
result.path0 = './LM_code/BikeRental/result/'
load(paste0(data.path0,'BikeRental.Rdata'))

n = nrow(Y)

matplot(t(10^(Y)), type = 'l', lty = 1, col = rainbow(n),
        xlab = "Hour", ylab = "Bike rental count",
        main = "Bike sharing data")


