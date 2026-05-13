

# source('./__sources/R/utils.R')
# sourceCpp('./__sources/cpp/QM.cpp')
# sourceCpp('./__sources/cpp/QM_CV.cpp')
# sourceCpp('./__sources/cpp/QM_GCV.cpp')
# sourceCpp('./__sources/cpp/QM_kfold.cpp')
# source('./__sources/R/QM.R')
# source('./__sources/R/QM_oracle.R')
# source('./__sources/R/QM_CV.R')
# source('./__sources/R/QM_GCV.R')
# source('./__sources/R/QM_kfold.R')
# source('./__sources/R/QM_oracle_CV.R')



lambda.list = 0.01 * seq(2,10,2)
Xdim.max.list = c(1,2,3,4)

lambda.list = seq(0.01,0.1,length.out=19)
Xdim.max.list = 1:6


# tau = 0.2
h = NULL
c.h = 1.0
penalty = 'MCP'
gamma = 3.7

cv.const = 2
kfold = 5
seed = 3

kernel = 'Gaussian'
# kernel = 'Uniform'
# kernel = 'Epanechnikov'
# kernel = 'Triangular'


model0 = LM.CV(data$X,data$Y,'Euclid','AIC',penalty,gamma,lambda.list,Xdim.max.list,c(100,200),cv.const=cv.const)
Yhat0 = predict(model0,datatest$X)

model1 = QM.CV(data$X,data$Y,tau,h,kernel,'AIC',penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const,c.h=c.h)
Yhat1 = predict(model1,datatest$X)

model2 = QM.CV(data$X,data$Y,tau,h,kernel,'BIC',penalty,gamma,lambda.list,Xdim.max.list,cv.const=cv.const,c.h=c.h)
Yhat2 = predict(model2,datatest$X)

model3 = QM.GCV(data$X,data$Y,datanew$X,datanew$Y,tau,h,kernel,penalty,gamma,lambda.list,Xdim.max.list,c.h=c.h)
Yhat3 = predict(model3,datatest$X)

# model4 = QM.kfold(data$X,data$Y,kfold,seed,tau,h,kernel,penalty,gamma,lambda.list,Xdim.max.list,c.h=c.h)
model4 = model3
Yhat4 = predict(model4,datatest$X)

model5 = QM.oracle.CV(data$X,data$Y,tau,h,kernel,proper.indices,'AIC',Xdim.max.list,cv.const=cv.const,c.h=c.h)
Yhat5 = predict(model5,datatest$X)

model6 = QM.oracle.CV(data$X,data$Y,tau,h,kernel,proper.indices,'BIC',Xdim.max.list,cv.const=cv.const,c.h=c.h)
Yhat6 = predict(model6,datatest$X)

model7 = QM.oracle.GCV(data$X,data$Y,datanew$X,datanew$Y,proper.indices,tau,h,kernel,Xdim.max.list,c.h=c.h)
Yhat7 = predict(model7,datatest$X)

# model8 = QM.oracle.kfold(data$X,data$Y,proper.indices,kfold,seed,tau,h,kernel,Xdim.max.list,c.h=c.h)
model8 = model7
Yhat8 = predict(model8,datatest$X)


model0$runtime
model1$runtime;model2$runtime;model3$runtime;model4$runtime
model5$runtime;model6$runtime;model7$runtime;model8$runtime
model0$parameter.list
model1$parameter.list;model2$parameter.list;model3$parameter.list;model4$parameter.list
model5$Xdim.max;model6$Xdim.max;model7$Xdim.max;model8$Xdim.max
proper.indices;model0$proper.indices
model1$proper.indices;model2$proper.indices;model3$proper.indices;model4$proper.indices


model0$Ymu;
model1$beta0;model2$beta0;model3$beta0;model4$beta0
model5$beta0;model6$beta0;model7$beta0;model8$beta0
model0$beta.norm[model0$proper.indices]
model1$beta.norm[model1$proper.indices]
model2$beta.norm[model2$proper.indices]
model3$beta.norm[model3$proper.indices]
model4$beta.norm[model4$proper.indices]
model5$beta.norm[model5$proper.indices]
model6$beta.norm[model6$proper.indices]
model7$beta.norm[model7$proper.indices]
model8$beta.norm[model8$proper.indices]


mean(dist.manifold(datatest$Xbeta,Yhat0,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat1,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat4,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat5,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat6,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat7,Yspace)^2)
mean(dist.manifold(datatest$Xbeta,Yhat8,Yspace)^2)


median(dist.manifold(datatest$Xbeta,Yhat0,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat1,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat2,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat3,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat4,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat5,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat6,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat7,Yspace))
median(dist.manifold(datatest$Xbeta,Yhat8,Yspace))


quantile(datatest$error,tau)
quantile(datatest$Y-Yhat0,tau) 
quantile(datatest$Y-Yhat1,tau)
quantile(datatest$Y-Yhat2,tau)
quantile(datatest$Y-Yhat3,tau)
quantile(datatest$Y-Yhat4,tau)
quantile(datatest$Y-Yhat5,tau)
quantile(datatest$Y-Yhat6,tau)
quantile(datatest$Y-Yhat7,tau)
quantile(datatest$Y-Yhat8,tau)









