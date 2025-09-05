

lambda.list = c(0.1,0.2,0.3,0.4)
Xdim.max.list = c(1,2,3,4)
R.list = c(100,200)

# lambda.list = seq(0.001,0.4,length.out=40)
# Xdim.max.list = 1:6
# R.list = c(100,200)

penalty = 'SCAD'
gamma = 3.7

cv.const = 2

model1 = LM.CV(data$X,data$Y,Yspace,'AIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
Yhat = predict(model1,datatest$X)

model2 = LM.CV(data$X,data$Y,Yspace,'BIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
Yhat2 = predict(model2,datatest$X)

model3 = LM.CV(data$X,data$Y,Yspace,'ABIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,cv.const=cv.const)
Yhat3 = predict(model3,datatest$X)

model4 = LM.oracle.CV(data$X,data$Y,Yspace,proper.indices,'AIC',Xdim.max.list,cv.const=cv.const)
Yhat4 = predict(model4,datatest$X)

model5 = LM.oracle.CV(data$X,data$Y,Yspace,proper.indices,'BIC',Xdim.max.list,cv.const=cv.const)
Yhat5 = predict(model5,datatest$X)

model6 = LM.oracle.CV(data$X,data$Y,Yspace,proper.indices,'ABIC',Xdim.max.list,cv.const=cv.const)
Yhat6 = predict(model6,datatest$X)


model1$runtime;model2$runtime;model3$runtime;model4$runtime;model5$runtime
model1$parameter.list;model2$parameter.list;model3$parameter.list
model4$Xdim.max;model5$Xdim.max;model6$Xdim.max
proper.indices;model1$proper.indices;model2$proper.indices;model3$proper.indices

model1$beta.norm[model1$proper.indices]
model2$beta.norm[model2$proper.indices]
model3$beta.norm[model3$proper.indices]
model4$beta.norm[model4$proper.indices]
model5$beta.norm[model5$proper.indices]
model6$beta.norm[model6$proper.indices]

mean(dist.manifold(datatest$Ymu,Yhat,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat4,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat5,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat6,Yspace)^2)

mean(dist.manifold(datatest$ExpXbeta,Yhat,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat4,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat5,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat6,Yspace)^2)








