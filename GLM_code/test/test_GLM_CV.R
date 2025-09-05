
# sourceCpp('GLM_CV.cpp')

lambda.list = c(0.1,0.2,0.3,0.4)
Xdim.max.list = c(1,2,3,4)
R.list = c(100,200)

# lambda.list = seq(0.005,0.2,length.out=40)
# Xdim.max.list = 1:6
# R.list = c(100,200)

penalty = 'MCP'
gamma = 3.7
phi = 1
phi2 = 3
phi3 = 0.5

phi = 0.5

cv.const = 2


model1 = GLM.CV(data$X,data$Y,link,'AIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,phi,cv.const=cv.const)
Yhat = predict(model1,datatest$X)

model2 = GLM.CV(data$X,data$Y,link,'BIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,phi2,cv.const=cv.const)
Yhat2 = predict(model2,datatest$X)

model3 = GLM.CV(data$X,data$Y,link,'ABIC',penalty,gamma,lambda.list,Xdim.max.list,R.list,phi3,cv.const=cv.const)
Yhat3 = predict(model3,datatest$X)

model4 = GLM.oracle.CV(data$X,data$Y,link,proper.indices,'AIC',Xdim.max.list,cv.const=cv.const)
Yhat4 = predict(model4,datatest$X)

model5 = GLM.oracle.CV(data$X,data$Y,link,proper.indices,'BIC',Xdim.max.list,cv.const=cv.const)
Yhat5 = predict(model5,datatest$X)

model6 = GLM.oracle.CV(data$X,data$Y,link,proper.indices,'ABIC',Xdim.max.list,cv.const=cv.const)
Yhat6 = predict(model6,datatest$X)



model1 = GLM.CV(data$X,data$Y,link,'AIC','LASSO',gamma,lambda.list,Xdim.max.list,R.list,phi,cv.const=cv.const)
Yhat = predict(model1,datatest$X)

model2 = GLM.CV(data$X,data$Y,link,'AIC','SCAD',gamma,lambda.list,Xdim.max.list,R.list,phi,cv.const=cv.const)
Yhat2 = predict(model1,datatest$X)

model3 = GLM.CV(data$X,data$Y,link,'AIC','MCP',gamma,lambda.list,Xdim.max.list,R.list,phi,cv.const=cv.const)
Yhat3 = predict(model1,datatest$X)




model1$runtime;model2$runtime;model3$runtime;model4$runtime;model5$runtime
model1$parameter.list;model2$parameter.list;model3$parameter.list
model4$Xdim.max;model5$Xdim.max;model6$Xdim.max
proper.indices;model1$proper.indices;model2$proper.indices;model3$proper.indices

model1$beta0;model2$beta0;model3$beta0;model4$beta0
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
mean(dist.manifold(Yhat,Yhat2,Yspace)^2)
mean(dist.manifold(Yhat,Yhat3,Yspace)^2)
mean(dist.manifold(Yhat,Yhat4,Yspace)^2)




(nrow(model1$parameter.list)-2)%/%3




