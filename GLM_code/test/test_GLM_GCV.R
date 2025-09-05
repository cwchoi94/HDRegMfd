
# sourceCpp('GLM_GCV.cpp')

lambda.list = c(0.1,0.2,0.3,0.4)
Xdim.max.list = c(1,2,3,4)
R.list = c(1,2,3,4)

lambda.list = seq(0.005,0.2,length.out=40)
Xdim.max.list = 1:6
R.list = c(100,200)

penalty = 'SCAD'
gamma = 3.7
phi = 1
phi2 = 3
phi3 = 0.5


model1 = GLM.GCV(data$X,data$Y,datanew$X,datanew$Y,link,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi)
Yhat = predict(model1,datatest$X)

model2 = GLM.GCV(data$X,data$Y,datanew$X,datanew$Y,link,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi2)
Yhat2 = predict(model2,datatest$X)

model3 = GLM.GCV(data$X,data$Y,datanew$X,datanew$Y,link,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi3)
Yhat3 = predict(model3,datatest$X)

model4 = GLM.oracle.GCV(data$X,data$Y,datanew$X,datanew$Y,link,proper.indices,Xdim.max.list)
Yhat4 = predict(model4,datatest$X)


model1$runtime;model2$runtime;model3$runtime;model4$runtime
model1$parameter.list;model2$parameter.list;model3$parameter.list
proper.indices;model1$proper.indices;model2$proper.indices;model3$proper.indices;model4$proper.indices

model1$beta0;model2$beta0;model3$beta0;model4$beta0
model1$beta.norm[model1$proper.indices]
model2$beta.norm[model2$proper.indices]
model3$beta.norm[model3$proper.indices]
model4$beta.norm[model4$proper.indices]

mean(dist.manifold(datatest$Ymu,Yhat,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$Ymu,Yhat4,Yspace)^2)
mean(dist.manifold(Yhat,Yhat2,Yspace)^2)
mean(dist.manifold(Yhat,Yhat3,Yspace)^2)
mean(dist.manifold(Yhat,Yhat4,Yspace)^2)








