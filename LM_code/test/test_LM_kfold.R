
# sourceCpp('./__sources/cpp/LM_kfold.cpp')

lambda.list = c(0.1,0.2,0.3,0.4)
Xdim.max.list = c(1,2,3,4)
R.list = c(1,2,3,4)

lambda.list = seq(0.001,0.4,length.out=40)
Xdim.max.list = 1:6
R.list = c(4,8,12,16,20)

penalty = 'SCAD'
gamma = 3.7
phi = 1
phi2 = 3
phi3 = 0.5
kfold = 5
seed = 2


model1 = LM.kfold(data$X,data$Y,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi)
Yhat = predict(model1,datatest$X)

model2 = LM.kfold(data$X,data$Y,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi2)
Yhat2 = predict(model2,datatest$X)

model3 = LM.kfold(data$X,data$Y,Yspace,kfold,seed,penalty,gamma,lambda.list,Xdim.max.list,R.list,phi3)
Yhat3 = predict(model3,datatest$X)

model4 = LM.oracle.kfold(data$X,data$Y,Yspace,proper.indices,kfold,seed,Xdim.max.list)
Yhat4 = predict(model4,datatest$X)


model1$runtime;model2$runtime;model3$runtime;model4$runtime
model1$parameter.list;model2$parameter.list;model3$parameter.list
proper.indices;model1$proper.indices;model2$proper.indices;model3$proper.indices;model4$proper.indices

model1$beta.norm[model1$proper.indices]
model2$beta.norm[model2$proper.indices]
model3$beta.norm[model3$proper.indices]
model4$beta.norm[model4$proper.indices]

mean(dist.manifold(datatest$Y,Yhat,Yspace)^2)
mean(dist.manifold(datatest$Y,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$Y,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$Y,Yhat4,Yspace)^2)

mean(dist.manifold(datatest$ExpXbeta,Yhat,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat2,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat3,Yspace)^2)
mean(dist.manifold(datatest$ExpXbeta,Yhat4,Yspace)^2)









