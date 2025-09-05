
# sourceCpp('GLM.cpp')

lambda = 0.02
Xdim.max = 6
R = 100

penalty = 'MCP'
gamma = 3.7
phi = 1
phi2 = 0.2



object1 = GLM(data$X,data$Y,link,penalty,gamma,lambda,Xdim.max,R,phi)
Yhat = predict(object1,data$X)

object2 = GLM(data$X,data$Y,link,penalty,gamma,lambda,Xdim.max,R,phi2)
Yhat2 = predict(object2,data$X)

object3 = GLM.oracle(data$X,data$Y,link,proper.indices,Xdim.max)
Yhat3 = predict(object3,data$X)


object1$runtime;object2$runtime;object3$runtime
proper.indices;object1$proper.indices;object2$proper.indices;object3$proper.indices
object1$beta0;object2$beta0;object3$beta0
object1$beta.norm[object1$proper.indices]
object2$beta.norm[object2$proper.indices]
object3$beta.norm[object3$proper.indices]



mean(dist.manifold(data$Ymu,Yhat,Yspace)^2)
mean(dist.manifold(data$Ymu,Yhat2,Yspace)^2)
mean(dist.manifold(data$Ymu,Yhat3,Yspace)^2)
mean(dist.manifold(Yhat,Yhat2,Yspace)^2)
mean(dist.manifold(Yhat,Yhat3,Yspace)^2)






