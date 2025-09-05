


# sourceCpp('./__sources/cpp/LM.cpp')

lambda = 0.1
Xdim.max = 100
R = 1

penalty = 'MCP'
gamma = 3.7
phi = 1
phi2 = 5


object1 = LM(data$X,data$Y,Yspace,penalty,gamma,lambda,Xdim.max,R,phi)
Yhat = predict(object1,data$X)

object2 = LM(data$X,data$Y,Yspace,penalty,gamma,lambda,Xdim.max,R,phi2)
Yhat2 = predict(object2,data$X)

object3 = LM.oracle(data$X,data$Y,Yspace,proper.indices,Xdim.max)
Yhat3 = predict(object3,data$X)


object1$runtime;object2$runtime;object3$runtime
proper.indices;object1$proper.indices;object2$proper.indices;object3$proper.indices
object1$beta.norm[object1$proper.indices]
object2$beta.norm[object2$proper.indices]
object3$beta.norm[object3$proper.indices]


mean(dist.manifold(data$Y,Yhat,Yspace)^2)
mean(dist.manifold(data$Y,Yhat2,Yspace)^2)
mean(dist.manifold(data$Y,Yhat3,Yspace)^2)
mean(dist.manifold(Yhat,Yhat2,Yspace)^2)
mean(dist.manifold(Yhat,Yhat3,Yspace)^2)






