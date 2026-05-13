
# source('initialize.R')

# source('./__sources/R/utils.R')
# sourceCpp('./__sources/cpp/QM.cpp')
# source('./__sources/R/QM.R')
# source('./__sources/R/QM_oracle.R')

lambda = 0.01
Xdim.max = 4

penalty = 'MCP'
gamma = 3.7
h = NULL



kernel = 'Gaussian'
c.h = 1.0
phi0 = 1e-4
c.phi = 1.1
max.iter = 500
threshold = 1e-8

Xorg = data$X
Yorg = data$Y


object0 = LM(data$X,data$Y,'Euclid',penalty,gamma,lambda,Xdim.max,1e5)
Yhat0 = predict(object0,datatest$X)

object1 = QM(data$X,data$Y,tau,h,'Gaussian',penalty,gamma,lambda,Xdim.max,c.h)
Yhat1 = predict(object1,datatest$X)

object2 = QM(data$X,data$Y,tau,h,'Uniform',penalty,gamma,lambda,Xdim.max,c.h)
Yhat2 = predict(object2,datatest$X)

object3 = QM(data$X,data$Y,tau,h,'Epanechnikov',penalty,gamma,lambda,Xdim.max,c.h)
Yhat3 = predict(object3,datatest$X)

object4 = QM(data$X,data$Y,tau,h,'Triangular',penalty,gamma,lambda,Xdim.max,c.h)
Yhat4 = predict(object4,datatest$X)

object5 = QM.oracle(data$X,data$Y,proper.indices,tau,h,'Gaussian',Xdim.max,c.h)
Yhat5 = predict(object5,datatest$X)

object6 = QM.oracle(data$X,data$Y,proper.indices,tau,h,'Uniform',Xdim.max,c.h)
Yhat6 = predict(object6,datatest$X)

object7 = QM.oracle(data$X,data$Y,proper.indices,tau,h,'Epanechnikov',Xdim.max,c.h)
Yhat7 = predict(object7,datatest$X)

object8 = QM.oracle(data$X,data$Y,proper.indices,tau,h,'Triangular',Xdim.max,c.h)
Yhat8 = predict(object8,datatest$X)




object0$runtime;
object1$runtime;object2$runtime;object3$runtime;object4$runtime
object5$runtime;object6$runtime;object7$runtime;object8$runtime
proper.indices;object0$proper.indices
object1$proper.indices;object2$proper.indices;object3$proper.indices;object4$proper.indices
object5$proper.indices;object6$proper.indices;object7$proper.indices;object8$proper.indices
object0$beta.norm[object0$proper.indices]
object1$beta.norm[object1$proper.indices]
object2$beta.norm[object2$proper.indices]
object3$beta.norm[object3$proper.indices]
object4$beta.norm[object4$proper.indices]
object5$beta.norm[object5$proper.indices]
object6$beta.norm[object6$proper.indices]
object7$beta.norm[object7$proper.indices]
object8$beta.norm[object8$proper.indices]


quantile(data$error,tau);object0$Ymu
object1$beta0;object2$beta0;object3$beta0;object4$beta0
object5$beta0;object6$beta0;object7$beta0;object8$beta0


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





