

n = 200
n2 = 100

p = 9
dims = c(3,2,3,4,4,3,30,20,40)
spaces = c(rep('Euclid',2),rep('simplex',2),rep('sphere',2),'functional','Wasserstein','BayesHilbert')

mu.list = list()
mu.list[[1]] = rep(0,dims[1])
mu.list[[2]] = rep(0,dims[2])
mu.list[[3]] = inv.clr.simplex(rep(0,dims[3]+1))
mu.list[[4]] = inv.clr.simplex(rep(0,dims[4]+1))
mu.list[[5]] = c(rep(0,dims[5]),1)
mu.list[[6]] = c(rep(0,dims[6]),1)
mu.list[[7]] = rep(0,100)
mu.list[[8]] = rep(0,150)
mu.list[[9]] = inv.clr.BayesHilbert(rep(0,120))

Xdata = covariates.generate(n,spaces,dims,mu.list,0.5)
Xdatanew = covariates.generate(n2,spaces,dims,mu.list,0.5)


pca.list = PCA.RiemannSpace.list(Xdata)
scores = predict(pca.list,Xdata)



