

n = 200
n2 = 100

Xdata = list()
Xdatanew = list()

# real
Xdata[[1]] = matrix(rnorm(3*n),n,3)
Xdatanew[[1]] = matrix(rnorm(3*n2),n2,3)

Xdata[[2]] = matrix(rnorm(2*n),n,2)
Xdatanew[[2]] = matrix(rnorm(2*n2),n2,2)


# simplex
basis3 = basis.simplex(1:3)
A = matrix(rnorm(2*n),n,2)
Anew = matrix(rnorm(2*n2),n2,2)

Xdata[[3]] = RieExp.simplex(c(1/3,1/3,1/3),A %*% basis3)
Xdatanew[[3]] = RieExp.simplex(c(1/3,1/3,1/3),Anew %*% basis3)

basis4 = basis.simplex(1:4)
A = matrix(rnorm(3*n),n,3)
Anew = matrix(rnorm(3*n2),n2,3)

Xdata[[4]] = RieExp.simplex(c(1/4,1/4,1/4,1/4),A %*% basis4)
Xdatanew[[4]] = RieExp.simplex(c(1/4,1/4,1/4,1/4),Anew %*% basis4)


# sphere

mu5 = c(1/2,1/2,1/2,1/2)
basis5 = basis.sphere(mu5)
A = matrix(rnorm(3*n),n,3)
Anew = matrix(rnorm(3*n2),n2,3)

Xdata[[5]] = RieExp.sphere(mu5,A %*% basis5)
Xdatanew[[5]] = RieExp.sphere(mu5,Anew %*% basis5)


mu6 = c(1/2,1/2,1/sqrt(2))
basis6 = basis.sphere(mu6)
A = matrix(rnorm(2*n),n,2)
Anew = matrix(rnorm(2*n2),n2,2)

Xdata[[6]] = RieExp.sphere(mu6,A %*% basis6)
Xdatanew[[6]] = RieExp.sphere(mu6,Anew %*% basis6)


# functional

basis7 = basis.functional(1:100,30)
A = matrix(rnorm(30*n),n,30)
Anew = matrix(rnorm(30*n2),n2,30)

Xdata[[7]] = RieExp.functional(rep(0,100),A %*% basis7)
Xdatanew[[7]] = RieExp.functional(rep(0,100),Anew %*% basis7)


# Wasserstein

t = seq(0,1,length.out=102)[2:101]
mu8 = qnorm(t)

basis8 = basis.Wasserstein(mu8,40)
A = matrix(rnorm(40*n),n,40)
Anew = matrix(rnorm(40*n2),n2,40)

Xdata[[8]] = RieExp.Wasserstein(mu8,A %*% basis8)
Xdatanew[[8]] = RieExp.Wasserstein(mu8,Anew %*% basis8)


# BayesHilbert

mu9 = inv.clr.BayesHilbert(rep(0,100))[1,]

basis9 = basis.BayesHilbert(mu9,20)
A = matrix(rnorm(20*n),n,20)
Anew = matrix(rnorm(20*n2),n2,20)

Xdata[[9]] = RieExp.BayesHilbert(mu9,A %*% basis9)
Xdatanew[[9]] = RieExp.BayesHilbert(mu9,Anew %*% basis9)



Xdata[['p']] = 9
Xdatanew[['p']] = 9

Xdata[['spaces']] = c(rep('Euclid',2),rep('simplex',2),rep('sphere',2),'functional','Wasserstein','BayesHilbert')
Xdatanew[['spaces']] = c(rep('Euclid',2),rep('simplex',2),rep('sphere',2),'functional','Wasserstein','BayesHilbert')




pca.list = PCA.manifold.list(Xdata)
scores = predict(pca.list,Xdata)


pca.list2 = PCA.manifold.list(Xdata,alpha=1)
scores2 = predict(pca.list2,Xdata)

sapply(scores,ncol)
sapply(scores2,ncol)






