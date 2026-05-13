max.cv.iter=20
cv.threshold=1e-10
cv.const=2
phi0=1e-4
c.phi=1.1
max.iter=500
threshold=1e-10




Xorg = data$X
Yorg = data$Y
cv.type = 'AIC'


# check validility of inputs
Check.penalty(penalty)
Check.kernel.QM(kernel)
Check.cv.type(cv.type)

if ((penalty=='SCAD') & (gamma<2)){
  gamma = 3.7
} else if ((penalty=='MCP') & (gamma<1)){
  gamma = 3
}

# define basic parameters
n = nrow(Yorg)
p = Xorg[['p']]
if(is.null(lambda.list)){lambda.list=c(0)}
Ymu = 0
Yspace = 'Euclid'

# PCA for X
pca = PCA.manifold.list(Xorg)
X = predict(pca,Xorg)
Xdims = sapply(X,ncol)
if(is.null(Xdim.max.list)){Xdim.max.list = c(max(sapply(X,ncol)))}

# compute the default value of h
# In the CV function, the default value of h is computed in the 'QM_CV' function.
if (is.null(h)){h = c(-1.0)}



################################################
################################################
################################################


source('./__sources/R/utils.R')
sourceCpp('./__sources/cpp/QM2.cpp')
sourceCpp('./__sources/cpp/QM_CV2.cpp')


lambda.list = seq(0.005,0.1,length.out=20)[-(1:4)][1:10]
Xdim.max.list = 1:6

lambda.list = c(0.035,0.04)
Xdim.max.list = c(1,6)

# Use QM_CV function to obtain the optimal parameters
result = QM_CV2(X,Yorg,lambda.list,Xdim.max.list,cv.type,tau,h,kernel,penalty,gamma,cv.const,max.cv.iter,cv.threshold)

parameter.list = result$parameter.list[which(rowMeans(result$parameter.list)!=0),,drop=FALSE]
colnames(parameter.list) = c('lambda','Xdim.max')
loss.list = result$loss.list[-which(sapply(result$loss.list,is.null))]


# apply QM with the optimal parameters
opt.start.time = Sys.time()

opt.lambda = result$opt.lambda
opt.Xdim.max = result$opt.Xdim.max

opt.lambda = 0.04
opt.Xdim.max = 6

sourceCpp('./__sources/cpp/QM2.cpp')

object = QM_each2(X,Yorg,opt.lambda,opt.Xdim.max,tau,h,kernel,penalty,gamma,phi0,c.phi,max.iter,threshold)

get_loss_CV_QM2(X, Yorg, opt.lambda, opt.Xdim.max, cv.type, tau, h, kernel, penalty, gamma, cv.const)







