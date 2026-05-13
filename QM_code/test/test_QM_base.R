

sourceCpp('./__sources/cpp/QM_base.cpp')




u_org = matrix(0.1+0.1*1:6,3,2)
u_org = matrix(rnorm(150,0,2),10,15)
tau = 0.2
h = 0.5

u = u_org/h


sQRloss(u_org,tau,h,"Gaussian")     - (tau-1/2)*u_org - h / 2 * (sqrt(2/pi) * exp(-(u^2) / 2) + u * (1-2*pnorm(-u)))
sQRloss(u_org,tau,h,"Uniform")      - (tau-1/2)*u_org - h / 2 * ((u^2 + 1)/2 * (abs(u)<1) + abs(u) * (abs(u)>=1))
sQRloss(u_org,tau,h,"Epanechnikov") - (tau-1/2)*u_org - h / 2 * ((3/4*u^2 - u^4 / 8 + 3/8) * (abs(u)<1) + abs(u) * (abs(u)>=1))
sQRloss(u_org,tau,h,"Triangular")   - (tau-1/2)*u_org - h / 2 * ((abs(u)^2 - abs(u)^3 / 3 + 1/3) * (abs(u)<1) + abs(u) * (abs(u)>=1))


sQRloss_diff(u_org,tau,h,"Gaussian")     - (tau-1/2) - 1 / 2 * (-sqrt(2/pi) * u * exp(-(u^2) / 2) + (1-2*pnorm(-u)) + 2 * u * dnorm(-u))
sQRloss_diff(u_org,tau,h,"Uniform")      - (tau-1/2) - 1 / 2 * ( u * (abs(u)<1) - (u<=-1) + (u>=1))
sQRloss_diff(u_org,tau,h,"Epanechnikov") - (tau-1/2) - 1 / 2 * ( (3/2*u - u^3 / 4) * (abs(u)<1) - (u<=-1) + (u>=1))
sQRloss_diff(u_org,tau,h,"Triangular")   - (tau-1/2) - 1 / 2 * ( (2*u + u^2) * (-1<u) * (u<=0) + (2*u-u^2) * (0<u) * (u<1) - (u<=-1) + (u>=1))



u_org = matrix(seq(-10,10,length.out=100),ncol=1)
tau = 1.1
h = 0.05

plot(sQRloss(u_org,tau,h,"Gaussian"),lty=1)
lines(sQRloss(u_org,tau,h,"Uniform"),lty=2)
lines(sQRloss(u_org,tau,h,"Epanechnikov"),lty=3)
lines(sQRloss(u_org,tau,h,"Triangular"),lty=4)







