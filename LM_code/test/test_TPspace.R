
space1 = 'functional'
space2 = 'simplex'

mu1 = rep(0,100)
mu2 = inv.clr.simplex(rep(0,5))[1,]

basis1 = basis.functional(mu1,3)
basis2 = basis.simplex(mu2)[1:3,]

element1 = basis1
element2 = basis2

z1 = make.tensor(element1,element2,space1,space2,mu1,mu2)

element1 = basis1[c(1,3,2,2),] + 2*basis1[c(2,1,2,3),]
element2 = basis2[c(1,3,2,3),]

z2 = make.tensor(element1,element2,space1,space2,mu1,mu2)

inner.tensor(z1,z2)
norm.tensor(z1)
norm.tensor(z2)
dist.tensor(z1,z2)


operator.tensor(z1,basis1[c(1,3),]) - basis2[c(1,3),]
operator.tensor(z2,basis1)[c(1,3),] - rbind(basis2[1,]+2*basis2[3,],2*basis2[1,]+3*basis2[2,]+basis2[3,],3*basis2[3,])[c(1,3),]






