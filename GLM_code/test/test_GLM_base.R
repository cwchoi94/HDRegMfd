



a = matrix(0.1+0.1*1:6,3,2)

Link(a,'binomial') - log(a/(1-a))
Link(a,'poisson') - log(a)
Link(a,'exponential') - (-1/a)

Inv_Link(Link(a,'binomial'),'binomial') - a
Inv_Link(Link(a,'poisson'),'poisson') - a
Inv_Link(Link(a,'exponential'),'exponential') - a


Psi(a,'binomial') - (log(1+exp(a)))
Psi(a,'poisson') - (exp(a))
Psi(-a,'exponential') - (-log(a))


Psi_1d(a,'binomial') - exp(a)/(1+exp(a))
Psi_1d(a,'poisson') - exp(a)
Psi_1d(a,'exponential') - (-1/a)

Psi_2d(a,'binomial') - exp(a)/(1+exp(a))^2
Psi_2d(a,'poisson') - exp(a)
Psi_2d(a,'exponential') - (1/a^2)





