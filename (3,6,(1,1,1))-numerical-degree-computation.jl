### (3,6,(1,1,1))
#### identifiable (6:1)
using HomotopyContinuation
@var t111, t131, t141, t151, t161, t211, t231, t241, t251, t261, t311, t331, t341, t351,t361, t110, t210, t310
t=[t111, t131, t141, t151, t161, t211, t231, t241, t251, t261, t311, t331, t341, t351, t361, t110, t210, t310]
@var m000111, m010110, m011001, m110100, m100101, m001110, m010101, m101100, m110010, m100011, m001101, m010011, m101010, m110001, m011100, m001011, m111000, m101001, m100110, m011010
m=[m000111, m010110, m011001, m110100, m100101, m001110, m010101, m101100, m110010, m100011, m001101, m010011, m101010, m110001, m011100, m001011, m111000, m101001, m100110, m011010]

par=[t110*t141*t151*t161+t210*t241*t251*t261+t310*t341*t351*t361,t110*t141*t151+t210*t241*t251+t310*t341*t351, t110*t131*t161+t210*t231*t261+t310*t331*t361, t111*t141+t211*t241+t311*t341, t111*t141*t161+t211*t241*t261+t311*t341*t361, t110*t131*t141*t151+t210*t231*t241*t251+t310*t331*t341*t351, t110*t141*t161+t210*t241*t261+t310*t341*t361, t111*t131*t141+t211*t231*t241+t311*t331*t341, t111*t151+t211*t251+t311*t351, t111*t151*t161+t211*t251*t261+t311*t351*t361, t110*t131*t141*t161+t210*t231*t241*t261+t310*t331*t341*t361,t110*t151*t161+t210*t251*t261+t310*t351*t361, t111*t131*t151+t211*t231*t251+t311*t331*t351, t111*t161+t211*t261+t311*t361,t110*t131*t141+t210*t231*t241+t310*t331*t341,t110*t131*t151*t161+t210*t231*t251*t261+t310*t331*t351*t361,t111*t131+t211*t231+t311*t331,t111*t131*t161+t211*t231*t261+t311*t331*t361,t111*t141*t151+t211*t241*t251+t311*t341*t351,t110*t131*t151+t210*t231*t251+t310*t331*t351]

@var a[1:length(t),1:length(m)+1]

# compute a point on the image, and a length(t)-1 codimensional linear space containing it.

seed_t = randn(ComplexF64,length(t))
seed_m = evaluate(par,t=>seed_t)

using LinearAlgebra
N = nullspace([transpose(seed_m) 1])
NN = randn(ComplexF64,length(t),size(N,2))*transpose(N)
linspace = NN*[m;1]

# Set up the parametrized system. I'll use a moving linear space: the coefficients are parameters
# I'll set m[1] = 1 for dehomogenization

polsys = System([m-par;a*[m;1]],variables = [t;m], parameters = a[:])

MR=monodromy_solve(polsys, [seed_t;seed_m], NN[:])