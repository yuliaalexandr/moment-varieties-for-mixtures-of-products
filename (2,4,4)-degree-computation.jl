using HomotopyContinuation

@var t110, t111, t112, t113, t122, t123, t131, t132, t133, t141, t142, t143, t210, t211, t212, t213, t222, t223, t231, t232, t233, t241, t242, t243

t = [t110, t111, t112, t113, t122, t123, t131, t132, t133, t141, t142, t143, t210, t211, t212, t213, t222, t223, t231, t232, t233, t241, t242, t243]

@var m3100, m3010, m3001, m2200, m2110, m2101, m2020, m2011, m2002, m1300, m1210, m1201, m1120, m1111, m1102, m1030, m1021, m1012, m1003, m0310, m0301, m0220, m0211, m0202, m0130, m0121, m0112, m0103, m0031, m0022, m0013

m = [m3100, m3010, m3001, m2200, m2110, m2101, m2020, m2011, m2002, m1300, m1210, m1201, m1120, m1111, m1102, m1030, m1021, m1012, m1003, m0310, m0301, m0220, m0211, m0202, m0130, m0121, m0112, m0103, m0031, m0022, m0013]

par = [t113+t213, t113*t131+t213*t231, t113*t141+t213*t241, t112*t122+t212*t222, t112*t131+t212*t231, t112*t141+t212*t241, t112*t132+t212*t232, t112*t131*t141+t212*t231*t241, t112*t142+t212*t242, t111*t123+t211*t223, t111*t122*t131+t211*t222*t231, t111*t122*t141+t211*t222*t241, t111*t132+t211*t232, t111*t131*t141+t211*t231*t241, t111*t142+t211*t242, t111*t133+t211*t233, t111*t132*t141+t211*t232*t241, t111*t131*t142+t211*t231*t242, t111*t143+t211*t243, t110*t123*t131+t210*t223*t231, t110*t123*t141+t210*t223*t241, t110*t122*t132+t210*t222*t232, t110*t122*t131*t141+t210*t222*t231*t241, t110*t122*t142+t210*t222*t242, t110*t133+t210*t233, t110*t132*t141+t210*t232*t241, t110*t131*t142+t210*t231*t242, t110*t143+t210*t243, t110*t133*t141+t210*t233*t241, t110*t132*t142+t210*t232*t242, t110*t131*t143+t210*t231*t243]

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
