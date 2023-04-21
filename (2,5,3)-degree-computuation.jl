using HomotopyContinuation

@var t110, t111, t112, t122, t131, t132, t141, t142, t151, t152, t210, t211, t212, t222, t231, t232, t241, t242, t251, t252

t = [t110, t111, t112, t122, t131, t132, t141, t142, t151, t152, t210, t211, t212, t222, t231, t232, t241, t242, t251, t252]

@var m21000, m20100, m20010, m20001, m12000, m11100, m11010, m11001, m10200, m10110, m10101, m10020, m10011, m10002, m02100, m02010, m02001, m01200, m01110, m01101, m01020, m01011, m01002, m00210, m00201, m00120, m00111, m00102, m00021, m00012

m = [m21000, m20100, m20010, m20001, m12000, m11100, m11010, m11001, m10200, m10110, m10101, m10020, m10011, m10002, m02100, m02010, m02001, m01200, m01110, m01101, m01020, m01011, m01002, m00210, m00201, m00120, m00111, m00102, m00021, m00012]

par = [t112+t212, t112*t131+t212*t231, t112*t141+t212*t241, t112*t151+t212*t251,t111*t122+t211*t222, t111*t131+t211*t231, t111*t141+t211*t241, t111*t151+t211*t251,t111*t132+t211*t232, t111*t131*t141+t211*t231*t241, t111*t131*t151+t211*t231*t251,t111*t142+t211*t242, t111*t141*t151+t211*t241*t251, t111*t152+t211*t252,t110*t122*t131+t210*t222*t231, t110*t122*t141+t210*t222*t241, t110*t122*t151+t210*t222*t251,t110*t132+t210*t232, t110*t131*t141+t210*t231*t241, t110*t131*t151+t210*t231*t251,t110*t142+t210*t242, t110*t141*t151+t210*t241*t251, t110*t152+t210*t252,t110*t132*t141+t210*t232*t241, t110*t132*t151+t210*t232*t251, t110*t131*t142+t210*t231*t242,t110*t131*t141*t151+t210*t231*t241*t251, t110*t131*t152+t210*t231*t252,t110*t142*t151+t210*t242*t251, t110*t141*t152+t210*t241*t252]

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
