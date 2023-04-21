#The code produces a generating function which allows us
#to extract all A-fibers  of degree 4, and for each such  (relevant)
#fiber we are solving the linear system and we are Grobner-reducing
#modulo the three known cubics. 


n := 4:
d := 4:
                                               interface(quiet=true):
with(Groebner):
faceideal := [
 m2020*m1111*m0202-m2110*m1021*m0202-m2020*m1201*m0112
 +m2200*m1021*m0112+m2110*m1201*m0022-m2200*m1111*m0022,
m2002*m1111*m0220-m2101*m1012*m0220-m2002*m1210*m0121
 +m2200*m1012*m0121+m2101*m1210*m0022-m2200*m1111*m0022,
m0202*m1111*m2020-m0202*m1120*m2011-m0211*m1102*m2020
 +m0211*m1120*m2002 +m0220*m1102*m2011-m0220*m1111*m2002]:
mvar := convert(indets(faceideal),list);
GB := Basis(faceideal,tdeg(mvar[])):
 
 
para := 
{m0013 = t110*t120*t131*t143+t210*t220*t231*t243, m0022 = t110*t120*t132*t142+t210*t220*t232*t242, m0031 = t110*t120*t133*t141+t210*t220*t233*t241, 
m0103 = t110*t121*t130*t143+t210*t221*t230*t243, m0112 = t110*t121*t131*t142+t210*t221*t231*t242, m0121 = t110*t121*t132*t141+t210*t221*t232*t241, 
m0130 = t110*t121*t133*t140+t210*t221*t233*t240, m0202 = t110*t122*t130*t142+t210*t222*t230*t242, m0211 = t110*t122*t131*t141+t210*t222*t231*t241, 
m0220 = t110*t122*t132*t140+t210*t222*t232*t240, m0301 = t110*t123*t130*t141+t210*t223*t230*t241, m0310 = t110*t123*t131*t140+t210*t223*t231*t240, 
m1003 = t111*t120*t130*t143+t211*t220*t230*t243, m1012 = t111*t120*t131*t142+t211*t220*t231*t242, m1021 = t111*t120*t132*t141+t211*t220*t232*t241, 
m1030 = t111*t120*t133*t140+t211*t220*t233*t240, m1102 = t111*t121*t130*t142+t211*t221*t230*t242, m1111 = t111*t121*t131*t141+t211*t221*t231*t241, 
m1120 = t111*t121*t132*t140+t211*t221*t232*t240, m1201 = t111*t122*t130*t141+t211*t222*t230*t241, m1210 = t111*t122*t131*t140+t211*t222*t231*t240, 
m1300 = t111*t123*t130*t140+t211*t223*t230*t240, m2002 = t112*t120*t130*t142+t212*t220*t230*t242, m2011 = t112*t120*t131*t141+t212*t220*t231*t241, 
m2020 = t112*t120*t132*t140+t212*t220*t232*t240, m2101 = t112*t121*t130*t141+t212*t221*t230*t241, m2110 = t112*t121*t131*t140+t212*t221*t231*t240, 
m2200 = t112*t122*t130*t140+t212*t222*t230*t240, m3001 = t113*t120*t130*t141+t213*t220*t230*t241, m3010 = t113*t120*t131*t140+t213*t220*t231*t240, 
m3100 = t113*t121*t130*t140+t213*t221*t230*t240}:


die := rand(-30..30): with(linalg):
P1 := 0: P2 := 0:  P3 := 0:

L := 0:  for j from 1 to n do L := L + cat(x,j): od:
Ld := expand(L^d): X := sort(convert(indets(L),list)):

partit := {}:

for i from 1 to nops(Ld) do
id := []:
for j from 1 to n do id := [id[],degree(op(i,Ld),X[j])]: od:
partit := partit union {sort(id)}:
od:


M1 := []: M2 := []:

for i from 1 to nops(Ld) do
id := []:
for j from 1 to n do id := [id[],degree(op(i,Ld),X[j])]: od:
if (not(member(d,id))) then
  M1 := [M1[],cat(m,id[])]:
 ff := 1:  for jj from 1 to n do  ff := ff*cat(t,jj,id[jj]):  od:
 M2 := [M2[],ff]:
P1 := P1 + ff:
P2 := P2 + ff*cat(m,id[]):
if (member(3,id)) then  P3 := P3 + ff*cat(m,id[]): fi:
fi:
od:

print(` `); var := sort(convert(indets(M2),list)):

lprint(M1); print(` `);
lprint(M2); print(` `); var;
tvar := convert(indets(para) minus indets(M1),list):


P2;
# P3;
Q := expand(P2^4):
nops(Q);
C := {coeffs(Q,var)}:

for cc in C do
if (not(type(cc,monomial)) and (nops(cc) > 10)) then
big := cc/content(cc):
zvar := []:
ansatz := 0:
for i from 1 to nops(big) do:
ansatz := ansatz + cat(z,i)*op(i,big)/lcoeff(op(i,big)):
zvar := [zvar[],cat(z,i)]: od:


{coeffs(collect(expand(subs(para,ansatz)),tvar,distributed),tvar)}:
winner := subs(solve(%),ansatz):
winner2 := NormalForm(winner,GB,tdeg(mvar[])):
# print(winner,winner2);
top := {coeffs(collect(expand(winner2),zvar,distributed),zvar)}:
if (top <> {0}) then
print(` `);
lprint(nops(big));
for p in top  do sfp := sort(factor(p)): lprint(sfp,`,`); od:
fi:
fi:
od:
time();