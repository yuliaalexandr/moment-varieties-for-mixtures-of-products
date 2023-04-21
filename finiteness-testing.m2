--This files tests if for a fixed partition increasing n by 1 produces any new generators


--la=(2,1,1)
--start with n=6
n=6; d=4; la={2,1,1,0,0,0};  
las=toList(set(permutations(la))); len=length(las); L={}; for i from 0 to len-1 do L=append(L,m_(las_i)); R=QQ[L]
S=QQ[toList(x_(1,1) ..< y_(n+1,la_0+1))]
T= QQ[join(generators(S),generators(R)), MonomialOrder=>Eliminate length(generators(S))]
M={}; for i from 0 to len-1 do (prod=1; for j from 0 to n-1 do 
    (if las_i_j!=0 then prod=prod*x_(j+1,las_i_j);); M=append(M,prod);)
L=apply(toList({length(generators S)..length(generators T)-1}#0), i-> (generators T)#i)
I=ideal(L-M); I=ideal(selectInSubring(1,gens gb(I)));
toString mingens I
betti mingens I
codim I, degree I


---let's go from to n=7:
Gs={}; for j from 0 to n do
(gs=apply(generators(R),  i->m_(insert(j,0,(baseName i)#1)));
Gs=join(Gs,gs);)
bigR=QQ[toList(set(Gs))];


JJ={};
for j from 0 to n do
(
zeros={}; for i from 1 to length(generators(S)) do zeros=append(zeros,0);
ts=apply(generators(R), i->m_(insert(j,0,(baseName i)#1)));
addZeros=map(bigR,T,join(zeros,ts));
G=first entries mingens addZeros(I);
JJ=join(JJ,G);
)

length(toList(set(JJ)));
II=ideal(toList(set(JJ)));
betti mingens II

II:ideal(m_{2,1,1,0,0,0,0}) == II --true! so indeed n=7 does not produce any new generators!


---------------------------------------------------------------------------------------------

