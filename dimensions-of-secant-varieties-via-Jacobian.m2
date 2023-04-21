-- This file computes the dimensions of the secant varieties sigma_r(M_{n,d}) and sigma_r(M_{n,lambda}) 
-- using the Jacobian matrix of the parametrization (first approach) and Terrachini lemma (second approach).


--in general:

r=2; n=5; d=2; --these parameters will be changed depending on moment variety
las=compositions(n, d); len=length(las); L={}; for i from 0 to len-1 do L=append(L,m_(las_i)); R=QQ[L]
R=QQ[L]
--define parameters:
S=QQ[toList(x_(1,1,1) ..< y_(r+1,n+1,d+1))]

--create parametrization for the secant variety:
M=apply(toList({0..length(L)-1}_0), i->0);
for k from 1 to r do
{Mx={}; for i from 0 to len-1 do (prod=1; for j from 0 to n-1 do 
    (if las_i_j!=0 then prod=prod*x_(k,j+1,las_i_j);); Mx=append(Mx,prod););
M=M+Mx;}

nums={}; nums=apply(length(generators(S)), i -> random(1,100));
toNums=map(QQ,S,nums);
A=toNums(jacobian(ideal(M)));
--rank of this Jacobiam matrix will be one more than the projective dimension of the variety:
rank A 

--if terminates, the below will give the implicit description of the moment variety:
f=map(S,R,M); 
I=ker f; 


------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
--another approach (Terracini lemma), may be faster:

r=2; n=5; d=4;
las=compositions(n, d); len=length(las); L={}; for i from 0 to len-1 do L=append(L,m_(las_i)); R=QQ[L]
R=QQ[L]
S=QQ[toList(x_(1,1) ..< y_(n+1,d+1))]
M={}; for i from 0 to len-1 do (prod=1; for j from 0 to n-1 do 
    (if las_i_j!=0 then prod=prod*x_(j+1,las_i_j);); M=append(M,prod);)

nums={}; nums=apply(length(generators(S)), i -> random(1,100));
toNums=map(QQ,S,nums); J=jacobian(ideal(M));
A=toNums(J);
rank A

for k from 1 to r-1 do
{
nums={}; nums=apply(length(generators(S)), i -> random(1,100));
toNums=map(QQ,S,nums);
A1=toNums(J); 
A=A||A1;
}

--rank of this matrix will be one more than the projective dimension of the variety:
A;
rank A 


------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
--computing dimension of the moment variety for a partition la:

r=2; n=5; d=2; --these parameters will vary
la={1,1,0,0,0}; --this partition will vary
las=toList(set(permutations(la))); len=length(las); L={}; for i from 0 to len-1 do L=append(L,m_(las_i)); 
R=QQ[L]; S=QQ[toList(x_(1,1,1) ..< y_(r+1,n+1,d+1))];

M=apply(toList({0..length(L)-1}_0), i->0);
for k from 1 to r do
{Mx={}; for i from 0 to len-1 do (prod=1; for j from 0 to n-1 do 
    (if las_i_j!=0 then prod=prod*x_(k,j+1,las_i_j);); Mx=append(Mx,prod););
M=M+Mx;}

nums={}; nums=apply(length(generators(S)), i -> random(1,100));
toNums=map(QQ,S,nums);
A=toNums(jacobian(ideal(M)));
rank A

f=map(S,R,M); I=ker f; 
(dim I-1, degree I)



nums={}; nums=apply(length(generators(S)), i -> random(1,100));
toNums=map(QQ,S,nums);
A=toNums(jacobian(ideal(M)))
rank A



