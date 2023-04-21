-- This file invokes the Macaulay2 packages "FourTiTwo" and "Polyhedra" to compute 
-- the ideals of the toric ideals for M_{n,d} and M_{n,lambda}.
-- We also compute the betti table of the ideals and degree of the projective varieties here.
-- It was used for some of the computations reported in Section 4: "Toric Combinatorics".
-- Various results were copied and pasted here into the commments below.

needsPackage "FourTiTwo"
needsPackage "Polyhedra"

----------------------------------------------------------------------
----------------------------------------------------------------------

constantList = (leng,elt) -> (--for convenience
    apply(leng, i-> elt)
    ) 


howMany = (elt,lis) -> (--for convenience
    #positions(lis, l -> l == elt)
    )


allLists = (n,k) -> (--for convenience...check this is RIGHT?!
    myLexSubsets = select(toList(constantList(k,1)..constantList(k,n)), s -> (s==sort(s) and #set s == k));
    flatten apply(myLexSubsets, s -> permutations s)
    )

where = (elt,lis) -> (
    apply(positions(lis, i -> (i==elt)), j->j+1)
    )

isIn = (elt, lis) -> (
    unique sort (lis | {elt}) == unique sort lis
    )

makeRowLabels = (n,lambda) -> (
    flatten apply(toList(1..n), i -> apply(unique lambda, p -> {i,p}))
    )


makeColLabels = (n,lambda) -> (
    myColLabels = allLists(n,#lambda); --I shouldnt have deleted this
    unique apply(myColLabels, s -> sort flatten flatten apply(#lambda, i -> constantList(lambda#i,s#i)))
    )
    

makeRowLabelsFull = (n,d) -> (
    toList({1,1}..{n,d})
    )
     
        
makeColLabelsFull = (n,d) -> (
    unique apply(toList(constantList(d,1)..constantList(d,n)), s -> sort s) 
    )

    
makeMatrix = (rowLabels,colLabels) -> (--rowLabels=makeRowLabels(n,lambda); colLabels=makeColLabels(n,lambda); or can use Full variants
     matrix apply(rowLabels, r -> apply(colLabels, c -> (if (howMany(r#0,c))==(r#1) then 1 else 0)))
    )


makeRing = (colLabels,kk) -> (--colLabels=makeColLabels(n,lambda); or can use Full variants
    kk[apply(colLabels, c -> m_c)]
    )


makeIdeal = (A,R) -> (--A=makeMatrix(rowLabels,colLabels); R=makeRing(colLabels,kk);
    toBinomial(toricMarkov(A),R)
    )

----------------------------------------------------------------------
----------------------------------------------------------------------


n=3
d=6
kk=QQ

rowLabelsFull = makeRowLabelsFull(n,d);
colLabelsFull = makeColLabelsFull(n,d); 
AFull = makeMatrix(rowLabelsFull,colLabelsFull);
RFull = makeRing(colLabelsFull,kk);
IFull = makeIdeal(AFull,RFull);
(dim IFull, codim IFull)
degree IFull
JFull = mingens IFull;
I2Full = ideal JFull;
betti I2Full
betti res I2Full

PFull = convexHull AFull;
fVector PFull









-----------------
-----------------
n=3
d=1,2,3,4,5,6,7 (Full case)

--d=1: all of P^2, fVector {3,3,1}
--d=2: all of P^5, fVector {6, 15, 20, 15, 6, 1}
--d=3: in P^9, dim 8, degree 3, fVector {10, 45, 118, 196, 211, 145, 60, 13, 1}, betti table as follows
             0 1
o68 = total: 1 1
          0: 1 .
          1: . .
          2: . 1
	  
--d=4: in P^{14}, dim 11, degree 19, fVector {15, 105, 447, 1272, 2518, 3513, 3435, 2298, 1001, 258, 32, 1}, betti table as follows 
             0 1 2 3
o82 = total: 1 4 6 3
          0: 1 . . .
          1: . . . .
          2: . 4 . .
          3: . . 3 .
          4: . . 3 3


--d=5: in P^{20}, dim 14, degree 139, didnt attempt fVector, betti table as follows 
             0  1   2   3   4  5  6
o96 = total: 1 23 106 195 180 87 18
          0: 1  .   .   .   .  .  .
          1: .  .   .   .   .  .  .
          2: . 11   .   .   .  .  .
          3: . 12  69  45   3  .  .
          4: .  .  37 138 123 27  .
          5: .  .   .  12  54 60 18



--d=6: in P^{27}, dim 17, degree 1428, didnt attempt fVector, betti table as follows 	  
              0  1   2    3    4     5     6    7    8    9  10
o108 = total: 1 84 779 3099 7155 11234 12690 9999 5001 1364 154
           0: 1  .   .    .    .     .     .    .    .    .   .
           1: .  .   .    .    .     .     .    .    .    .   .
           2: . 24   .    .    .     .     .    .    .    .   .
           3: . 54 346  320   84    16     1    .    .    .   .
           4: .  6 412 2100 3070  1648   438   60    3    .   .
           5: .  .  12  535 2999  5688  4038 1178  147    3   .
           6: .  .   9  144 1002  3882  8213 8751 4809 1310 138
           7: .  .   .    .    .     .     .   10   42   51  16


--d=7: in P^{35}, dim 20, degree 14922, mingen degrees as follows
              0   1
o119 = total: 1 367
           0: 1   .
           1: .   .
           2: .  46
           3: . 168
           4: . 135
           5: .  18



-----------------
-----------------
n=4
d=1,2,3,4,5,6,7 (Full case)

--d=1: all of P^3, fVector {4,6,4,1}

--d=2: in P^9, dim 7, degree 4, fVector {10, 42, 96, 130, 106, 50, 12, 1}, betti table as follows 
              0 1 2
o146 = total: 1 2 1
           0: 1 . .
           1: . 2 .
           2: . . 1

--d=3: in P^{19}, dim 11, degree 47, fVector {20, 166, 754, 2108, 3858, 4768, 4017, 2284, 848, 192, 23, 1}, betti table as follows
              0  1   2   3   4   5   6   7  8
o160 = total: 1 18 122 302 406 403 284 106 16
           0: 1  .   .   .   .   .   .   .  .
           1: . 14  12   .   .   .   .   .  .
           2: .  4 110 280 246  56   4   .  .
           3: .  .   .  22 160 347 280 106 16	   
	   
--d=4: in P^{34}, dim 15, degree 1072, didnt attempt fVector, mingen degrees as follows 
              0  1
o173 = total: 1 80
           0: 1  .
           1: . 52
           2: . 28




















n=5
lambda={2,1}
kk=QQ

rowLabels = makeRowLabels(n,lambda)
colLabels = makeColLabels(n,lambda)--currently this fails
A = makeMatrix(rowLabels,colLabels)

R = makeRing(colLabels, kk)
I = makeIdeal(A, R);
(dim I, codim I)
degree I
J = mingens I;
I2 = ideal J;
betti I2
betti res I2

--Here are computations for the singular locus.  They seem to be costly
--temp = jacobian(I2);
--sing = ideal(temp) + I2; 
--sing2 = radical(sing);
--codim sing2



P = convexHull A;
fVector P


--This is an alternative way to get the degree.  Generally it seems slower
--(volume P)*((dim P)!)











----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {1}
n=1,2,3,4,5,6,7,8,9,10

--these are all of P^{n-1} respectively, fVector is {binom(n,1),binom(n,2),...,binom(n,n)}

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {1,1}
n=2,3,4,5,6,7,8,9,10


--n=2: all of P^0, fVector {1}


--n=3: all of P^2, fVector {3,3,1}


--n=4: in P^5, dim 3, degree 4, fVector {6,12,8,1}, betti table as follows
--       0 1 2
--total: 1 2 1
--    0: 1 . .
--    1: . 2 .
--    2: . . 1              


--n=5: in P^9, dim 4, degree 11, fVector {10,30,30,10,1}, betti table as follows
--       0  1  2  3  4 5
--total: 1 10 20 26 20 5
--    0: 1  .  .  .  . .
--    1: . 10 15  .  . .
--    2: .  .  5 26 20 5
	   
	   
--n=6: in P^14, dim 5, degree 26, fVector {15,60,80,45,12,1}, betti table as follows
--       0  1   2   3   4    5   6   7  8 9
--total: 1 30 121 330 756 1044 819 376 90 7
--    0: 1  .   .   .   .    .   .   .  . .
--    1: . 30 106 114   .    .   .   .  . .
--    2: .  .  15 216 756 1044 819 376 90 6
--    3: .  .   .   .   .    .   .   .  . 1  
	   
	   
--n=7: in P^20, dim 6, degree 57, fVector {21,105,175,140,63,14,1}, mingen degrees as follows
--       0  1
--total: 1 70
--    0: 1  .
--    1: . 70


--n=8: in P^{27}, dim 7, degree 120, fVector {28,168,336,350,224,84,16,1}, mingen degrees as follows
--       0   1
--total: 1 140
--    0: 1   .
--    1: . 140


--n=9: in P^{35}, dim 8, degree 247, fVector {36,252,588,756,630,336,108,18,1}, mingen degrees as follows 
--       0   1
--total: 1 252
--    0: 1   .
--    1: . 252


--n=10: in P^{44}, dim 9, degree 502, fVector {45,360,960,1470,1512,1050,480,135,20,1}, mingen degrees as follows
--       0   1
--total: 1 420
--    0: 1   .
--    1: . 420


--CONCLUSIONS/CONJECTURES for lambda = {1,1}: 
--1) in P^{binom(n,2)-1}
--2) dim n-1 if n>2
--3) degree from 2nd column of Euler's triangle (A000295)
--4) fVector given by (A271238O)
--5) minimally generated by 2*binom(n,4) quadrics
--6) Cohen-Macaulay




----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {2,1}
n=2,3,4,5,6,7,8,9,10


--n=2: all of P^1, fVector {2,1}


--n=3: in P^5, dim 4, degree 3, fVector {6,15,18,9,1}, betti table as follows
--       0 1
--total: 1 1
--    0: 1 .
--    1: . .
--    2: . 1


--n=4: in P^{11}, dim 6, degree 16, fVector {12,54,110,108,52,12,1}, betti table as follows
--       0  1  2  3  4 5
--total: 1 10 40 56 30 5
--    0: 1  .  .  .  . .
--    1: .  6  .  .  . .
--    2: .  4 40 56 30 4
--    3: .  .  .  .  . 1



--n=5: in P^{19}, dim 8, degree 65, fVector {20,130,370,550,480,260,85,15,1}, betti table as follows
--       0  1   2    3    4    5    6    7    8   9  10 11
--total: 1 40 390 1903 4675 6770 6245 3898 1881 720 150 11
--    0: 1  .   .    .    .    .    .    .    .   .   .  .
--    1: . 30  70    .    .    .    .    .    .   .   .  .
--    2: . 10 320 1903 4675 6765 6110 3325 1101 330  60  5
--    3: .  .   .    .    .    5  135  573  780 390  90  5
--    4: .  .   .    .    .    .    .    .    .   .   .  1



--n=6: in P^{29}, dim 10, degree 246, fVector {30,255,930,1860,2362,2040,1215,490,126,18,1}, mingen degrees as follows
--        0   1
-- total: 1 110
--     0: 1   .
--     1: .  90
--     2: .  20



--n=7: in P^{41}, dim 12, degree 917, fVector {42,441,1960,4935,8239,9800,8533,5467,2541,826,175,21,1}, mingens degrees as follows
--       0   1
--total: 1 245
--    0: 1   .
--    1: . 210
--    2: .  35


--n=8: in P^{55}, dim 14, degree 3424, didn't attempt fVector or betti res, mingen degrees as follows 
--       0   1
--total: 1 476
--    0: 1   .
--    1: . 420
--    2: .  56


--n=9: in P^{71}, dim 16, degree 12861, mingen degrees as follows
--       0   1
--total: 1 840
--    0: 1   .
--    1: . 756
--    2: .  84


--n=10: in P^{89}, dim 18, degree 48610, mingen degrees as follows
--       0    1
--total: 1 1380
--    0: 1    .
--    1: . 1260
--    2: .  120


--CONCLUSIONS/CONJECTURES for lambda = {2,1}: 
--1) in P^{n*(n-1)-1}
--2) dim 2*(n-1) if n>2
--3) degrees are 1,3,16,65,246,917,3424,12861,48610 for n=2,3,4,5,6,7,8,9,10 (not in OEIS)
--4) fVector: number of edges is  (4*(n-1)^3 + (n-1)^2 - 3*(n-1))/2 (A172073); number of facets is 3*n (n > 2)
--5) minimally generated by 6*binom(n+3,4) quadrics (A033487) and binom(n,3) cubics
--6) Cohen-Macaulay 


----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {1,1,1}
n=3,4,5,6,7,8,9,10

--n=3: all of P^0, fVector {1}



--n=4: all of P^3, fVector {4,6,4,1}



--n=5: in P^9, dim 4, degree 11, fVector {10,30,30,10,1}, betti table as follows
--       0  1  2  3  4 5
--total: 1 10 20 26 20 5
--    0: 1  .  .  .  . .
--    1: . 10 15  .  . .
--    2: .  .  5 26 20 5



--n=6: in P^{19}, dim 5, degree 66, fVector {20,90,120,60,12,1}, mingen degrees as follows
--       0  1
--total: 1 69
--    0: 1  .
--    1: . 69


--n=7: in P^{34}, dim 6, degree 302, fVector {35,210,350,245,84,14,1}, mingen degrees as follows
--       0   1
--total: 1 273
--    0: 1   .
--    1: . 273


--n=8: in P^{55}, dim 7, degree 1191, fVector {56,420,840,770,392,112,16,1}, mingen degrees as follows
--       0   1
--total: 1 812
--    0: 1   .
--    1: . 812


--n=9: in P^{83}, dim 8, degree 4293, fVector {84,756,1764,2016,1386,588,144,18,1}, mingen degrees as follows
--       0    1
--total: 1 2016
--    0: 1    .
--    1: . 2016


--n=10: was taking a bit long so I stopped it


--CONCLUSIONS/CONJECTURES for lambda = {1,1,1}: 
--1) in P^{binom(n,3)-1}
--2) dim n-1 if n>3
--3) degrees are third column of Euler's triangle (A008292)
--4) minimally generated by quadrics
--5) number of edges in fVector here (lambda = {2,1,1}) matches number of quadrics in {2,1} ideal w/ same n??
--- also number of vertices in fVector here matches number of cubics in {2,1} ideal w/ same n
--6) Cohen-Macaulay I guess (...is the toric variety smooth?)



----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {2,1,1}
n=3,4,5,6,7,8,9,10


--n=3: all of P^2, fVector {3,3,1}


--n=4: in P^{11}, dim 6, degree 16, fVector {12,54,110,108,52,12,1}, betti table as follows 
--       0  1  2  3  4 5
--total: 1 10 40 56 30 5
--    0: 1  .  .  .  . .
--    1: .  6  .  .  . .
--    2: .  4 40 56 30 4
--    3: .  .  .  .  . 1


--n=5: in P^{29}, dim 8, degree 355, fVector {30,240,720,1035,810,370,100,15,1}, mingen degrees as follows 
--       0   1
--total: 1 130
--    0: 1   .
--    1: . 110
--    2: .  20


--n=6: in P^{59}, dim 10, degree 4326, fVector {60,690,2690,5175,5876,4320,2130,700,147,18,1}, mingen degrees as follows
--       0   1
--total: 1 705
--    0: 1   .
--    1: . 645
--    2: .  60



--n=7: too slow


--TOADD CONJECTURES/CONJECTURES for lambda = {2,2,1}

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {3,2,1}
n=3,4,5,6,7,8,9,10

--n=3: in P^5, dim 4, degree 3, fVector {6,15,18,9,1}, betti table as follows 
--       0 1
--total: 1 1
--    0: 1 .
--    1: . .
--    2: . 1



--n=4: in P^{23}, dim 9, degree 352, fVector {24,240,978,1968,2176,1392,528,120,16,1}, mingen degrees as follows
--       0   1
--total: 1 178
--    0: 1   .
--    1: .  18
--    2: . 160


--n=5: in P^{59}, dim 12, degree 19145, fVector took more than 5 mins, mingen degrees as follows
--       0    1
--total: 1 1360
--    0: 1    .
--    1: .  360
--    2: . 1000


--TOADD CONCLUSION/CONJECTURES for lambda = {3,2,1}


----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {1,1,1,1}
n=4,5,6,7,8,9


--n=4: all of P^0, fVector {1}


--n=5: all of P^4, fVector {5,10,10,5,1}


--n=6: in P^{14}, dim 5, degree 26, fVector {15,60,80,45,12,1}, betti table as follows 
--       0  1   2   3   4    5   6   7  8 9
--total: 1 30 121 330 756 1044 819 376 90 7
--    0: 1  .   .   .   .    .   .   .  . .
--    1: . 30 106 114   .    .   .   .  . .
--    2: .  .  15 216 756 1044 819 376 90 6
--    3: .  .   .   .   .    .   .   .  . 1



--n=7: in P^{34}, dim 6, degree 302, fVector {35,210,350,245,84,14,1}, mingen degrees as follows (betti res had "too many heap sections")
--       0   1
--total: 1 273
--    0: 1   .
--    1: . 273




--n=8: in P^{69}, dim 7, degree 2416, fVector {70,560,1120,980,448,112,16,1}, mingen degrees as follows
--       0    1
--total: 1 1378
--    0: 1    .
--    1: . 1378



--n=9: in P^{126}, dim 8, degree 15619, fVector {126,1260,2940,3150,1890,672,144,18,1}, mingen degrees as follows (this took awhile)
--       0    1
--total: 1 5094
--    0: 1    .
--    1: . 5094


----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {2,1,1,1}
n=4,5,6

--n=4: all of P^3, fVector {4,6,4,1}


--n=5: in P^{19}, dim 8, degree 65, fVector {20,130,370,550,480,260,85,15,1}, betti table as follows
--       0  1   2    3    4    5    6    7    8   9  10 11
--total: 1 40 390 1903 4675 6770 6245 3898 1881 720 150 11
--    0: 1  .   .    .    .    .    .    .    .   .   .  .
--    1: . 30  70    .    .    .    .    .    .   .   .  .
--    2: . 10 320 1903 4675 6765 6110 3325 1101 330  60  5
--    3: .  .   .    .    .    5  135  573  780 390  90  5
--    4: .  .   .    .    .    .    .    .    .   .   .  1



--n=6: in P^{59}, dim 10, degree 4326, fVector {60,690,2690,5175,5876,4320,2130,700,147,18,1}, mingen degrees as follows
--       0   1
--total: 1 705
--    0: 1   .
--    1: . 645
--    2: .  60




--n=7: taking too long



----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {2,2,1,1}
n=4,5,6

--n=4: in P^5, dim 3, degree 4, fVector {6,12,8,1}, betti table as follows
--       0 1 2
--total: 1 2 1
--    0: 1 . .
--    1: . 2 .
--    2: . . 1


--n=5: in P^{29}, dim 8, degree 355, fVector {30,240,720,1035,810,370,100,15,1}, mingen degrees as follows
--       0   1
--total: 1 130
--    0: 1   .
--    1: . 110
--    2: .  20


--n=6: in P^{89}, dim 10, degree 29646, fVector {90,1260,5220,9945,10530,6840,2880,810,153,18,1}, mingen degrees as follows
--       0    1
--total: 1 1830
--    0: 1    .
--    1: . 1710
--    2: .  120


----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {3,2,1,1}
n=4,5

--n=4: in P^{11}, dim 6, degree 16, fVector {12,54,110,108,52,12,1}, betti table as follows
--       0  1  2  3  4 5
--total: 1 10 40 56 30 5
--    0: 1  .  .  .  . .
--    1: .  6  .  .  . .
--    2: .  4 40 56 30 4
--    3: .  .  .  .  . 1



--n=5: in P^{59}, dim 12, degree 19145, didnt do fVector, mingen degrees as follows
--       0    1
--total: 1 1360
--    0: 1    .
--    1: .  360
--    2: . 1000

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
lambda = {4,3,2,1}
n=4,5

--n=4: in P^{23}, dim 9, degree 352, fVector {24,240,978,1968,2176,1392,528,120,16,1}, mingen degrees as follows
--       0   1
--total: 1 178
--    0: 1   .
--    1: .  18
--    2: . 160


--n=5: in P^{119}, taking too long




















-- This part of the file was used to develop the proof of Theorem 21.
-- In effect here we apply the map phi to the ideal J below.




lambda = {3,2,1}
n=5
rowLabels = makeRowLabels(n,lambda)
colLabels = makeColLabels(n,lambda)--currently this fails
A = makeMatrix(rowLabels,colLabels)

R = makeRing(colLabels, kk)
I = makeIdeal(A, R);
(dim I, codim I)
degree I
J = mingens I;
I2 = ideal J;
betti I2
LinEqns = ideal apply(colLabels, s->(
	mult3 = select(toList(1..n), i -> howMany(i,s)==3); 
	mult2 = select(toList(1..n), i -> howMany(i,s)==2); 
	mult1 = select(toList(1..n), i -> howMany(i,s)==1);
	t = sort flatten {mult1, mult3, mult3, mult2, mult2, mult2};
	m_s - m_t));
I3 = LinEqns + I2; 
	








kk=QQ
n=5
lambda1 = {2,1,1}

lambda1 = sort lambda1;
rowLabels1 = makeRowLabels(n,lambda1);
colLabels1 = makeColLabels(n,lambda1);
A1 = makeMatrix(rowLabels1,colLabels1);
R1 = makeRing(colLabels1, kk);
I1 = makeIdeal(A1, R1);
J1 = ideal mingens I1;

lambda2 = toList(1..(#lambda1))--antisymmetrized
rowLabels2 = makeRowLabels(n,lambda2);
colLabels2 = makeColLabels(n,lambda2);
A2 = makeMatrix(rowLabels2,colLabels2);
R2 = makeRing(colLabels2, kk);
I2 = makeIdeal(A2, R2);
J2 = ideal mingens I2;
symSubs = apply(colLabels2, s->(
	mu = {};
	for i in unique lambda1 do(
	    posi = where(i, lambda1);
	    mu = mu | flatten apply(unique select(s, j ->isIn(howMany(j,s),posi)), k -> constantList(i,k));
	    );
	mu = sort mu;
	use R1;
	m_(mu)));
phi = map(R1,R2,symSubs);
J3 = phi(J2);
J3 == J1



--a specific example of the linear equations
--LinEqns = ideal apply(colLabels, s->(
--	mult3 = select(toList(1..n), i -> howMany(i,s)==3); 
--	mult2 = select(toList(1..n), i -> howMany(i,s)==2); 
--	mult1 = select(toList(1..n), i -> howMany(i,s)==1);
--	t = sort flatten {mult1, mult3, mult3, mult2, mult2, mult2};
--	m_s - m_t));





checkSymStatement = (n,lambda1) -> (
kk=QQ;
lambda1 = sort lambda1;
rowLabels1 = makeRowLabels(n,lambda1);
colLabels1 = makeColLabels(n,lambda1);
A1 = makeMatrix(rowLabels1,colLabels1);
R1 = makeRing(colLabels1, kk);
I1 = makeIdeal(A1, R1);
J1 = ideal mingens I1;
lambda2 = toList(1..(#lambda1));--antisymmetrized
rowLabels2 = makeRowLabels(n,lambda2);
colLabels2 = makeColLabels(n,lambda2);
A2 = makeMatrix(rowLabels2,colLabels2);
R2 = makeRing(colLabels2, kk);
I2 = makeIdeal(A2, R2);
J2 = ideal mingens I2;
symSubs = apply(colLabels2, s->(
	mu = {};
	for i in unique lambda1 do(
	    posi = where(i, lambda1);
	    mu = mu | flatten apply(unique select(s, j ->isIn(howMany(j,s),posi)), k -> constantList(i,k));
	    );
	mu = sort mu;
	use R1;
	m_(mu)));
phi = map(R1,R2,symSubs);
J3 = phi(J2);
J3 == J1
    )

checkSymStatement(2,{1,1})
checkSymStatement(3,{1,1})
checkSymStatement(4,{1,1})
checkSymStatement(5,{1,1})
checkSymStatement(6,{1,1})
checkSymStatement(3,{1,1,1})
checkSymStatement(4,{1,1,1})
checkSymStatement(5,{1,1,1})
checkSymStatement(3,{2,1,1})
checkSymStatement(4,{2,1,1})
checkSymStatement(5,{2,1,1})
checkSymStatement(4,{1,1,1,1})
checkSymStatement(4,{2,1,1,1})
checkSymStatement(4,{2,2,1,1})
checkSymStatement(4,{3,2,1,1})
checkSymStatement(5,{3,2,2,1,1})--this and above all true...this one took awhile


