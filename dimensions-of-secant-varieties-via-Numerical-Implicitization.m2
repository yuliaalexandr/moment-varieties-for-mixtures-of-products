-- This file invokes the Macaulay2 package "NumericalImplicitization' to compute 
-- the dimensions of the secant varieties sigma_r(M_{n,d}) and sigma_r(M_{n,lambda}).
-- It was used for some of the computations reported in Section 5: "Secant Varieties".





debug needsPackage "NumericalImplicitization"

----

myNumericalDim = method(Options => options numericalSourceSample)
myNumericalDim (Matrix, Ideal, Point) := List => opts -> (F, I, p) -> (
    (F, I, p) = checkRings(F, I, {p});
    p0 := 1/norm(2, matrix p#0)*(matrix p#0);
    dF := sub(transpose jacobian F, p0);
    if I == 0 then (
		mapDim := numericalNullity(dF, true);
		return {#gens ring I - mapDim#0, mapDim#1#0}
	);
    sourceJacobian := sub(transpose jacobian I, p0);
    sourceDim := numericalNullity(sourceJacobian, true);
	mapDim = numericalNullity(sourceJacobian || dF, true);
    {sourceDim#0 - mapDim#0, sourceDim#1#0, mapDim#1#0}
)
myNumericalDim (Matrix, Ideal) := ZZ => opts -> (F, I) -> myNumericalDim(F, I, first numericalSourceSample(I, Software => opts.Software), opts)

---

howMany = (elt,lis) -> (--for convenience
    #positions(lis, l -> l == elt)
    )

----

ambientDimLambda = (n,lambda) -> (
    lift((product toList((n-#lambda+1)..n))/(product apply(unique lambda, i-> (howMany(i,lambda))!)),ZZ)-1
    )

----

expectedDimSecantLambda = (n,r,lambda) -> (--need n >= length lambda >= 3
    min((#unique lambda)*(n-1)*r + (r-1),ambientDimLambda(n,lambda)) --projective dim, this needs a min
    )


expectedFillingRankLambda = (n,lambda) -> (
    ceiling((ambientDimLambda(n,lambda) + 1)/((#unique lambda)*(n-1)+1))
    )

----

dimSecantLambda = (n,r,lambda) -> (
ell=#lambda;
k = max lambda;
R=CC[x_(1,1,1)..x_(n,k,r)]; --first index is for ambient dimension, second index is for moment order, third index is for component
myMatrices = apply(toList(1..r), r1 -> matrix apply(toList(1..n), n1 -> apply(toList(1..k), k1 -> x_(n1,k1,r1))));
allOnes = apply(toList(1..ell), k1 ->1);
allNs =  apply(toList(1..ell), k1 ->n);
myProducts = select(toList(allOnes..allNs), s -> (#set s == ell));
F=matrix {sum apply(myMatrices, M -> toList set apply(myProducts, s -> product(apply(toList(0..(ell-1)), ell1 -> M_(s#(ell1)-1,lambda#(ell1)-1)))))};--check the set! maybe safer to drop  
affineDim=myNumericalDim(F, ideal 0_R);
{affineDim#0-1, affineDim#1} --subtract 1 to go to projective dimension
)

----

ambientDimFull = (n,d) -> (
    binomial(n-1+d, d)-1 --subtract 1 to go to projective dimension
    )

----

dimSecantFull = (n,r,d) -> (
S = CC[y_1..y_n];
myVars = flatten entries vars S;
myMonomials = flatten entries basis(d,S);
myExponents = apply(myMonomials, m -> apply(myVars, v -> degree(v,m)));--this was all to get these exponent vectors
R = CC[x_(1,1,1)..x_(n,d,r)];
F = matrix {sum apply(toList(1..r), r1->apply(myExponents, i-> product apply(toList(1..n), n1-> (if i#(n1-1)==0 then 1 else x_(n1,i#(n1-1),r1)))))};
affineDim = myNumericalDim(F, ideal 0_R);
{affineDim#0-1, affineDim#1}
)

----

P^2 d=4

(dimSecantFull(3,1,4))#0 -- r=1 --> projective dim 11
(dimSecantFull(3,2,4))#0 -- r=2 --> projective dim 14
(dimSecantFull(3,3,4))#0 -- r=3 --> projective dim 
(dimSecantFull(3,4,5))#0 


P^3 d=4
(dimSecantFull(4,1,4))#0 --15
(dimSecantFull(4,2,4))#0 --27 --could be trouble?
(dimSecantFull(4,3,4))#0 --34
(dimSecantFull(4,4,4))#0


P^4 d=3
(dimSecantFull(5,1,3))#0 --14
(dimSecantFull(5,2,3))#0 --24
(dimSecantFull(5,3,3))#0 --34


P^4 d=4
(dimSecantFull(5,1,4))#0
(dimSecantFull(5,2,4))#0
(dimSecantFull(5,3,4))#0
(dimSecantFull(5,4,4))#0
(dimSecantFull(5,5,4))#0
(dimSecantFull(5,6,4))#0







myPartitions = {
{1,1,1},
{2,1,1},
{3,2,1},
{1,1,1,1},
{2,1,1,1},
{2,2,1,1},
{3,2,1,1},
{4,3,2,1},
{1,1,1,1,1},
{2,1,1,1,1},
{2,2,1,1,1},
{3,2,1,1,1},
{3,2,2,1,1},
{4,3,2,1,1},
{5,4,3,2,1}
}


maxNs = apply(myPartitions, lambda -> (last select(30, i -> (ambientDimLambda(i,lambda) <= 1000))))--200 too small for {5,4,3,2,1}


expectedFillingRanks = apply(#myPartitions, i -> expectedFillingRankLambda(maxNs#i, myPartitions#i))


elapsedTime results = flatten apply(#myPartitions, i -> (
	lambda = myPartitions#i;
	n = maxNs#i;
	rfill = expectedFillingRanks#i;
	dims1=dimSecantLambda(n,rfill-1,lambda);
	dims2=dimSecantLambda(n,rfill,lambda);
	expect1=expectedDimSecantLambda(n,rfill-1,lambda);
	expect2=expectedDimSecantLambda(n,rfill,lambda);
	print {{(n,rfill-1,lambda),expect1==dims1#0,expect1,dims1#0}, {(n,rfill,lambda),expect2==dims2#0,expect2,dims2#0}};
	{{(n,rfill-1,lambda),expect1==dims1#0,expect1,dims1#0,dims1#1}, {(n,rfill,lambda),expect2==dims2#0,expect2,dims2#0,dims2#1}}
	)
    );













elapsedTime joe = flatten apply(myPartitions, lambda -> (apply(toList((#lambda+1,1)..(10,10)), nr -> (--there are a bunch of counterexamples when n = #lambda for r = 1 and r=2 although sometimes it's still the expected dim
		n=nr#0;
		r=nr#1;
		myDim = dimSecantLambda(n,r,lambda);
		expectedDim = expectedDimSecantLambda(n,r,lambda);
		if myDim#0 == expectedDim then (n,r,lambda,true,expectedDim) else (n,r,lambda,myDim#0,false,expectedDim,myDim#0)
		))))








