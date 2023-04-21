-- This Macaulay2 file was used to form a masked Hankel flattening matrix.
-- It related to Section 6: "Implicitization".



n=4
lambda={1,1,1,1} --does not need to be ordered, should be all positive integers
flattenPoint = 1 --a positive integer strictly less than the length of lambda
--r = 1 --a positive integer

kk=ZZ/101;--e.g. QQ or ZZ/101
rowLambda = lambda_{0..flattenPoint-1};
colLambda = lambda_{flattenPoint..#lambda-1};
constantList = (leng,a) -> (apply(leng, i-> a)) --for convenience
lexSubsets = (n,k) -> (
    select(toList(constantList(k,1)..constantList(k,n)), s -> (s==sort(s) and #set s == k))
    ) --built-in subsets command doesnt return lex ordering which is annoying
rowLabels = lexSubsets(n,#rowLambda);
rowLabels = apply(rowLabels, s -> flatten flatten apply(#rowLambda, i -> constantList(rowLambda#i,s#i))); 
colLabels = lexSubsets(n,#colLambda);
colLabels = apply(colLabels, s -> flatten flatten apply(#colLambda, i -> constantList(colLambda#i,s#i))); 
myIndices = unique flatten apply(rowLabels, r -> apply(colLabels, c -> (sort(r | c)))); 
myVars = apply(myIndices, i -> m_i);
R = kk[myVars];
missingVars = (flatten entries vars R)_(positions(myIndices, i-> (#set i < #lambda)));
visibleVars = (flatten entries vars R)_(positions(myIndices, i-> (#set i == #lambda)));
M = matrix apply(rowLabels, r -> apply(colLabels, c -> (m_(sort(r | c)))));
--I = minors(r+1,M);
--eliminate(missingVars,I)
-----S=kk[apply(select(myIndices, i-> (#set i == #lambda)), s -> a_s)];
-----ker(map(R/I,S,visibleVars),SubringLimit=>1)
MwithZeros = matrix apply(rowLabels, r -> apply(colLabels, c -> (if #unique(r | c) < #lambda then 0 else m_(sort(r | c)))));



