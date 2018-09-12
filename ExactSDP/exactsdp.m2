restart
needsPackage "NumericalAlgebraicGeometry"   
solveSDP = method()
solveSDP (Matrix,RingElement) := (M',obj') -> (
    RM := ring M';
    x := symbol x;
    n := numcols M';
    assert(numrows M' == n);
    R := QQ (monoid [x_1..x_(binomial(n+1,2)),gens RM]);
    M := sub(M',R);
    obj := sub(obj',R);
    X := genericSymmetricMatrix(R,R_0,n);
    MX := M*X;
    aux := trace MX - obj;
    -- compute linear equations for primal
    I1 := ideal apply(gens RM,k->diff(sub(k,R),aux));
    -- complementary slackness
    I2 := ideal flatten entries(MX);
    -- iterate through ranks
    -- replace by Pataki rank, or check if zerodim
    for r from 2 to n-1 list  (
    	I3 := minors(r+1,M) + minors(n-r+1,X);
    	I := I1+I2+I3;
    	print (r, dim I, degree I);
	sols := solveSystem(first entries gens gb I);
	-- sols := elements(S);
	apply( sols, s -> apply( { obj, vars ring M, M, X }, k-> sub(k, matrix s)))
	)
    )        
-- The "NAG" package runs forever, since the system is not square ;(
-- Perhaps using a GB is better?
-- solveSystem(flatten entries gens gb I)
-- elements(S)

end


restart
load "exactsdp.m2"
QQ[a,b,c]
M = matrix {{1,a,3-b},{a,5,c},{3-b,c,9+a}} ;
obj = a+b+c ;
solveSDP(M,obj)

-- solve twice, just to check
QQ[a,b]
M = matrix {{1,a,b},{a,1,0},{b,0,1}}
obj = a+b
solveSDP(M,obj)



solveSDP(matrix {{1,a,b},{a,1,0},{b,0,1}},a+b)







 