restart
needsPackage "SemidefiniteProgramming"
needsPackage "NumericalAlgebraicGeometry"   

solsSDP = method()
solsSDP (Matrix,RingElement) := (M',obj') -> (
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
    return I1+I2;
    for r from 1 to n-1 list  (
    	I3 := minors(r+1,M) + minors(n-r+1,X);
    	I := I1+I2+I3;
    	print (r, dim I, degree I);
        sys := first entries gens gb I;
        sols := solveSystem(sys);
        apply( sols, s -> apply( { obj, vars R, M, X }, k-> sub(k, matrix s)))
        )
    )        
-- The "NAG" package runs forever, since the system is not square ;(
-- Perhaps using a GB is better?
-- solveSystem(flatten entries gens gb I)

--QQ[a,b,c]
--M = matrix {{1,a,3-b},{a,5,c},{3-b,c,9+a}} ;
--obj = a+b+c ;
--I = solveSDP(M,obj)

--QQ[a,b]
--M = matrix {{1,a,b},{a,1,0},{b,0,1}}
--obj = a+b
--solveSDP(M,obj)

