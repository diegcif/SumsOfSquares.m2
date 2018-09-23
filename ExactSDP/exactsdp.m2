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

sdpIdeal = method( Options => {CoefficientRing => null, Square=>false} )
sdpIdeal (Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    kk := o.CoefficientRing;
    if kk===null then kk = ring C;
    local x; x = symbol x;
    local y; y = symbol y;
    n := numRows C;
    n2 := binomial(n+1,2);
    m := numRows b;
    R := kk(monoid [x_0..x_(n2-1),y_0..y_(m-1)]);
    (C,A,b) = toRing_R (C,A,b);
    -- primal/dual matrices
    X := genericSymmetricMatrix(R,R_0,n);
    y = drop(gens R, n2);
    Z := C - sum(for i to m-1 list y_i * A_i);
    -- primal feasibility
    I1 := ideal(for i to m-1 list trace(A_i*X)-b_(i,0));
    -- complementary slackness
    I2 := if o.Square then ideal smat2vec entries(Z*X + X*Z)
        else ideal flatten entries(X*Z);
    return (I1+I2, X, Z);
    )
sdpIdeal (Matrix, Sequence, Matrix, ZZ) := o -> (C,A,b,r) -> (
    (I,X,Z) = sdpIdeal(C,A,b);
    Ir := minors(r+1,Z) + minors(n-r+1,X);
    return (I+Ir, X, Z);
    )

toRing = (R,C,A,b) -> (
    C = promote(C,R);
    A = apply(A, Ai -> promote(Ai,R));
    b = promote(b,R);
    return (C,A,b);
    )

refineSDP = (C,A,b,X0,y0) -> (
    (J,X,Z) := sdpIdeal(C,A,b,Square=>true);
    pt := smat2vec entries X0 | flatten entries y0;
    pt' := first refine (polySystem J, {pt});
    if pt'#SolutionStatus==RefinementFailure then(
        print "refinement failed";
        return (X0,y0) );
    L := coordinates pt';
    m := numColumns y0;
    X1 := matrix vec2smat drop(L,-m);
    y1 := matrix transpose {take(L,-m)};
    return (X1,y1);
    )

C = matrix(QQ,{{1,0,3},{0,5,0},{3,0,9}})
A1 = -matrix{{0,1,0},{1,0,0},{0,0,1}}
A2 = -matrix{{0,0,-1},{0,0,0},{-1,0,0}}
A3 = -matrix{{0,0,0},{0,0,1},{0,1,0}}
b = -matrix{{1},{1},{1}}
A = (A1,A2,A3)
(J,X,Z) = sdpIdeal(C,A,b)

-- refine SDP solution
C = matrix {{1,0},{0,2}};
A1 = matrix {{0,1},{1,0}};
A = sequence(A1)
b = matrix {{-1}};
(X0,y0,Q0) = solveSDP(C,A,b);
(X1,y1) = refineSDP(C,A,b,X0,y0)
