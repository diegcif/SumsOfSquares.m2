newPackage(
    "SOS",
    Version => "1.5", 
    Date => "November 17, 2017",
    Authors => {
     {Name => "Helfried Peyrl", 
      Email => "peyrl@control.ee.ethz.ch",
      HomePage => "https://scholar.google.com/citations?user=cFOV7nYAAAAJ&hl=de"},
     {Name => "Pablo A. Parrilo",
      Email => "parrilo@mit.edu",
      HomePage => "http://www.mit.edu/~parrilo/"},
     {Name => "Diego Cifuentes",
      Email => "diegcif@mit.edu",
      HomePage => "http://www.mit.edu/~diegcif/"}
    },
    Headline => "Sum-of-Squares Package",
    DebuggingMode => true,
    Configuration => {"CSDPexec"=>"csdp","SDPAexec"=>"sdpa","DefaultSolver"=>null},
    AuxiliaryFiles => true,
    -*
    The following two settings make use of a cached version of example output.
    When the documentation changes, the developer has to exchange true and
    false below, create the new example files, commit them, and change it back.
    *-
    UseCachedExampleOutput => false,
    CacheExampleOutput => true,
    PackageImports => {"SimpleDoc","FourierMotzkin"},
    PackageExports => {}
)

export {
--Types
    "SOSPoly",
    "SDPResult",
--Methods/Functions
    "sosPoly",
    "solveSOS",
    "sosdecTernary",
    "sumSOS",
    "LDLdecomposition",
    "solveSDP",
    "sosInIdeal",
    "lowerBound",
    "roundPSDmatrix",
    "checkSolver",
    "smat2vec",
    "vec2smat",
    "recoverSolution",
--only for debugging
--Method options
    "GramMatrix",
    "MomentMatrix",
    "RndTol",
    "UntilObjNegative",
    "Solver",
    "TraceObj",
    "Scaling"
}

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

makeGlobalPath = (fname) -> (
    tmp := temporaryFileName();
    r := run( "which '" | fname | "' > " | tmp);
    if r>0 then(
        print("Warning: " | fname | " executable was not found.");
        return;
        );
    fname = replace("\n","",get tmp);
    if first fname != "/" then fname = currentDirectory() | fname;
    return "'" | fname | "'";
    )

csdpexec = makeGlobalPath ((options SOS).Configuration)#"CSDPexec"
sdpaexec = ((options SOS).Configuration)#"SDPAexec"
defaultSolver = ((options SOS).Configuration)#"DefaultSolver"

--##########################################################################--
-- TYPES
--##########################################################################--

SOSPoly = new Type of HashTable

-- constructor for sos decomposition
sosPoly = method()
sosPoly (Ring, List, List) := SOS => (R, polys, coeffs) -> (
    new SOSPoly from {
        ring => R,
        gens => polys,
        coefficients => coeffs
        }
    )
sosPoly (List, List) := SOS => (polys, coeffs) -> sosPoly(ring polys#0,polys,coeffs)

ring SOSPoly := S -> S#ring

gens SOSPoly := o -> S -> S#gens

coefficients SOSPoly := o -> S -> S#coefficients

length SOSPoly := S -> #(S#gens)

substitute (SOSPoly,Ring) := (S,R) ->
    sosPoly(for g in S#gens list sub(g,R), S#coefficients)

net SOSPoly := S -> (
    if #gens S == 0 then return "0";
    return "coeffs:"||net coefficients S||"gens:"||net gens S;
    )

Number * SOSPoly := (a,S) -> (
    if a<0 then error "scalar must be nonnegative";
    if a==0 then return sosPoly(ring S, {}, {});
    return sosPoly(ring S, gens S, a * coefficients S);
    )

SOSPoly + SOSPoly := (S,S') -> (
    R := ring S;
    if R =!= ring S' then error "cannot add elements of different rings";
    return sosPoly(R,S#gens|S'#gens, S#coefficients|S'#coefficients);
    )

SOSPoly * SOSPoly := (g1,g2)-> (
    if g1#ring =!= g2#ring then error "cannot multiply elements of different rings";
    q1:=for i in g1#gens list(
        for j in g2#gens list i*j);
    q2:=for i in g1#coefficients list(
        for j in g2#coefficients list i*j);
    return sosPoly(g1#ring, flatten(q1),flatten(q2));
    )

SOSPoly ^ ZZ := (p1,D)->(
    if D<=0 then error "power should be a positive integer.";
    if odd D then error "power should be an even integer.";
    p2 := (sumSOS p1)^(D//2);
    return sosPoly(ring p1,{p2},{1});
    )

SOSPoly == RingElement := (S, f) -> (
    if ring S=!=ring f then
        error "Cannot compare elements of different rings. Try to use 'sub'.";
    return sumSOS S == f;
    )

RingElement == SOSPoly := (f, S) -> S == f

SOSPoly == SOSPoly := (S, S') -> S == sumSOS S'

sumSOS = method()

sumSOS (List, List) := (g,d) -> sum for i to #g-1 list g_i^2 * d_i

sumSOS SOSPoly := a -> sum for i to #(a#gens)-1 list a#gens_i^2 * a#coefficients_i

clean(RR,SOSPoly) := (tol,s) -> (
    if s===null then return (,);
    R := ring s;
    kk := coefficientRing R;
    if kk === QQ then tol=0.;
    g := gens s;
    d := coefficients s;
    I := positions(d, di -> di>tol);
    d = d_I;
    g = g_I;
    if kk =!= QQ then g = clean_tol \ g;
    return sosPoly(R,g,d);
    )

SDPResult = new Type of HashTable

sdpResult = (mon,Q,X,tval,mult) -> (
    new SDPResult from {
        Monomials => mon,
        GramMatrix => Q,
        MomentMatrix => X,
        "tval" => tval,
        "mult" => mult
        }
    )

net SDPResult := sol -> (
    mat2str := M -> 
        if M===null then "null" 
        else numRows M | "x" | numColumns M | " matrix over " | toString ring M;
    str := {
        {"MomentMatrix", mat2str sol#MomentMatrix},
        {"GramMatrix", mat2str sol#GramMatrix},
        {"Monomials", mat2str sol#Monomials}
        };
    mult := sol#"mult";
    if mult=!=null then
        str = append(str,{"Multipliers",mat2str mult});
    return netList(str,HorizontalSpace=>1,Alignment=>Center)
    )

readSdpResult = sol -> (sol#Monomials, sol#GramMatrix, sol#MomentMatrix, sol#"tval")

--##########################################################################--
-- METHODS
--##########################################################################--

verbose = (s,o) -> if o.Verbose then print s

--###################################
-- basicLinearAlgebra
--###################################

isExactField = kk -> (
    try (kk = ring kk);
    kk = ultimate(coefficientRing,kk);
    return precision 1_kk == infinity;
    )

linsolve = (A,b) -> (
    if isExactField A then return try solve(A,b);
    tol := 1e-12;
    x := solve(A,b,ClosestFit=>true);
    if norm(A*x-b) > tol then return;
    return x;
    )

truncatedSVD = (A,tol) -> (
    -- truncate small (or big) singular values
    (S, U, Vt) := SVD A;
    idx := positions(S, s->(s > abs tol));
    if tol>0 then(
        S = take(S,#idx);
        U = submatrix(U,,idx);
        Vt = submatrix(Vt,idx,);
    )else(
        S = drop(S,#idx);
        U = submatrix'(U,,idx);
        Vt = submatrix'(Vt,idx,); );
    return (S,U,Vt);
    )

kernelGens = A -> (
    if isExactField A then return gens kernel A;
    tol := 1e-12;
    (S,U,Vt) := truncatedSVD(A,-tol);
    return transpose Vt;
    )

zeros = (kk,m,n) -> map(kk^m,kk^n,{})

smat2vec = method( Options => {Scaling => 1} )
smat2vec(List) := o -> A -> (
    n := #A;
    v := for i to n-1 list
        for j from i to n-1 list 
            if i==j then A#i#j else o.Scaling*A#i#j;
    return flatten v;
    )
smat2vec(Matrix) := o -> A -> matrix(ring A, apply(smat2vec(entries A,o), a->{a}))

vec2smat = method( Options => {Scaling => 1} )
vec2smat(List) := o -> v -> (
    N := #v;
    n := (-1 + round sqrt(1+8*N))//2;
    ct := -1;
    L := for i to n-1 list (toList(i:0) |
        for j from i to n-1 list (ct = ct+1; ct));
    A := table(toList(0..n-1), toList(0..n-1), (i,j) -> 
        if i==j then v_(L#i#j) 
        else if i<j then v_(L#i#j)/(o.Scaling)
        else v_(L#j#i)/(o.Scaling) );
    return A;
    )
vec2smat(Matrix) := o -> v -> matrix(ring v, vec2smat(flatten entries v,o))

--###################################
-- solveSOS
--###################################

sosPoly(Matrix,Matrix) := (mon,Q) -> (
    if mon===null or Q===null then return;
    (L,D,P,err) := PSDdecomposition(Q);
    if err != 0 then (
        print "Gram Matrix is not positive semidefinite";
        return 
        );
    n := numRows Q;
    g := toList flatten entries (transpose mon * P * L);
    d := for i to n-1 list D_(i,i);
    idx := positions (d, i->i!=0);
    d = d_idx;
    g = g_idx;
    return sosPoly(ring mon,g,d);
    )
sosPoly(SDPResult) := sol -> if sol#GramMatrix=!=null then sosPoly(sol#Monomials, sol#GramMatrix)

-- internal way to call solveSOS
rawSolveSOS = method(
     Options => {RndTol => 3, Solver=>null, Verbose => false, TraceObj => false} )
 
rawSolveSOS(Matrix,Matrix,Matrix) := o -> (F,objP,mon) -> (
    -- f is a polynomial to decompose
    -- mon is a vector of monomials
    -- objFcn is a linear objective function for the SDP solver

    checkInputs := (mon) -> (
        if numColumns mon > 1 then error("Monomial vector should be a column.");
        isMonomial := max(length \ terms \ flatten entries mon)==1;
        if not isMonomial then error("Vector must consist of monomials.");
        );

    checkInputs(mon);
    kk := coefficientRing ring F;
         
    -- build SOS model --     
    (C,Ai,p0,V,A,B,b) := createSOSModel(F,mon,Verbose=>o.Verbose);

    ndim := numRows C;
    np := numRows objP;

    obj := 
        if o.TraceObj then
            map(RR^(#Ai),RR^1,(i,j)-> -trace Ai#i)
        else
            -(transpose V * objP); 
    if obj==0 then verbose( "Solving SOS feasibility problem...", o)
    else verbose("Solving SOS optimization problem...", o);

    (my,X,Q) := solveSDP(C, Ai, obj, Solver=>o.Solver, Verbose=>o.Verbose);
    if Q===null then return sdpResult(mon,Q,X,,);
    y := -my;
    pvec0 := flatten entries(p0 + V*y);

    if not isExactField kk then return sdpResult(mon,Q,X,pvec0,);
    if o.RndTol==infinity then
        return sdpResult(changeMatrixField(RR,mon),Q,X,pvec0,);

    -- rational rounding --
    (ok,Qp,pVec) := roundSolution(pvec0,Q,A,B,b,o.RndTol);
    if ok then return sdpResult(mon,Qp,X,pVec,);
    print "rounding failed, returning real solution";
    return sdpResult(changeMatrixField(RR,mon),Q,X,pvec0,);
    )

rawSolveSOS(Matrix,Matrix) := o -> (F,objP) -> (
    mon := choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,,,);
    return rawSolveSOS(F,objP,mon,o);
    )
rawSolveSOS(Matrix) := o -> (F) -> 
    rawSolveSOS(F,zeros(QQ,numRows F-1,1),o)

-- This is the main method to decompose a polynomial as a 
-- sum of squares using an SDP solver.
solveSOS = method(
     Options => {RndTol => 3, Solver=>null, Verbose => false, TraceObj => false} )

solveSOS(RingElement,RingElement,Matrix) := o -> (f,objFcn,mon) -> (
    (F,objP) := parameterVector(f,objFcn);
    return rawSolveSOS(F,objP,mon,o);
    )
solveSOS(RingElement,Matrix) := o -> (f,mon) -> 
    solveSOS(f,0_(ring f),mon,o)

solveSOS(RingElement,RingElement) := o -> (f,objFcn) -> (
    (F,objP) := parameterVector(f,objFcn);
    mon := choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,mon,,);
    return rawSolveSOS(F,objP,mon,o);
    )
solveSOS(RingElement) := o -> (f) -> 
    solveSOS(f,0_(ring f),o)


changeRingField = (kk,R) -> kk(monoid[gens R])

changeMatrixField = (kk, M) -> (
    -- M is a matrix whose entries are polynomials whose coefficient
    -- ring should be changed.
    R := changeRingField(kk, ring M);
    return matrix applyTable(entries M, m -> toRing(R,m));
    )

toRing = method ()
toRing (Ring, RingElement) := (S,f) -> (
    -- maps f to ring S
    R := ring f;
    kk := coefficientRing R;
    phi := map(S,R);
    -- QQ => RR
    if kk===QQ then return phi(f);
    -- RR => QQ
    if not instance(kk, RealField) then error "Expecting conversion from real here";
    (mon,coef) := coefficients f;
    mon = matrix {liftMonomial_S \ flatten entries mon};
    K := 2^(precision kk);
    coef = matrix(QQ, {for c in flatten entries coef list round(K*sub(c,RR))/K});
    f' := (mon*transpose coef)_(0,0);
    return f';
    )

toRing (Ring, SOSPoly) := (S, s) -> (
    -- maps s to Ring S
    R := ring s;
    kk := coefficientRing R;
    -- QQ => RR
    if kk===QQ then 
        return sosPoly (S, (x -> sub (x, S)) \ gens s,
            (q -> sub (q, kk)) \ coefficients s);
    -- RR => QQ
    if not (instance (kk, RealField) and coefficientRing S===QQ) then 
        error "Error: only conversion between real and rational coefficient fields is implemented.";
	g' := toRing_S \ gens s;
	K := 2^(precision kk);
	c' := for c in coefficients s list round(K*sub(c,RR))/K;
	return sosPoly (S, g', c')
    )

liftMonomial = (S,f) -> (
    -- maps monomial f to ring S
    n := numgens S;
    e := first exponents f;
    e = e_(toList(0..n-1)); -- ignore some variables
    return S_e;
    )

roundSolution = {Verbose=>false} >> o -> (pvec0,Q,A,B,b,RndTol) -> (
    -- round and project --
    Qnum := matrix applyTable(entries Q, a -> round(a*2^52)/2^52);

    dhi := 52;
    d := RndTol;
    np := numColumns B;
    
    while (d < dhi) do (
        verbose("rounding step #" | d, o);
        if np!=0 then (
           pVec := map(QQ^np,QQ^1,(i,j) -> round(pvec0_i*2^d)/2^d);
           bPar := b - B*pVec;
           ) 
        else bPar= b;

        (ok,Qp) := roundPSDmatrix(Qnum,A,bPar,d,Verbose=>o.Verbose);
        if ok then break else d = d + 1;
        );
    pVec = if np!=0 then flatten entries pVec else null;
    return (ok,Qp,pVec);
    )

createSOSModel = method(
    Options => {Verbose => false} )
createSOSModel(RingElement,Matrix) := o -> (f,v) -> (
    F := parameterVector(f);
    return createSOSModel(F,v);
    )
createSOSModel(Matrix,Matrix) := o -> (F,v) -> (
    kk := coefficientRing ring F;
    np := numRows F - 1;
    n := numRows v;
    
    -- monomials in vvT
    vvT := entries(v* transpose v);
    mons := g -> set first entries monomials g;
    K1 := toList \\ sum \\ mons \ flatten vvT;

    -- monomials in F and not in vvT
    lmf := sum \\ mons \ flatten entries F;
    K2 := toList(lmf - K1);
    K := K1 | K2;
    
    -- Linear constraints: b
    b := map(kk^#K, kk^1, (i,j) -> coefficient(K#i,F_(0,0)) );
    
    -- Linear constraints: A, B
    coeffMat := (x,A) -> applyTable(A, a -> coefficient(x,a));
    A := matrix(kk, for i to #K1-1 list smat2vec(coeffMat(K1_i, vvT),Scaling=>2) );
    A = A || zeros(kk,#K2,n*(n+1)//2);
    
    -- Consider search-parameters:
    B := map(kk^#K, kk^np, (i,j) -> -coefficient(K#i, F_(j+1,0)) );
    
    (C,Ai,p0,V) := getImageModel(A,B,b);
    
    return (C,Ai,p0,V,A,B,b);
    )

getImageModel = (A,B,b) -> (
    -- given the affine subspace {Aq + Bp = b}
    -- find a parametrization of the form
    -- Q = C + sum_i ti Ai
    -- p = p0 + V t

    -- compute the C matrix
    n1 := numColumns A;
    n2 := numColumns B;
    AB := A|B;
    x := linsolve(AB,b);
    c := x^{0..n1-1};
    p0 := x^{n1..n1+n2-1};
    C := vec2smat(c);
    
    -- compute the A_i matrices     
    W := - kernelGens AB;
    r := numColumns W;
    U := W^{0..n1-1};
    V := W^{n1..n1+n2-1};
    Ai := toSequence for k to r-1 list vec2smat(U_{k});

    return (C,Ai,p0,V);
    )

parameterVector = method()
parameterVector(RingElement,RingElement) := (f,objFcn) -> (
    -- given a polynomial f = f_0 + \sum_i p_i * f_i
    -- the method returns the vector with the f_i's
    R := ring f;
    kk := coefficientRing R;
    if isField kk then (
        if objFcn!=0 then 
            error "Objective must be zero if there are no parameters.";
        return (matrix{{f}},zeros(kk,0,1));
        );
    if first degree f > 1 then 
        error("Polynomial should depend affinely on the parameters.");
    p := gens R;
    F := matrix for t in {1_R}|p list {coefficient(t,f)};
    if degree objFcn > {1,0} then 
        error("Objective should be a linear function of the parameters.");
    kk = coefficientRing kk;
    objP := matrix for t in p list {sub(coefficient(t,objFcn),kk)};
    return (F,objP);
    )
parameterVector(RingElement) := (f) -> first parameterVector(f,0_(ring f))

choosemonp = method(
    Options => {Verbose => false} )
choosemonp(RingElement) := o -> (f) -> (
    F := parameterVector(f);
    mon := choosemonp(F);
    if mon===null then return;
    return sub(mon,ring f);
    )
choosemonp(Matrix) := o -> (F) -> (
     R := ring F;
     if F==0 then error "Expected nonzero inputs.";
     if isQuotientRing R then error("Monomial vector must be provided in quotient rings.");
     n:= numgens R;
     mons := g -> set first entries monomials g;
     lm0 := mons F_(0,0);
     filterVerts := (verts) -> (
         -- only consider those without parameters (this is a hack!)
         return select(verts, v -> member(R_v,lm0));
         );
     lmf := sum \\ mons \ flatten entries F;
     falt := sum lmf;
     
     -- Get exponent-points for Newton polytope:
     points := substitute(matrix (transpose exponents falt),QQ);
     maxdeg := first degree falt;
     mindeg := floor first min entries (transpose points*matrix map(ZZ^n,ZZ^1,i->1));
     maxdegs := apply(entries points, i-> max i);
     mindegs := apply(entries points, i-> min i);
     
     -- Regard exponent-points in a possible subspace
     numpoints := numColumns points;
     shift := first entries transpose points;
     V := matrix transpose apply(entries transpose points, i -> i - shift);
     basV := mingens image V;
     basVdim := numgens image basV;
     if basVdim != n then T := id_(QQ^n)//basV else T = id_(QQ^n);
     basVtrans := kernelGens transpose basV;
     
     -- Compute Newton polytope:
     liftedpts := T*V || map (QQ^1,QQ^(size falt),i->1);
     dualpolytope := transpose substitute(first fourierMotzkin(liftedpts),QQ);
     bidual := first fourierMotzkin transpose dualpolytope;
     polytope := basV * matrix drop(entries bidual,-1);
     polytope = matrix transpose apply(entries transpose polytope, i -> i + shift);
     polytope = sub(polytope,ZZ);
     oddverts := select(entries transpose polytope, i->any(i,odd));
     if #filterVerts(oddverts)>0 then(
         print("Newton polytope has odd vertices. Terminate.");
         return;
         );

     -- Get candidate points
     cp := pointsInBox(mindeg,maxdeg,mindegs,maxdegs);
     verbose("#candidate points: " | #cp, o);
     -- Only the even ones
     cpf := select(cp,i-> all(i,even)); 
     verbose("#even points: " | #cpf, o);
     -- Drop points that do not live on the subspace: 
     cpf2 := select(cpf,i-> matrix{i-shift}*basVtrans==0);
     verbose("#points in subspace of exponent-points: " | #cpf2, o);
     
     -- Find points within the polytope:
     lexponents := select(cpf2, i-> 
           max flatten entries (dualpolytope * ((T * transpose matrix {i-shift})||1)) <=0)/2;
     isInteger := l -> denominator l == 1;
     assert all(flatten lexponents, isInteger );
     lexponents = apply(lexponents, i -> numerator \ i);
     lmSOS := for i in lexponents list R_i;
     verbose("#points inside Newton polytope: " | #lmSOS, o);

     if #lmSOS==0 then return;
     return matrix transpose {lmSOS};
     )

pointsInBox = (mindeg,maxdeg,mindegs,maxdegs) -> (
    -- integer vectors within specified bounds
    n := #mindegs;
    -- Get candidate points
    local x; x= symbol x;
    R0 := QQ(monoid[x_0..x_(n-1)]);
    mon := flatten entries basis(mindeg,maxdeg,R0);
    e := apply (mon, i -> flatten exponents i);
    -- Only those within the box of degrees[mindegs:maxdegs]:
    e = select(e,i-> all(i-mindegs,j->j>=0) and all(maxdegs-i,j->j>=0)); 
    return e;
    )

project2linspace = (A,b,x0) -> (
     -- cast into QQ (necessary class to compute inverse)
     A2 := promote (A,QQ);
     -- ugly hack to convert b into a matrix if it is a scalar in QQ/ZZ:
     b2 := promote (matrix{{b}},QQ);
     x02 := promote (x0,QQ);

     -- compute projection:
     xp := x02 - transpose(A2)*((A2*x02-b2)//(A2*transpose(A2)))
     )

roundPSDmatrix = {Verbose=>false} >> o -> (Q,A,b,d) -> (
     verbose("Rounding precision: " | d, o);
     Q0 := matrix (applyTable (entries Q, i -> round(i*2^d)/2^d) );
     x0 := smat2vec(Q0);
     t := timing (xp := project2linspace(A,b,x0););
     verbose("Time needed for projection: " | net t#0, o);
     Q = vec2smat(xp);

     t = timing((L,D,P,Qpsd) := PSDdecomposition(Q););
     verbose("Time needed for LDL decomposition: " | net t#0, o);
     if Qpsd == 0 then (true, Q) else (false,Q)
     )

PSDdecomposition = (A) -> (
    kk := ring A;
    if isExactField kk then
        return LDLdecomposition(A);
    if kk=!=RR and not instance(kk,RealField) then
        error "field must be QQ or RR";
    tol := 1e-9;
    (e,V) := eigenvectors(A,Hermitian=>true);
    err := if all(e, i -> i > -tol) then 0 else 1;
    e = max_0 \ e;
    D := diagonalMatrix e;
    P := id_(kk^(numRows A));
    return (V,D,P,err);
    )
    
LDLdecomposition = (A) -> (
     kk := ring A;
     if kk=!=QQ and kk=!=RR and not instance(kk,RealField) then
        error "field must be QQ or RR";
     if transpose A != A then error("Matrix must be symmetric.");
     tol := if isExactField kk then 0 else 1e-9;

     n := numRows A;
     Ah := new MutableHashTable; map (kk^n,kk^n,(i,j)->Ah#(i,j) = A_(i,j));
     v := new MutableList from for i to n-1 list 0_kk;
     d := new MutableList from for i to n-1 list 0_kk;
     piv := new MutableList from toList(0..n-1);
     err := 0;

     for k from 0 to n-1 do (
      q := k + maxPosition apply(k..n-1, i->Ah#(i,i));
      -- Symmetric Matrix Permutation:
      tmp := piv#q; piv#q = piv#k; piv#k = tmp;
      for i to n-1 do (tmp := Ah#(i,q); Ah#(i,q) = Ah#(i,k); Ah#(i,k) = tmp;);
      for i to n-1 do (tmp := Ah#(q,i); Ah#(q,i) = Ah#(k,i); Ah#(k,i) = tmp;);

      --  positive semidefinite?
      if Ah#(k,k) < -tol then (err = k+1; break;);
      if abs(Ah#(k,k))<=tol then 
          if any(0..n-1, i->abs(Ah#(i,k))>tol) then (
               err = k+1; break;);

      -- Perform LDL factorization step:
      if Ah#(k,k) > 0 then (
                 for i to k-1 do (v#i = Ah#(k,i)*Ah#(i,i));
           Ah#(k,k) = Ah#(k,k) - sum for i to k-1 list Ah#(k,i)*v#i;
           if Ah#(k,k) < -tol then (err = k+1; break;);
           if Ah#(k,k) > 0 then
             for i from k+1 to n-1 do
                 Ah#(i,k) = (Ah#(i,k)-sum for j to k-1 list Ah#(i,j)*v#j) / Ah#(k,k);
      );
     );

     A = map(kk^n,kk^n,(i,j)-> if i>j then Ah#(i,j) else if i==j then 1_kk else 0_kk);
     D := map(kk^n,kk^n,(i,j)->if i==j then Ah#(i,j) else 0_kk);
     P := submatrix(id_(kk^n),toList piv);
     (A,D,P,err)
)

--###################################
-- SOS IDEAL
--###################################

makeMultiples = (h, D, homog) -> (
    -- h is a list of polynomials
    -- multiplies each hi with monomials up to degree D
    if #h==0 then return ({},{});
    if D < max\\first\degree\h then
        error "increase degree bound";
    R := ring h#0;
    -- compute monomials
    mon := for i to #h-1 list (
        di := D - first degree h#i;
        b := if homog then basis(di,R) else basis(0,di,R);
        flatten entries b
        );
    H := for i to #h-1 list h#i * mon#i;
    return (flatten H, mon);
    )

sosInIdeal = method(
     Options => {RndTol => 3, Solver=>"CSDP", Verbose => false} )
sosInIdeal (Ring, ZZ) := o -> (R,D) -> (
    -- find sos polynomial in a quotient ring
    if odd D then error "D must be even";
    mon := if isHomogeneous R then transpose basis(D//2,R)
        else transpose basis(0,D//2,R);
    (mon',Q,X,tval) := readSdpResult solveSOS (0_R, mon, o);
    if Q===null or Q==0 or (not isExactField Q and norm Q<1e-6) then (
        print("no sos polynomial in degree "|D);
        return sdpResult(mon,,X,,);
        );
    return sdpResult(mon,Q,X,tval,);
    )
sosInIdeal (Matrix,ZZ) := o -> (h,D) -> (
    -- h is a row vector of polynomials
    -- D is a degree bound
    -- returns sos polynomial in <h>

    -- The output is an SDPResult with the multipliers that
    -- express the SOS in terms of the generators.
    
    if numRows h > 1 then error "h must be a row vector";
    if odd D then error "D must be even";
    homog := isHomogeneous h;
    (H,m) := makeMultiples(first entries h, D, homog);
    F := matrix transpose {{0}|H};
    (mon,Q,X,tval) := readSdpResult rawSolveSOS (F, o);
    if Q===null or Q==0 or (not isExactField Q and norm Q<1e-6) then (
        print("no sos polynomial in degree "|D);
        return sdpResult(mon,,X,,);
        );
    mult := getMultipliers(m,tval,ring mon);
    return sdpResult(mon,Q,X,tval,mult);
    )

getMultipliers = (mon,tval,S) -> (
    if #mon==0 then return zeros(S,0,1);
    mon = applyTable(mon, i -> sub(i,S));
    k := -1;
    mult := matrix(S, for m in mon list
        {sum for i in m list( k=k+1; i*tval#k)} );
    return mult;
    )

sosdecTernary = method(
     Options => {RndTol => 3, Solver=>"CSDP", Verbose => false} )
sosdecTernary(RingElement) := o -> (f) -> (
    -- Implements Hilbert's algorithm to write a non-negative ternary
    -- form as sos of rational functions.
    -- Returns two lists of SOSPolys, the numerator and the denomenator polys
    if numgens ring f =!= 3 then error "polynomial must involve 3 variables";
    if not isHomogeneous f then error "polynomial must be homogeneous";
    fi := f;
    S := {};
    di := first degree fi;
    while di > 4 do(
        sol := sosInIdeal(matrix{{fi}},2*di-4,o);
        Si := sosPoly sol;
        if Si===null then return (,);
        mult := sol#"mult";
        fi = mult_(0,0);
        if fi==0 then return (,);
        di = first degree fi;
        S = append(S,Si);
        );
    (mon,Q,X,tval) := readSdpResult rawSolveSOS matrix{{fi}};
    if Q===null or Q==0 or (not isExactField Q and norm Q<1e-6) then return (,);
    Si = sosPoly(mon,Q);
    if Si===null then return (,);
    S = append(S,Si);
    nums := for i to #S-1 list if odd i then continue else S#i;
    dens := for i to #S-1 list if even i then continue else S#i;
    return (nums, dens);
    )

--###################################
-- SOS OPTIMIZATION
--###################################

recoverSolution = (mon,X) -> (
    if X===null then return {};
    e := eigenvalues(X,Hermitian=>true);
    if e#(-1)<=0 or e#0/e#(-1) < -1e-9 then 
        error "Moment matrix is not positive semidefinite";
    i0 := position(flatten entries mon, i -> i==1);
    if i0===null then
        error "The monomial vector must contain 1";
    if e#(-2) > 1e-4 then 
        print "Moment matrix is not rank one, solution might not be correct.";
    sol := for i to numRows mon-1 list (
        y := mon_(i,i0);
        if sum degree y!=1 then continue;
        y => X_(i,i0) );
    return sol;
    )

-- Unconstrained minimization 
-- sos lower bound for the polynomial f
lowerBound = method(
     Options => {RndTol => 3, Solver=>null, Verbose => false} )
lowerBound(RingElement) := o -> (f) -> lowerBound(f,-1,o)
lowerBound(RingElement,ZZ) := o -> (f,D) -> lowerBound(f,zeros(ring f,1,0),D,o)

-- Minimize a polynomial on an algebraic set
lowerBound(RingElement,Matrix,ZZ) := o -> (f,h,D) -> (
    -- Lasserre hierarchy for the problem
    -- min f(x) s.t. h(x)=0
    numdens := (f) -> (
        R := ring f;
        (num,den) := (f, 1_R);
        if isField R then(
            (num,den) = (numerator f, denominator f);
            R = last R.baseRings;
            );
        return (R,num,den)
        );
    checkInputs := (D,num,den,h,R) -> (
        if numRows h > 1 then error "h must be a row vector";
        if D<0 then(
            if numColumns h>0 or isQuotientRing R then
                error "Degree bound must be provided"
        )else(
            if odd D then error "degree bound must be even";
            );
        );
    
    (R,num,den) := numdens(f);
    checkInputs(D,num,den,h,R);
    -- prepare input
    (H,m) := makeMultiples(flatten entries h, D, false);
    F := matrix transpose {{num,-den}|H};
    objP := matrix{{-1}} || zeros(ZZ,#H,1);

    -- call solveSOS
    o' := new OptionTable from
        {RndTol=>o.RndTol, Solver=>o.Solver, Verbose=>o.Verbose};
    mon := if isQuotientRing R then transpose basis(0,D//2,R)
        else choosemonp (F,Verbose=>o.Verbose);
    if mon===null then return (,sdpResult(,,,,));
    (mon',Q,X,tval) := readSdpResult rawSolveSOS(F,objP,mon,o');
    if tval=!=null then(
        bound := tval#0;
        mult := getMultipliers(m,drop(tval,1),ring mon');
    )else (bound,mult) = (,);
    return (bound,sdpResult(mon,Q,X,tval,mult));
    )

--###################################
-- SDP SOLVER
--###################################

solveSDP = method(
     Options => {Solver=>null, Verbose => false} )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C,sequence A,b,o)

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C,sequence A,b,y,o)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    if numRows b=!=#A then error "Bad matrix dimensions.";
    (C,A,b) = toReal(C,A,b);
    (ok,y,X,Z) := (,,,);
    (ok,y,X,Z) = sdpNoConstraints(C,A,b);
    if ok then return (y,X,Z);
    solver := chooseSolver o;
    if solver == "M2" then(
        (ok,y,X,Z) = trivialSDP(C,A,b);
        if ok then return (y,X,Z)
        else (y,Z) = simpleSDP(C,A,b,Verbose=>o.Verbose)
        )
    else if solver == "CSDP" then
        (y,X,Z) = solveCSDP(C,A,b,Verbose=>o.Verbose)
    else if solver == "SDPA" then
        (y,X,Z) = solveSDPA(C,A,b,Verbose=>o.Verbose)
    else
        error "unknown SDP solver";
    ntries := 10;
    (y,Z) = findNonZeroSolution(C,A,b,o,y,Z,ntries);
    return (y,X,Z);
    )

findNonZeroSolution = (C,A,b,o,y,Z,ntries) -> (
    if not(C==0 and b==0 and y==0) then return (y,Z);
    print "Zero solution obtained. Trying again.";
    m := numRows b;
    badCoords := set();
    iszero := a -> norm a < 1e-8;
    for i to ntries-1 do(
        if #badCoords==m then break;
        b' := map(RR^m,RR^1, (j,l) -> 
            if member(j,badCoords) then 0 else random(RR)-.5 );
        (y',X',Z') := solveSDP(C,A,b',o);
        if Z'=!=null and not iszero Z' then return (y',Z');
        if X'===null and y'=!=null and (transpose b' * y')_(0,0) < -.1 then
            badCoords = badCoords + set select(0..m-1, j -> not iszero y'_(j,0));
        );
    return (y,Z);
    )

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y0) -> (
    (C,A,b) = toReal(C,A,b);
    y0 = promote(y0,RR);
    (ok,y,X,Z) := (,,,);
    (ok,y,X,Z) = sdpNoConstraints(C,A,b);
    if ok then return (y,X,Z);
    if chooseSolver o != "M2" then return solveSDP(C,A,b,o);
    (ok,y,X,Z) = trivialSDP(C,A,b);
    if ok then return (y,X,Z);
    (y,Z) = simpleSDP(C,A,b,y0,Verbose=>o.Verbose);
    return (y,,Z);
    )

chooseSolver = o -> (
    if o.Solver=!=null then return o.Solver;
    if defaultSolver=!=null then return defaultSolver;
    if csdpexec=!=null then return "CSDP";
    if sdpaexec=!=null then return "SDPA";
    return "M2";
    )

toReal = (C,A,b) -> (
    C = promote(C,RR);
    A = apply(A, Ai -> promote(Ai,RR));
    b = promote(b,RR);
    return (C,A,b);
    )

sdpNoConstraints = (C,A,b) -> (
    tol := 1e-10;
    if #A==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=-tol then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, y0, 0*C, C);
        )else(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
    )

-- check trivial cases
trivialSDP = (C,A,b) -> (
    if #A==0 or b==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=0 then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, y0, 0*C, C);
        )else if #A==0 then(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
    )

--simpleSDP

simpleSDP = method(
    TypicalValue => Matrix,
    Options => {UntilObjNegative => false, Verbose => false} )

simpleSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    R := RR;
    n := numRows C;

    -- try to find strictly feasible starting point --
    (y,Z) := (,);
    lambda := min eigenvalues (C, Hermitian=>true);
    if lambda > 0 then
        y = zeros(R,#A,1)
    else(
        verbose("Computing strictly feasible solution...", o);
        y =  zeros(R,#A,1) || matrix{{lambda*1.1}};
        obj :=  zeros(R,#A,1) || matrix{{-1_R}};
        (y,Z) = simpleSDP(C,append(A,id_(R^n)), obj, y, UntilObjNegative=>true, Verbose=>o.Verbose);
        if y===null then return (,);
        y = transpose matrix {take (flatten entries y,numRows y - 1)};
        );
    verbose("Computing an optimal solution...", o);
    return simpleSDP(C, A, b, y, o);
    )

simpleSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y) -> (
    print "Running M2 solver";
    R := RR;
    n := numgens target C;

    m := numgens target y;
    mu := 1_R;
    theta := 10_R;
    iter := 1;
    NewtonIterMAX := 40;

    verbose("#It:       b'y      dy'Hdy   mu   alpha", o);

    while mu > 0.000001 do (
        mu = mu/theta;
        while true do (
            S := C - sum toList apply(0..m-1, i-> y_(i,0) * A_i);
            try Sinv := solve(S, id_(target S)) else (
                print "slack matrix is singular, terminate";
                return (,) );
            -- compute Hessian:
            H := map(R^m,R^m,(i,j) -> trace(Sinv*A_i*Sinv*A_j));
            if H==0 then (
                print "Hessian is zero";
                return (,) );
            -- compute gradient:
            g := map(R^m,R^1,(i,j) -> b_(i,0)/mu + trace(Sinv*A_i));
            
            -- compute damped Newton step:
            dy := -solve(H,g,ClosestFit=>true);
            alpha := backtrack(S, -sum for i to m-1 list matrix(dy_(i,0) * entries A_i));
            if alpha===null then return (,);
            y = y + transpose matrix {alpha* (flatten entries dy)};
            lambda := (transpose dy*H*dy)_(0,0);
            obj := transpose b * y;
            
            -- print some information:
            verbose(iter | ":  " | net obj | "    " | net lambda | "    " | net mu | "    " | net alpha, o);

            iter = iter + 1;
            if iter > NewtonIterMAX then (
                verbose("Warning: exceeded maximum number of iterations", o);
                break);
            if o.UntilObjNegative and (obj_(0,0) < 0) then break;
            if lambda < 0.4 then break;
            ); 
        );
    Z := C - sum(for i to #A-1 list y_(i,0) * A_i);
    return (y,Z);
    )     

backtrack = args -> (
     S0 := args#0;
     R := ring S0;
     dS := args#1;
     alpha := 1_R;
     BacktrackIterMAX := 100;
     S :=  matrix( alpha * entries dS) + S0;
     
     cnt := 1;     
     while min eigenvalues(S,Hermitian=>true) <= 0 do (
      cnt = cnt + 1;
      alpha = alpha / sqrt(2_R);
      S = S0 + matrix( alpha * entries dS);
      if cnt > BacktrackIterMAX then (
          print ("line search did not converge.");
          return null );
      );
     return alpha;
     )


--solveCSDP

solveCSDP = method( Options => {Verbose => false} )
solveCSDP(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    -- CSDP expects the file fparam to be in the working directory.
    -- That's why we need to change directory before executing csdp.
    if csdpexec===null then error "csdp executable not found";
    n := numColumns C;
    fin := getFileName ".dat-s";
    (dir,fin1) := splitFileName(fin);
    fparam := dir | "param.csdp";
    fout := getFileName "";
    fout2 := getFileName "";
    writeSDPA(fin,C,A,b);
    writeCSDPparam(fparam);
    print("Executing CSDP on file " | fin);
    r := run("cd " | dir | " && " | csdpexec | " " | fin1 | " " | fout | ">" | fout2);
    if r == 32512 then error "csdp executable not found";
    print("Output saved on file " | fout);
    (y,X,Z) := readCSDP(fout,fout2,n,o.Verbose);
    y = checkDualSol(C,A,y,Z,o.Verbose);
    return (y,X,Z);
)

getFileName = (ext) -> (
     filename := temporaryFileName() | ext;
     while fileExists(filename) do filename = temporaryFileName();
     return filename
)

splitFileName = (fname) -> (
    s := separate("/",fname);
    dir := demark("/",drop(s,-1))|"/";
    file := last s;
    return (dir,file);
    )

writeSDPA = (fin,C,A,b) -> (
    digits := 16;
    m := length A;
    n := numColumns C;
    A = prepend(C,A);
    f := openOut fin;
    inputMatrix := l -> (
        a := -A_l;
        pref := toString l | " 1 ";
        for i to n-1 do
            for j from i to n-1 do
                if a_(i,j)!=0 then
                    f << pref | toString(i+1) | " " | toString(j+1) | " " | format(digits,a_(i,j)) << endl;
    );
    f << "*SDPA file generated by SOSm2" << endl;    
    f << toString m << " =mdim" << endl;
    f << "1 =nblocks" << endl;
    f << toString n << endl;
    f << demark(" ", toString\flatten entries b) << endl;
    for l to m do(
        inputMatrix(l);
    );
    f << close;
)

writeCSDPparam = (fparam) -> (
    f := openOut fparam;
    f << "axtol=1.0e-8" << endl;
    f << "atytol=1.0e-8" << endl;
    f << "objtol=1.0e-8" << endl;
    f << "pinftol=1.0e8" << endl;
    f << "dinftol=1.0e8" << endl;
    f << "maxiter=100" << endl;
    f << "minstepfrac=0.90" << endl;
    f << "maxstepfrac=0.97" << endl;
    f << "minstepp=1.0e-8" << endl;
    f << "minstepd=1.0e-8" << endl;
    f << "usexzgap=1" << endl;
    f << "tweakgap=0" << endl;
    f << "affine=0" << endl;
    f << "printlevel=1" << endl;
    f << "perturbobj=0" << endl;
    f << "fastmode=0" << endl;
    f << close;
    )

readCSDP = (fout,fout2,n,Verbose) -> (
    sdpa2matrix := s -> (
        e := for i in s list (i_2-1,i_3-1) => i_4;
        e' := for i in s list (i_3-1,i_2-1) => i_4;
        return map(RR^n, RR^n, e|e');
        );
    readLine := l -> for s in separate(" ",l) list if s=="" then continue else value s;
    --READ SOLUTIONS
    tmp := getFileName "";
    r := run("cat " | fout | " | tr -d + > " | tmp);
    L := lines get openIn tmp;
    y := transpose matrix{readLine L_0};
    S := readLine \ drop(L,1);
    S1 := select(S, l -> l_0==1);
    S2 := select(S, l -> l_0==2);
    Z := sdpa2matrix(S1); -- slack matrix
    X := sdpa2matrix(S2); -- dual solution
    -- READ STATUS
    text := get openIn fout2;
    s := select(lines text, l -> match("Success",l));
    if #s==0 then( print "SDP failed"; return (,,) );
    s = first s;
    print if Verbose then text else s;
    if match("SDP solved",s) then null
    else if match("primal infeasible",s) then X=null
    else if match("dual infeasible",s) then (y=null;Z=null;)
    else error "unknown message";
    return (y,X,Z);
)

checkDualSol = (C,A,y,Z,Verbose) -> (
    if y===null then return;
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    if norm(Z-C+yA)<1e-5 then return y;
    if Verbose then print "updating dual solution";
    AA := transpose matrix(RR, smat2vec \ entries \ toList A);
    bb := transpose matrix(RR, {smat2vec entries(C-Z)});
    y = solve(AA,bb,ClosestFit=>true);
    return y;
    )

--solveSDPA

solveSDPA = method( Options => {Verbose => false} )
solveSDPA(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    if sdpaexec===null then error "sdpa executable not found";
    n := numColumns C;
    fin := getFileName ".dat-s";
    fout := getFileName "";
    writeSDPA(fin,C,A,b);
    print("Executing SDPA on file " | fin);
    r := run(sdpaexec | " " | fin | " " | fout | "> /dev/null");
    if r == 32512 then error "sdpa executable not found";
    if r == 11 then error "Segmentation fault running sdpa.";
    print("Output saved on file " | fout);
    (y,X,Z) := readSDPA(fout,n,o.Verbose);
    return (y,X,Z);
)

readSDPA = (fout,n,Verbose) -> (
    readVec := l -> (
        l = replace("([{} +])","",l);
        for s in separate(",",l) list if s=="" then continue else value s
    );
    readMatrix := ll -> 
        matrix for l in ll list readVec l;
    text := get openIn fout;
    L := lines text;
    --READ SOLUTIONS
    y := null; X := null; Z := null;
    i := position(L, l -> match("xVec =",l));
    if i=!=null then 
        y = transpose matrix {readVec L#(i+1)};
    i = position(L, l -> match("xMat =",l));
    if i=!=null then 
        Z = matrix for j to n-1 list readVec L#(i+j+2);
    i = position(L, l -> match("yMat =",l));
    if i=!=null then 
        X = matrix for j to n-1 list readVec L#(i+j+2);
    --READ STATUS
    if Verbose then print text;
    s := first select(L, l -> match("phase.value",l));
    if match("pdOPT",s) or match("pdFEAS",s) then 
        print "SDP solved"
    else if match("dUNBD",s) then(
        print "dual infeasible";
        y=null;Z=null; )
    else if match("pUNBD",s) then(
        print "primal infeasible";
        X=null; )
    else( 
        print("Warning: Solver returns unknown message!!! " |s);
        );
    return (y,X,Z);
    )

--###################################
-- Methods for testing
--###################################

checkSolver = method()
checkSolver(String,String) := (solver,fun) -> (
    checkMethod := hashTable {
        "solveSDP" => checkSolveSDP,
        "solveSOS" => checkSolveSOS,
        "sosdecTernary" => checkSosdecTernary,
        "sosInIdeal" => checkSosInIdeal,
        "lowerBound" => checkLowerBound
        };
    if checkMethod#?fun then 
        return checkMethod#fun(solver);
    if fun != "AllMethods" then
        error "No test implemented for this function";
    T := for f in keys checkMethod list(
        print "################################";
        print("checking method "|f);
        print "################################";
        t := checkMethod#f(solver);
        {f, testsString t}
        );
    print "################################";
    print("Summary");
    print netList T;
    )
checkSolver(String,Function) := (solver,fun) -> checkSolver(solver,toString fun)
checkSolver(String) := (solver) -> checkSolver(solver,"AllMethods")

-- A method to inform about the results of the tests in one function
testsString = t -> concatenate apply(t, i -> if i then " ✓ " else " ✘ ")
informAboutTests = t -> (
    print("Test Results: " | testsString t);
    )

--checkSolveSDP
checkSolveSDP = solver -> (
    tol := .001;
    equal := (y0,y) -> y=!=null and norm(y0-y)<tol*(1+norm(y0));
    local C; local b; local A; local A1; local A2; local A3; 
    local y0; local y; local X; local Z; local yopt;
    ---------------TEST0---------------
    C = matrix{{0,2,0,0,0,0},{2,0,0,0,0,0},
     {0,0,10,0,0,0},{0,0,0,10,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    A1 = matrix{{-1,0,0,0,0,0},{0,0,0,0,0,0},
     {0,0,1,0,0,0},{0,0,0,0,0,0},{0,0,0,0,-1,0},{0,0,0,0,0,0}};
    A2 = matrix{{0,0,0,0,0,0},{0,-1,0,0,0,0},
     {0,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1}};
    A = (A1,A2);
    y0 = matrix{{7},{9}};
    b = matrix{{1},{1}};
    (y,X,Z) = solveSDP(C,A,b,y0,Solver=>solver);
    yopt = matrix{{2.},{2.}};
    t0 := equal(yopt,y);
    ---------------TEST1---------------
    C = matrix {{2,1,-1},{1,0,0},{-1,0,5}};
    A1 = matrix {{0,0,1/2},{0,-1,0},{1/2,0,0}};
    A2 = matrix {{1,0,0},{0,1,0},{0,0,1}};
    A = (A1,A2);
    b = matrix {{0},{-1}};
    y0 = matrix {{0},{-.486952}};
    (y,X,Z) = solveSDP(C,A,b,y0,Solver=>solver);
    yopt = matrix{{1.97619},{.466049}};
    t1 := equal(yopt,y);
    ---------------TEST2---------------
    C = matrix{{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
    A1 = matrix{{-1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    A2 = matrix{{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
    A = (A1,A2);
    b = matrix{{1},{0}};
    (y,X,Z) = solveSDP(C,A,b,Solver=>solver);
    yopt = matrix{{0.},{4.}};
    t2 := equal(yopt,y); 
    ---------------TEST3---------------
    -- solution not strictly feasible
    C = matrix {{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
    A1 = matrix {{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
    A = sequence A1;
    b = matrix {{1}};
    (y,X,Z) = solveSDP(C,A,b,Solver=>solver);
    yopt = 4.;
    t3 := equal(yopt,y);
    ---------------TEST4---------------
    -- zero objective function
    C = matrix(RR, {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    A1 = matrix(RR, {{1, 3/2, 3/2}, {3/2, 0, 1/2}, {3/2, 1/2, 0}});
    A2 = matrix(RR, {{0, 1/2, 3/2}, {1/2, 0, 3/2}, {3/2, 3/2, 1}});
    A3 = matrix(RR, {{0, 0, 1/2}, {0, -1, 0}, {1/2, 0, 0}});
    A = (A1,A2,A3);
    b = matrix(RR, {{0}, {0}, {0}});
    (y,X,Z) = solveSDP(C, A, b, Solver=>solver);
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    t4 := norm(Z-C+yA)<1e-5;
    -----------------------------------
    test := {t0,t1,t2,t3,t4};
    informAboutTests test;
    -- trivial cases
    (y,X,Z) = solveSDP (matrix{{1,0},{0,-1}},(),zeros(QQ,0,1),Solver=>solver);
    assert(y===null and X===null);
    (y,X,Z) = solveSDP (matrix{{1,0},{0,1}},(),zeros(QQ,0,1),Solver=>solver);
    assert(y==0);
    return test;
    )

--checkSolveSOS
checkSolveSOS = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local w; w= symbol w;
    local t; t= symbol t;
    isGram := (f,mon,Q) -> (
        if Q===null then return false;
        e := eigenvalues(Q,Hermitian=>true);
        tol := 1e-8;
        if min e < -tol then return false;
        S := ring mon;
        if isExactField S then
            return f == transpose(mon)*Q*mon;
        return norm(sub(f,S) - transpose(mon)*Q*mon) < tol;
        );
    isGramParam := (f,mon,Q,tval) ->
        if tval===null then false else isGram(sub(f,t=>tval#0),mon,Q);
    ---------------GOOD CASES1---------------
    -- Test 0
    R := QQ[x,y];
    f := 4*x^4+y^4;
    (mon,Q,X,tval) := readSdpResult solveSOS(f,Solver=>solver);
    t0 := isGram(f,mon,Q);

    -- Test 1
    f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
    (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
    t1 := isGram(f,mon,Q);

    -- Test 2
    R = QQ[x,y,z];
    f = x^4+y^4+z^4-4*x*y*z+x+y+z+3;
    (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
    t2 := isGram(f,mon,Q);
    
    -- Test 3
    R = QQ[x,y,z,w];
    f = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
    (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
    t3 := isGram(f,mon,Q);

    -- Test 4 (parametric)
    R = QQ[x][t];
    f = (t-1)*x^4+1/2*t*x+1;
    (mon,Q,X,tval) = readSdpResult solveSOS (f,Solver=>solver);
    t4 := isGramParam(f,mon,Q,tval);

    ---------------BAD CASES1---------------
    -- Test 5
    R = QQ[x,y][t];
    f = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1; --Motzkin
    (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver); 
    t5 := ( Q === null );

    -- Test 6
    (mon,Q,X,tval) = readSdpResult solveSOS(f-t,-t, Solver=>solver); 
    t6 := ( Q === null );

    results := {t0,t1,t2,t3,t4,t5,t6};
    informAboutTests (results);
    return results
    )

-- check sosdecTernary
checkSosdecTernary = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;

    cmp := (f,p,q) -> (
        if p===null then return false;
        d := product(sumSOS\p) - f*product(sumSOS\q);
        if isExactField f then return d==0;
        return norm(d) < 1e-4;
        );

    -- Test 0
    R:= QQ[x,y,z];
    f := x^2 + y^2 +z^2;
    (p,q) := sosdecTernary (f, Solver=>solver);
    t0 := cmp(f,p,q);

    -- Test 1
    R = QQ[x,y,z];
    f = x^4*y^2 + x^2*y^4 + z^6 - 4*x^2 *y^2 * z^2;
    (p,q) = sosdecTernary (f, Solver=>solver);
    t1 := (p===null);

    -- Test 2
    R = RR[x,y,z];
    f = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2; --Motzkin
    (p,q) = sosdecTernary (f, Solver=>solver);
    t2 := cmp(f,p,q);

    results := {t0,t1,t2};
    informAboutTests (results);
    return results
    )


-- check sosInIdeal
checkSosInIdeal = solver -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local sol; local s; local mult;
    cmp := (h,s,mult) -> (
        if s===null then return false;
        h = sub(h,ring s);
        d := (h*mult)_(0,0) - sumSOS s;
        if isExactField h then return d==0;
        return norm(d)<1e-4;
        );

    -- Test 0
    R:= QQ[x];
    h:= matrix {{x+1}};
    sol = sosInIdeal (h,2, Solver=>solver);
    s = sosPoly sol;
    mult = sol#"mult";
    t0 := cmp(h,s,mult);
    
    -- Test 1 (similar to test 0)
    R= RR[x];
    h= matrix {{x+1}};
    sol = sosInIdeal (h,4, Solver=>solver);
    s = sosPoly sol;
    mult = sol#"mult";
    t1 := cmp(h,s,mult);

    -- Test 2:
    R = RR[x,y,z];
    h = matrix {{x-y, x+z}};
    sol = sosInIdeal (h,2, Solver=>solver);
    s = sosPoly sol;
    mult = sol#"mult";
    t2 := cmp(h,s,mult);

    -- Test 3: (similar to test 2)
    R = RR[x,y,z];
    h = matrix {{x-y, x+z}};
    sol = sosInIdeal (h,6, Solver=>solver);
    s = sosPoly sol;
    mult = sol#"mult";
    t3 := cmp(h,s,mult);

    results := {t0,t1,t2,t3};
    informAboutTests (results);
    return results
    )


-- check lowerBound
checkLowerBound = solver -> (
    tol := 0.001;
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    equal := (a,b) -> (
        if a===null then return false;
        d := if abs(b)<1 then abs(a-b) else abs(a-b)/abs(b);
        return d < tol;
        );
    cmp := (f,h,bound,mon,Q,mult) -> (
        if Q===null then return false;
        d := f - bound + (h*mult - transpose mon * Q * mon)_(0,0);
        if isExactField h then return d==0;
        return norm(d)<1e-4;
        );

    --------------UNCONSTRAINED1--------------
    --- Test 0
    R := QQ[x];
    f := (x-1)^2 + (x+3)^2;
    (bound,sol) := lowerBound(f, Solver=>solver);
    t0 := equal(bound,8);

    -- Test 1
    R = RR[x,y];
    f = (x-pi*y)^2 + x^2 + (y-4)^2;
    (bound,sol) = lowerBound(f, Solver=>solver);
    t1 := equal(bound,16*pi^2/(2+pi^2));

    -- Test 2
    R = QQ[x,z];
    f = x^4+x^2+z^6-3*x^2*z^2;
    (bound,sol) = lowerBound (f,Solver=>solver,RndTol=>infinity);
    t2 := equal(bound,-.17798);

    -- Test 3 (rational function)
    R = QQ[x];
    f = (x^2-x)/(x^2+1);
    (bound,sol) = lowerBound(f, Solver=>solver, RndTol=>infinity);
    t3 := equal(bound,1/2-1/sqrt(2));

    ---------------CONSTRAINED1---------------
    --- Test 4
    R = RR[x,y];
    f = y;
    h := matrix {{y-pi*x^2}};
    (bound,sol) = lowerBound (f, h, 4, Solver=>solver);
    (mon,Q,X,tval) := readSdpResult sol;
    mult := sol#"mult";
    t4 := equal(bound,0) and cmp(f,h,bound,mon,Q,mult);

    -- Test 5
    R = QQ[x,y,z];
    f = z;
    h = matrix {{x^2 + y^2 + z^2 - 1}};
    (bound,sol) = lowerBound (f, h, 4, Solver=>solver);
    (mon,Q,X,tval) = readSdpResult sol;
    mult = sol#"mult";
    t5 := equal(bound,-1) and cmp(f,h,bound,mon,Q,mult);

    -----------------QUOTIENT1-----------------
    -- Test 6
    R = QQ[x,y];
    I := ideal (x^2 - x);
    S := R/I;
    f = sub(x-y,S);
    h = matrix {{sub(y^2 - y,S)}};
    (bound,sol) = lowerBound(f, h, 2, Solver=>solver);
    (mon,Q,X,tval) = readSdpResult sol;
    mult = sol#"mult";
    t6 := equal(bound,-1) and cmp(f,h,bound,mon,Q,mult);
    
    results := {t0,t1,t2,t3,t4,t5,t6};
    informAboutTests (results);
    return results
    )

--##########################################################################--
-- Documentation and Tests
--##########################################################################--

beginDocumentation()

load "./SOS/SOSdoc.m2"

TEST /// --sosPoly and sumSOS
    R = QQ[x,y,z]
    coeff1={3,1,1,1/4,1}
    pol1={-(1/2)*x^3*y-(1/2)*x*y^3+x*y*z^2, x^2*y^2-z^4, x^2*y*z-y*z^3,
      -x^3*y+x*y^3, x*y^2*z-x*z^3}
    p1=sosPoly(R,pol1,coeff1)
    p2=x^6*y^2 + 2*x^4*y^4 + x^2*y^6 - 2*x^4*y^2*z^2 - 2*x^2*y^4*z^2 - 
    3*x^2*y^2*z^4 + x^2*z^6 + y^2*z^6 + z^8
    assert(sumSOS(p1)===p2)
///

TEST /// --SOSmult
    R = QQ[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1,2,3})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3,1/2,1/4})
    assert(sumSOS(p1*p2)==sumSOS(p1)*sumSOS(p2))
    assert(sumSOS(p1^4)==sumSOS(p1)^4)

    equal = (f1,f2) -> norm(f1-f2) < 1e-8;
    R = RR[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1.32,1.47,12./7})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3.1,1.31,2.0})
    assert( equal(sumSOS(p1*p2),sumSOS(p1)*sumSOS(p2)) )
    assert( equal(sumSOS(p1^4),sumSOS(p1)^4) )
///

TEST /// --cleanSOS
    R = RR[x,y];
    s = sosPoly(R, {x^2+.0001*x+1,y}, {2,.0001})
    t1 = clean( .001, s )
    t2 = sosPoly(R, {x^2+1}, {2})
    assert (t1 == t2)
    
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,1/100000})
    t = clean( 0.001, s )
    assert (t == s)
///

TEST ///--substitute SOSPoly
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3})
    S = QQ[x,y,z]
    t1 = sosPoly(S, {x+1,y}, {2,3})
    t2 = sub (s, S)
    assert (t1 == t2)
///

TEST ///--toRing
    debug needsPackage "SOS"
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3});
    S = RR[x,y];
    s2 = toRing_S s;
    assert instance(coefficientRing ring s2, RealField)
    s3 = toRing_R s2;
    assert (s==s3)
    
    tol := 1e-10;
    f = 0.1*x_S^2 + y^2
    g = 1/10*(symbol x)_R^2 + (symbol y)_R^2
    -- comparison in rationals is complicated:
    residual = sum \\ abs \ (x -> lift (x,QQ)) \ flatten entries last coefficients (toRing_R f - g)
    assert (residual < tol)
    -- comparison in reals:
    assert (norm (toRing_S g - f) < tol)
///

TEST /// --sosdec
    R=QQ[x,y,z]
    Q=matrix{{1,-1/2,1},{-1/2,1,-1/2},{1,-1/2,1}}
    Q=promote(Q,QQ)
    mon=matrix{{x^3},{x^2*z},{y*z^2}}
    f=sosPoly(mon,Q)
    assert(f=!=null and sumSOS f==transpose mon * Q *mon)
///

TEST /// --choosemonp
    debug needsPackage "SOS"
    R = QQ[x,y];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f)
    assert( lmsos === null )
    R = QQ[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f-t)
    assert( ring lmsos===R and numRows lmsos == 6 )
    
    R = RR[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = choosemonp(f-t)
    assert( ring lmsos===R and numRows lmsos == 6 )
///

TEST /// --createSOSModel
    debug needsPackage "SOS"
    eval = (Q,v) -> (transpose v * Q * v)_(0,0)
    
    R = QQ[x][t];
    f = x^4 - 2*x + t;
    mon = matrix{{1},{x},{x^2}}
    (C,Ai,p0,V,A,B,b) = createSOSModel(f,mon)
    assert( eval(C,mon) == x^4 - 2*x )
    assert( #Ai==2 and all({0,1}, j -> eval(Ai#j,mon) == V_(0,j)) )
    
    equal = (f1,f2) -> norm(f1-f2) < 1e-8;
    R = RR[x][t];
    f = x^4 - 2*x + t;
    mon = matrix{{1},{x},{x^2}}
    (C,Ai,p0,V,A,B,b) = createSOSModel(f,mon)
    assert( equal(eval(C,mon), x^4 - 2*x) )
    assert( #Ai==2 and all({0,1}, j -> equal(eval(Ai#j,mon), V_(0,j))) )
    
    R = QQ[x,y][t]
    f = x^2+2*x + t*(y+1)
    mon = matrix{{1},{x}}
    (C,Ai,p0,V,A,B,b) = createSOSModel(f,mon) --infeasible
    assert(entries C=={{0,1},{1,1}} and #Ai==0 and p0==0)
///

TEST /// --LDLdecomposition
    A = matrix(QQ, {{5,3,5},{3,2,4},{5,4,10}})
    (L,D,P,err) = LDLdecomposition A
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    (L,D,P,err) = LDLdecomposition promote(A,RR)
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    
    V = random(QQ^12,QQ^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and L*D*transpose L == transpose P * A * P)

    equal = (f1,f2) -> norm(f1-f2) < 1e-6;
    V = random(RR^12,RR^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and equal(L*D*transpose L, transpose P * A * P))
///

TEST /// --roundPSDmatrix
    Q=matrix{{2.01,0,0},{0,1.1,0},{0,0,2}}
    A=matrix{{1,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,1}}
    b=matrix{{2},{1},{2}}
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,10)
    Qtrue = matrix(QQ,{{2,0,0},{0,1,0},{0,0,2}})
    assert(Qpsd==Qtrue and boolv)

    -- BUG: the next example fails with for d=1
    Q=matrix{{1,-0.75,-0.75},{-0.75,1,0.99},{-0.75,0.99,1}}
    A=matrix{{1,0,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,1}}
    b=matrix{{1},{1},{1},{1}}
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,10)
    Qtrue = matrix{{1, -3/4, 1}, {-3/4, 1, 507/512}, {1, 507/512, 1}};
    assert(Qpsd==Qtrue and not boolv)

    Q=matrix{{2,-1,-0.75},{-1,2,-1.1},{-0.75,-1.1,2}}
    A=matrix {{1, 3, 0, 0, 0, 0}, {2, 0, 5, 1, 0, 0}, {0, 0, 0, 4, -3, 1}, {6, 0, 1, 0, 0, -5}}
    b=matrix{{1},{4},{-3},{1}}
    (boolv,Qpsd)=roundPSDmatrix(Q,A,b,100)
    Qtrue = matrix{{27446399799074697971/15073547952809050112, -4124283948755215953/ 15073547952809050112, 8072814922052298793/45220643858427150336}, {-4124283948755215953/15073547952809050112, -24159897971001080447/ 45220643858427150336, 14488868131185674623/15073547952809050112}, {8072814922052298793/45220643858427150336, 14488868131185674623/ 15073547952809050112, 91377473489393942387/45220643858427150336}}
    assert(Qpsd==Qtrue and not boolv)
///

TEST ///--makeMultiples
    debug needsPackage "SOS"
    R = QQ[x,y,z]
    f1 = x + y
    f2 = x^2 + z^2
    h = {f1,f2}
    (H,m) = makeMultiples (h,3, false)
    assert (#H == 14 )
    assert( max(first\degree\H) == 3 )
    
    (H,m) = makeMultiples (h,3, true)
    assert(#H == 9)
    assert( unique(first\degree\H) == {3} )
///

TEST ///--recoverSolution
    R = RR[x,y];
    mon = matrix {{1},{x},{y}};
    X = matrix(RR, {{1,0,1},{0,0,0},{1,0,1}} );
    sol = recoverSolution(mon,X);
    assert(sol#0#1==0 and sol#1#1==1)
///

TEST /// --solveSDP
    test := checkSolver("M2","solveSDP")
    assert(test#0 and test#1 and test#2) --(test3 fails)
///

TEST /// --solveSOS
    test := checkSolver("M2","solveSOS")
    assert(all(test,identity))
///
