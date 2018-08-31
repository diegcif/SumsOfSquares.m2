newPackage(
    "SOS",
    Version => "2.0",
    Date => "August 22, 2018",
    Authors => {
     {Name => "Diego Cifuentes",
      Email => "diegcif@mit.edu",
      HomePage => "http://www.mit.edu/~diegcif/"},
     {Name => "Thomas Kahle",
      Email => "thomas.kahle@ovgu.de",
      HomePage => "https://thomas-kahle.de/"},
     {Name => "Pablo A. Parrilo",
      Email => "parrilo@mit.edu",
      HomePage => "http://www.mit.edu/~parrilo/"},
     {Name => "Helfried Peyrl", 
      Email => "peyrl@control.ee.ethz.ch",
      HomePage => "https://scholar.google.com/citations?user=cFOV7nYAAAAJ&hl=de"},
     {Name => "Special thanks: Ilir Dema, Nidhi Kaihnsa, Anton Leykin"}
    },
    Headline => "Sum-of-Squares Package",
    Configuration => {"CSDPexec"=>"csdp","MOSEKexec"=>"mosek","SDPAexec"=>"sdpa","DefaultSolver"=>null},
    AuxiliaryFiles => true,
    UseCachedExampleOutput => false,
    CacheExampleOutput => true,
--  DebuggingMode => true,
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
    "library",
--Method options
    "GramMatrix",
    "MomentMatrix",
    "Parameters",
    "RoundTol",
    "Solver",
    "TraceObj",
    "Scaling"
}

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

-- Solver executables
makeGlobalPath = (fname) -> (
    -- Turns a file name into a global path of the file
    -- Used to find global file names of external solvers
    tmp := temporaryFileName();
    r := run( "which '" | fname | "' > " | tmp);
    if r>0 then return;
    fname = replace("\n","",get tmp);
    if first fname != "/" then fname = currentDirectory() | fname;
    return "'" | fname | "'";
    )

-- Choose default solver
chooseDefaultSolver = execs -> (
    solvers := {"CSDP", "MOSEK", "SDPA"}; --sorted by preference
    found := for i to #solvers-1 list
        if execs#i=!=null then solvers#i else continue;
    print if #found>0 then "Solvers found: "|demark(", ",found)
        else "Warning: No external solver was found.";
    found = append(found,"M2");
    defaultSolver = ((options SOS).Configuration)#"DefaultSolver";
    if not member(defaultSolver,found) then
        defaultSolver = first found;
    print("Default solver: " | defaultSolver);
    return defaultSolver;
    )

csdpexec = makeGlobalPath ((options SOS).Configuration)#"CSDPexec"
mosekexec = makeGlobalPath ((options SOS).Configuration)#"MOSEKexec"
sdpaexec = makeGlobalPath ((options SOS).Configuration)#"SDPAexec"
defaultSolver = chooseDefaultSolver(csdpexec,mosekexec,sdpaexec)

-- SDP status
StatusFeas = "Status: SDP solved, primal-dual feasible"
StatusPFeas = "Status: SDP solved, primal feasible"
StatusDFeas = "Status: SDP solved, dual feasible"
StatusPInfeas = "Status: SDP solved, primal infeasible"
StatusDInfeas = "Status: SDP solved, dual infeasible"
StatusFailed = "Status: SDP failed"

-- Constants
MaxRoundTol = 32 --maximum rounding tolerance
HighPrecision = 1e-10 --e.g. for numerical linear algebra
MedPrecision = 1e-6 --e.g. for SDP solutions
LowPrecision = 1e-4

--##########################################################################--
-- TYPES
--##########################################################################--

SOSPoly = new Type of HashTable

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

SOSPoly == Matrix := (S, F) -> (
    if numRows F!=1 or numColumns F!=1 then
        error "matrices have different shapes"
    else S == F_(0,0)
    )

Matrix == SOSPoly := (F, S) -> S == F

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

sdpResult = (mon,Q,X,tval) -> (
    new SDPResult from {
        Monomials => mon,
        GramMatrix => Q,
        MomentMatrix => X,
        Parameters => tval,
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
    tval := sol#Parameters;
    if tval=!=null and numRows tval>0 then
        str = append(str,{"Parameters",mat2str tval});
    return netList(str,HorizontalSpace=>1,Alignment=>Center)
    )

-- Shortcut to extract keys from SDPResult:
readSdpResult = sol -> (sol#Monomials, sol#GramMatrix, sol#MomentMatrix, sol#Parameters)


--##########################################################################--
-- METHODS
--##########################################################################--

verbose = (s,o) -> if o.Verbose then print s

-- library of nonnegative polynomials
library = method()
library(String,List) := (name,X) -> (
    if #X<3 then error "Insufficient variables";
    x := X_0; y := X_1; z := X_2;
    if name=="Motzkin" then
        return  x^4*y^2 + x^2*y^4 - 3*x^2*y^2*z^2 + z^6;
    if name=="Robinson" then
        return x^6 + y^6 + z^6 - (x^4*y^2 + x^2*y^4 + x^4*z^2 + x^2*z^4 + y^4*z^2 + y^2*z^4) + 3*x^2*y^2*z^2;
    if name=="Schmuedgen" then
        return 200*(x^3 - 4*x*z^2)^2 + 200*(y^3 - 4*y*z^2)^2 + 
           (y^2 - x^2)*x*(x + 2*z)*(x^2 - 2*x*z + 2*y^2 - 8*z^2);
    if name=="Scheiderer" then
        return x^4 + x*y^3 + y^4 - 3*x^2*y*z - 4*x*y^2*z + 2*x^2*z^2 + x*z^3 + y*z^3 + z^4;
    if name=="Harris" then(
        (a,b,c,d,e) := (16,-36,20,57,-38);
        return a*( x^10 + y^10 + z^10)+ 
            b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
            c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) + 
            d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
            e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4);
        );
    if #X<4 then error "Insufficient variables";
    w := X_3;
    if name=="Lax-Lax" then
        return (x-y)*(x-z)*(x-w)*x+(y-x)*(y-z)*(y-w)*y+(z-x)*(z-y)*(z-w)*z+(w-x)*(w-y)*(w-z)*w+x*y*z*w;
    if name=="Choi-Lam" then
        return x^2*y^2 + y^2*z^2 + x^2*z^2 + w^4 - 4*x*y*z*w;
    error "Name was not recognized.";
    )
library(String,Ring) := (name,R) -> library(name,gens R)

--###################################
-- Transition QQ <=> RR
--###################################

isExactField = kk -> (
    try (kk = ring kk);
    kk = ultimate(coefficientRing,kk);
    return precision 1_kk == infinity;
    )

isZero = (tol,x) -> if isExactField x then x==0 else norm x<tol

-- rounds real number to rational
roundQQ = method()
roundQQ(ZZ,RR) := (d,x) -> round(x*2^d)/2^d
roundQQ(ZZ,Matrix) := (d,X) -> 
     matrix(QQ, applyTable (entries X, roundQQ_d ));

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
    prec := precision kk;
    coef = matrix(QQ, {for c in flatten entries coef list roundQQ(prec,sub(c,RR))});
    f' := (mon*transpose coef)_(0,0);
    return f';
    )

toRing (Ring, SOSPoly) := (S, s) -> (
    -- maps s to ring S
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
    prec := precision kk;
    c' := for c in coefficients s list roundQQ(prec,sub(c,RR));
    return sosPoly (S, g', c')
    )

liftMonomial = (S,f) -> (
    -- maps monomial f to ring S
    n := numgens S;
    e := first exponents f;
    e = e_(toList(0..n-1)); -- ignore some variables
    return S_e;
    )


--###################################
-- basicLinearAlgebra
--###################################

linsolve = (A,b) -> (
    -- This function becomes obsolete when solve
    -- has a threshold for infeasibility.
    if isExactField A then return try solve(A,b);
    tol := HighPrecision;
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
    -- Kernel up to a precision.  Becomes obselete when M2 gives the
    -- kernel of a numerical matrix.
    if isExactField A then return gens kernel A;
    tol := HighPrecision;
    (S,U,Vt) := truncatedSVD(A,-tol);
    return transpose Vt;
    )

zeros = (kk,m,n) -> map(kk^m,kk^n,{})

-- Store a symmetric matrix in a vector, avoiding duplication
smat2vec = method( Options => {Scaling => 1} )
smat2vec(List) := o -> A -> (
    n := #A;
    v := for i to n-1 list
        for j from i to n-1 list 
            if i==j then A#i#j else o.Scaling*A#i#j;
    return flatten v;
    )
smat2vec(Matrix) := o -> A -> matrix(ring A, apply(smat2vec(entries A,o), a->{a}))

-- Reverse of the above
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

PSDdecomposition = A -> (
    -- Factors a PSD matrix A = L D L^T
    -- with D diagonal
    kk := ring A;
    if isExactField kk then
        return LDLdecomposition(A);
    if kk=!=RR and not instance(kk,RealField) then
        error "field must be QQ or RR";
    tol := HighPrecision;
    (e,V) := eigenvectors(A,Hermitian=>true);
    err := if all(e, i -> i > -tol) then 0 else 1;
    e = max_0 \ e;
    D := diagonalMatrix e;
    P := id_(kk^(numRows A));
    return (V,D,P,err);
    )
    
LDLdecomposition = (A) -> (
    -- This implements Algorithm 4.2.2 from [Golub-VanLoan-2012]
    kk := ring A;
    if kk=!=QQ and kk=!=RR and not instance(kk,RealField) then
       error "field must be QQ or RR";
    tol := if isExactField kk then 0 else HighPrecision;

    n := numRows A;
    Ah := new MutableHashTable;
    for i to n-1 do for j to n-1 do Ah#(i,j) = A_(i,j);
    v := new MutableList from for i to n-1 list 0_kk;
    piv := new MutableList from toList(0..n-1);
    err := 0;

    permuteMat := (k,q) -> (  -- k<=q
        if k==q then return;
        tmp := piv#q; piv#q = piv#k; piv#k = tmp;
        for i to n-1 do (tmp := Ah#(i,q); Ah#(i,q) = Ah#(i,k); Ah#(i,k) = tmp;);
        for i to n-1 do (tmp := Ah#(q,i); Ah#(q,i) = Ah#(k,i); Ah#(k,i) = tmp;);
        );

    for k from 0 to n-1 do (
        q := k + maxPosition apply(k..n-1, i->Ah#(i,i));
        permuteMat(k,q);

        --  positive semidefinite?
        a := Ah#(k,k);
        if a < -tol then (err = k+1; break;);
        if a <= tol then(
            if any(k+1..n-1, i->abs(Ah#(i,k))>tol) then (
                err = k+1; break;);
            continue;
            );

        -- Schur complement
        for i from k+1 to n-1 do 
            v#i = Ah#(i,k);
        for i from k+1 to n-1 do(
            Ah#(i,k) = v#i/a;
            for j from k+1 to n-1 do
                Ah#(i,j) = Ah#(i,j) - (v#i*v#j)/a;
            );
    );

    L := map(kk^n,kk^n,(i,j)-> if i>j then Ah#(i,j) else if i==j then 1_kk else 0_kk);
    D := map(kk^n,kk^n,(i,j)->if i==j then Ah#(i,j) else 0_kk);
    P := submatrix(id_(kk^n),toList piv);

    return (L,D,P,err);
)

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
     Options => {RoundTol => 3, Solver=>null, Verbose => false, TraceObj => false} )
 
rawSolveSOS(Matrix,Matrix,Matrix) := o -> (F,objP,mon) -> (
    -- Consider a parametric problem f = f_0 + f_1 p_1 + ... + f_s p_s
    -- This minimizes a function over the p_i
    -- F is a column vector with the f_i
    -- objP is a column vector specifying the objective function
    -- mon is a vector of monomials
    -- (see parameterVector)

    kk := coefficientRing ring F;

    -- checkInputs --
    if numColumns mon > 1 then error("Monomial vector should be a column.");
    isMonomial := max(length \ terms \ flatten entries mon)==1;
    if not isMonomial then error("Vector must consist of monomials.");
         
    -- build SOS model --
    (C,Ai,p0,V,A,B,b) := createSOSModel(F,mon,Verbose=>o.Verbose);
    if C===null then return sdpResult(mon,,,);

    ndim := numRows C;
    np := numRows objP;

    obj := 
        if o.TraceObj then
            map(RR^(#Ai),RR^1,(i,j)-> trace Ai#i)
        else
            (transpose V * objP); 
    if obj==0 then verbose( "Solving SOS feasibility problem...", o)
    else verbose("Solving SOS optimization problem...", o);

    (X,my,Q) := solveSDP(C, Ai, obj, Solver=>o.Solver, Verbose=>o.Verbose);
    if Q===null then return sdpResult(mon,Q,X,);
    y := -my;
    pVec0 := (p0 + V*y);

    if not isExactField kk then return sdpResult(mon,Q,X,pVec0);
    if o.RoundTol > MaxRoundTol then(
        if o.RoundTol < infinity then
            print "Warning: RoundTol is too high. Rounding will be skipped.";
        return sdpResult(changeMatrixField(RR,mon),Q,X,pVec0);
        );

    -- rational rounding --
    (ok,Qp,pVec) := roundSolution(pVec0,Q,A,B,b,RoundTol=>o.RoundTol,Verbose=>o.Verbose);
    if ok then return sdpResult(mon,Qp,X,pVec);
    print "rounding failed, returning real solution";
    return sdpResult(changeMatrixField(RR,mon),Q,X,pVec0);
    )

-- Choose monomials internally:
rawSolveSOS(Matrix,Matrix) := o -> (F,objP) -> (
    mon := chooseMons (F,Verbose=>o.Verbose);
    if mon===null then return (,,,);
    return rawSolveSOS(F,objP,mon,o);
    )
rawSolveSOS(Matrix) := o -> (F) -> 
    rawSolveSOS(F,zeros(QQ,numRows F-1,1),o)

solveSOS = method(
     Options => {RoundTol => 3, Solver=>null, Verbose => false, TraceObj => false} )

solveSOS(RingElement,RingElement,Matrix) := o -> (f,objFcn,mon) -> (
    (F,objP) := parameterVector(f,objFcn);
    return rawSolveSOS(F,objP,mon,o);
    )
solveSOS(RingElement,Matrix) := o -> (f,mon) -> 
    solveSOS(f,0_(ring f),mon,o)

solveSOS(RingElement,RingElement) := o -> (f,objFcn) -> (
    (F,objP) := parameterVector(f,objFcn);
    mon := chooseMons (F,Verbose=>o.Verbose);
    if mon===null then return (,mon,,);
    return rawSolveSOS(F,objP,mon,o);
    )
solveSOS(RingElement) := o -> (f) -> 
    solveSOS(f,0_(ring f),o)

solveSOS(RingElement,RingElement,ZZ) := o -> (f,objFcn,D) -> (
    (F,objP) := parameterVector(f,objFcn);
    mon := chooseMons(F,D);
    return solveSOS(f,objFcn,mon,o);
    )
solveSOS(RingElement,ZZ) := o -> (f,D) -> 
    solveSOS(f,0_(ring f),D,o)

-- Main method to setup an SOS problem as an SDP problem
-- It is not exported, but there is (hidden) documentation for it.
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
    if x===null then(
        print "No Gram matrix exists. Terminate.";
        return (,,,) );
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

-- Choose monomial basis based on Newton polytope
chooseMons = method(
    Options => {Verbose => false} )
chooseMons(RingElement) := o -> (f) -> (
    F := parameterVector(f);
    mon := chooseMons(F);
    if mon===null then return;
    return sub(mon,ring f);
    )
chooseMons(Matrix) := o -> (F) -> (
    R := ring F;
    if F==0 then error "Expected nonzero inputs.";
    if isQuotientRing R then error("A monomial vector or degree bound must be provided in quotient rings.");
    n:= numgens R;
    monsPoly := g -> set first entries monomials g;
    monsList := G -> if #G>0 then sum(monsPoly\G) else {};
    filterVerts := (verts) -> (
        -- only consider those without parameters (this is a hack!)
        lmpars := monsList drop(flatten entries F,1);
        return select(verts, v -> not member(R_v,lmpars));
        );
    lmf := monsList flatten entries F;
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

-- Choose monomials, given a degree bound
chooseMons(Matrix,ZZ) := o -> (F,D) -> (
    if D<=0 or odd D then error "Expected even positive integer";
    R := ring F;
    mon := if isHomogeneous R and isHomogeneous F then basis(D//2,R)
        else basis(0,D//2,R);
    verbose("#monomials: " | numColumns mon, o);
    return transpose mon;
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

--###################################
-- Rational Rounding
--###################################

roundSolution = {RoundTol=>3,Verbose=>false} >> o -> (pVec0,Q,A,B,b) -> (
    -- round and project --
    d := o.RoundTol;
    np := numColumns B;
    pVec := null;
    
    print "Start rational rounding";
    while (d < MaxRoundTol) do (
        verbose("rounding step #" | d, o);
        if np!=0 then (
            pVec = roundQQ(d,pVec0);
            bPar := b - B*pVec;
            ) 
        else bPar= b;

        (ok,Qp) := roundPSDmatrix(Q,A,bPar,d);
        if ok then break else d = d + 1;
        );
    return (ok,Qp,pVec);
    )

project2linspace = (A,b,x0) -> (
     -- cast into QQ (necessary class to compute inverse)
     A2 := promote (A,QQ);
     -- convert b into a matrix if it is a scalar in QQ/ZZ:
     b2 := promote (matrix{{b}},QQ);
     x02 := promote (x0,QQ);

     -- compute projection:
     xp := x02 - transpose(A2)*((A2*x02-b2)//(A2*transpose(A2)))
     )

roundPSDmatrix = (Q,A,b,d) -> (
     Q0 := roundQQ(d,Q);
     x0 := smat2vec(Q0);
     xp := project2linspace(A,b,x0);
     Q = vec2smat(xp);

     (L,D,P,Qpsd) := LDLdecomposition(Q);
     if Qpsd == 0 then (true, Q) else (false,Q)
     )

--###################################
-- SOS in IDEAL
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
     Options => {RoundTol => 3, Solver=>"CSDP", Verbose => false} )
sosInIdeal (Ring, ZZ) := o -> (R,D) -> (
    -- find sos polynomial in a quotient ring
    if odd D then error "D must be even";
    mon := chooseMons(matrix{{0_R}}, D);
    sol := solveSOS (0_R, mon, o);
    (mon',Q,X,tval) := readSdpResult sol;
    if Q===null or isZero(MedPrecision,Q) then (
        print("no sos polynomial in degree "|D);
        return sdpResult(mon,,X,);
        );
    return sol;
    )
sosInIdeal (Matrix,ZZ) := o -> (h,D) -> (
    -- h is a row vector of polynomials
    -- D is a degree bound
    -- returns sos polynomial in <h>

    -- The output is an SDPResult and the multipliers that
    -- express the SOS in terms of the generators.
    -- We have them independent of Gröbner bases.
    
    if numRows h > 1 then error "h must be a row vector";
    if odd D then error "D must be even";
    homog := isHomogeneous h;
    (H,m) := makeMultiples(first entries h, D, homog);
    F := matrix transpose {{0}|H};
    sol := rawSolveSOS (F, o);
    (mon,Q,X,tval) := readSdpResult sol;
    if Q===null or isZero(MedPrecision,Q) then (
        print("no sos polynomial in degree "|D);
        return (sdpResult(mon,,X,),null);
        );
    tval = flatten entries tval;
    mult := getMultipliers(m,tval,ring mon);
    return (sol,mult);
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
     Options => {RoundTol => 3, Solver=>"CSDP", Verbose => false} )
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
        (sol,mult) := sosInIdeal(matrix{{fi}},2*di-4,o);
        Si := sosPoly sol;
        if Si===null then return (,);
        fi = mult_(0,0);
        if fi==0 then return (,);
        di = first degree fi;
        S = append(S,Si);
        );
    (mon,Q,X,tval) := readSdpResult rawSolveSOS matrix{{fi}};
    if Q===null or isZero(MedPrecision,Q) then return (,);
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

recoverSolution = method()
recoverSolution(Matrix,Matrix) := (mon,X) -> (
    if X===null then return {};
    e := eigenvalues(X,Hermitian=>true);
    if e#(-1)<=0 or e#0/e#(-1) < -HighPrecision then 
        error "Moment matrix is not positive semidefinite";
    i0 := position(flatten entries mon, i -> i==1);
    if i0===null then
        error "The monomial vector must contain 1";
    if e#(-2) > LowPrecision then 
        print "Moment matrix is not rank one, solution might not be correct.";
    sol := for i to numRows mon-1 list (
        y := mon_(i,i0);
        if sum degree y!=1 then continue;
        y => X_(i,i0) );
    return sol;
    )
recoverSolution(SDPResult) := sol -> recoverSolution(sol#Monomials,sol#MomentMatrix)

-- Unconstrained minimization 
-- sos lower bound for the polynomial f
lowerBound = method(
     Options => {RoundTol => 3, Solver=>null, Verbose => false} )
lowerBound(RingElement) := o -> (f) -> lowerBound(f,-1,o)
lowerBound(RingElement,ZZ) := o -> (f,D) -> drop(lowerBound(f,zeros(ring f,1,0),D,o),-1)

-- Minimize a polynomial on an algebraic set
lowerBound(RingElement,Matrix,ZZ) := o -> (f,h,D) -> (
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
        {RoundTol=>o.RoundTol, Solver=>o.Solver, Verbose=>o.Verbose};
    mon := if isQuotientRing R then chooseMons(F,D)
        else chooseMons (F,Verbose=>o.Verbose);
    if mon===null then return (,sdpResult(,,,),);
    sol := rawSolveSOS(F,objP,mon,o');
    (mon',Q,X,tval) := readSdpResult sol;
    (bound,mult) := (,);
    if tval=!=null then(
        tval = flatten entries tval;
        bound = tval#0;
        mult = getMultipliers(m,drop(tval,1),ring mon');
        );
    return (bound,sol,mult);
    )

--###################################
-- SOLVE SDP
--###################################

solveSDP = method(
     Options => {Solver=>null, Verbose => false} )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C,sequence A,b,o)

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C,sequence A,b,y,o)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    if numRows b=!=#A then error "Bad matrix dimensions.";
    (C,A,b) = toReal(C,A,b);
    (ok,X,y,Z) := (,,,);
    (ok,X,y,Z) = sdpNoConstraints(C,A);
    if ok then return (X,y,Z);
    solver := chooseSolver o;
    if solver == "M2" then(
        (ok,X,y,Z) = trivialSDP(C,A,b);
        if ok then return (X,y,Z)
        else (y,Z) = simpleSDP(C,A,b,Verbose=>o.Verbose)
        )
    else if solver == "CSDP" then
        (X,y,Z) = solveCSDP(C,A,b,Verbose=>o.Verbose)
    else if solver == "SDPA" then
        (X,y,Z) = solveSDPA(C,A,b,Verbose=>o.Verbose)
    else if solver == "MOSEK" then
        (X,y,Z) = solveMOSEK(C,A,b,Verbose=>o.Verbose)
    else
        error "unknown SDP solver";
    ntries := 6;
    (y,Z) = findNonZeroSolution(C,A,b,o,y,Z,ntries);
    return (X,y,Z);
    )

findNonZeroSolution = (C,A,b,o,y,Z,ntries) -> (
    -- Heuristic to trick the solver into returning a nonzero solution of an SDP
    -- by changing the objective.
    if y===null then return (y,Z);
    if not(C==0 and b==0 and y==0) then return (y,Z);
    print "Zero solution obtained. Trying again.";
    m := numRows b;
    badCoords := set();
    iszero := a -> norm a < MedPrecision;
    for i to ntries-1 do(
        if #badCoords==m then break;
        b' := map(RR^m,RR^1, (j,l) -> 
            if member(j,badCoords) then 0 else random(RR)-.5 );
        (X',y',Z') := solveSDP(C,A,b',o);
        if Z'=!=null and not iszero Z' then return (y',Z');
        if X'===null and y'=!=null and (transpose(-b') * y')_(0,0) < -.1 then
            badCoords = badCoords + set select(0..m-1, j -> not iszero y'_(j,0));
        );
    return (y,Z);
    )

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y0) -> (
    (C,A,b) = toReal(C,A,b);
    y0 = promote(y0,RR);
    (ok,X,y,Z) := (,,,);
    (ok,X,y,Z) = sdpNoConstraints(C,A);
    if ok then return (X,y,Z);
    if chooseSolver o != "M2" then return solveSDP(C,A,b,o);
    (ok,X,y,Z) = trivialSDP(C,A,b);
    if ok then return (X,y,Z);
    (y,Z) = simpleSDP2(C,A,b,y0,false,Verbose=>o.Verbose);
    return (,y,Z);
    )

chooseSolver = o -> if o.Solver=!=null then o.Solver else defaultSolver

toReal = (C,A,b) -> (
    C = promote(C,RR);
    A = apply(A, Ai -> promote(Ai,RR));
    b = promote(b,RR);
    return (C,A,b);
    )

-- Solve very simple SDP with no constraints
sdpNoConstraints = (C,A) -> (
    tol := HighPrecision;
    if #A==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=-tol then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, 0*C, y0, C);
        )else(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
    )

-- check trivial cases and solve them directly
trivialSDP = (C,A,b) -> (
    if #A==0 or b==0 then(
        lambda := min eigenvalues(C, Hermitian=>true);
        if lambda>=0 then(
            print "SDP solved";
            y0 := zeros(RR,#A,1);
            return (true, 0*C, y0, C);
        )else if #A==0 then(
            print "dual infeasible";
            return (true,,,);
            );
        );
    return (false,,,);
    )


-- Implementation of SDP in Macaulay2
-- Algorithm: Dual interior point method
-- see Boyd, Vandenberghe "Convex Optimization" pp. 618-619, pp. 463-466
simpleSDP = {Verbose => false} >> o -> (C,A,b) -> (
    print "Running M2 Solver";
    R := RR;
    n := numRows C;

    -- try to find strictly feasible starting point --
    local y; local Z;
    lambda := min eigenvalues (C, Hermitian=>true);
    if lambda > 0 then
        y = zeros(R,#A,1)
    else(
        verbose("Computing strictly feasible solution...", o);
        y =  zeros(R,#A,1) || matrix{{lambda*1.1}};
        obj :=  zeros(R,#A,1) || matrix{{1_R}};
        (y,Z) = simpleSDP2(C,append(A,id_(R^n)), obj, y, true, Verbose=>o.Verbose);
        if y===null then (
            print StatusFailed;
            return (,) );
        y = transpose matrix {take (flatten entries y,numRows y - 1)};
        );
    verbose("Computing an optimal solution...", o);
    (y,Z) = simpleSDP2(C, A, b, y, false, o);
    print if y=!=null then StatusDFeas else StatusFailed;
    return (y,Z);
    )


-- This second part solves given an interior starting point.
simpleSDP2 = {Verbose => false} >> o -> (C,A,mb,y,UntilObjNegative) -> (
    R := RR;
    n := numgens target C;
    b := -mb;

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
                print "Slack matrix is singular";
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
            if UntilObjNegative and (obj_(0,0) < 0) then break;
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

--###################################
-- Interface to CSDP
--###################################

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
    print("Executing CSDP");
    print("Input file: " | fin);
    runcmd("cd " | dir | " && " | csdpexec | " " | fin1 | " " | fout | ">" | fout2);
    print("Output file: " | fout);
    (X,y,Z) := readCSDP(fout,fout2,n,o.Verbose);
    y = checkDualSol(C,A,y,Z,o.Verbose);
    return (X,y,Z);
    )

runcmd = (cmd) -> (
    tmp := getFileName ".err";
    r := run(cmd | " 2> " | tmp);
    if r == 32512 then error "Executable not found.";
    if r == 11 then error "Segmentation fault.";
    if r>0 then (
        txt := get tmp;
        if #txt>0 then(
            print txt;
            error "Command could not be executed." );
        );
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

-- SDPA file format is a shared input format of SDPA and CSDP
writeSDPA = (fin,C,A,b) -> (
    digits := 16;
    formatD := format_digits;
    m := length A;
    n := numColumns C;
    A = prepend(C,A);
    f := openOut fin;
    smat2str := (a,pref) -> (
        s := "";
        for i to n-1 do
            for j from i to n-1 do
                if a_(i,j)!=0 then
                    s = s | pref | i+1 | " " | j+1 | " " | formatD a_(i,j) | "\n";
        return s;
        );
    f << "*SDPA file generated by SOSm2" << endl;    
    f << m << " =mdim" << endl;
    f << "1 =nblocks" << endl;
    f << n << endl;
    f << demark(" ", formatD \ flatten entries(-b)) << endl;
    for i to m do
        f << smat2str(-A#i, i|" 1 ");
    f << close;
    )

-- Writes parameter file for CSDP.
-- All but one are also defaults
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
    f << "perturbobj=0" << endl; -- This one is changed from the default
    f << "fastmode=0" << endl;
    f << close;
    )

-- read CSDP output
readCSDP = (fout,fout2,n,Verbose) -> (
    sdpa2matrix := s -> (
        e := for i in s list (i_2-1,i_3-1) => i_4;
        e' := for i in s list (i_3-1,i_2-1) => i_4;
        return map(RR^n, RR^n, e|e');
        );
    readLine := l -> for s in separate(" ",l) list if s=="" then continue else value s;
    --READ SOLUTIONS
    text := get fout;
    text = replace("\\+","",text);
    L := lines text;
    y := matrix(RR,transpose{readLine L_0});
    S := readLine \ drop(L,1);
    S1 := select(S, l -> l_0==1);
    S2 := select(S, l -> l_0==2);
    Z := sdpa2matrix(S1); -- slack matrix
    X := sdpa2matrix(S2); -- dual solution
    -- READ STATUS
    text = get fout2;
    s := select(lines text, l -> match("Success",l));
    if #s==0 then( print StatusFailed; return (,,) );
    s = first s;
    if Verbose then print text;
    if match("SDP solved",s) then 
        print StatusFeas
    else if match("primal infeasible",s) then(
        print StatusPInfeas;
        X=null; )
    else if match("dual infeasible",s) then (
        print StatusDInfeas;
        y=null;Z=null; )
    else 
        print("Warning: Solver returns unknown message!!! " |s);
    return (X,y,Z);
)

-- A heuristic to postprocess output of CSDP
checkDualSol = (C,A,y,Z,Verbose) -> (
    if y===null then return;
    yA := sum for i to #A-1 list y_(i,0)*A_i;
    if norm(Z-C+yA)<MedPrecision then return y;
    if Verbose then print "updating dual solution";
    AA := transpose matrix(RR, smat2vec \ entries \ toList A);
    bb := transpose matrix(RR, {smat2vec entries(C-Z)});
    y = solve(AA,bb,ClosestFit=>true);
    return y;
    )

--###################################
-- Interface to SDPA
--###################################

solveSDPA = method( Options => {Verbose => false} )
solveSDPA(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    if sdpaexec===null then error "sdpa executable not found";
    n := numColumns C;
    fin := getFileName ".dat-s";
    fout := getFileName "";
    writeSDPA(fin,C,A,b);
    print("Executing SDPA");
    print("Input file: " | fin);
    runcmd(sdpaexec | " " | fin | " " | fout | "> /dev/null");
    print("Output file: " | fout);
    (X,y,Z) := readSDPA(fout,n,o.Verbose);
    return (X,y,Z);
    )

readSDPA = (fout,n,Verbose) -> (
    readVec := l -> (
        l = replace("([{} +])","",l);
        for s in separate(",",l) list if s=="" then continue else value s
    );
    readMatrix := ll -> 
        matrix(RR, for l in ll list readVec l);
    text := get fout;
    L := lines text;
    --READ SOLUTIONS
    y := null; X := null; Z := null;
    i := position(L, l -> match("xVec =",l));
    if i=!=null then 
        y = transpose matrix(RR, {readVec L#(i+1)});
    i = position(L, l -> match("xMat =",l));
    if i=!=null then 
        Z = matrix(RR, for j to n-1 list readVec L#(i+j+2));
    i = position(L, l -> match("yMat =",l));
    if i=!=null then 
        X = matrix(RR, for j to n-1 list readVec L#(i+j+2));
    --READ STATUS
    if Verbose then print text;
    s := first select(L, l -> match("phase.value",l));
    if match("pdOPT|pdFEAS",s) then 
        print StatusFeas
    else if match("dFEAS",s) then 
        print StatusPFeas
    else if match("pFEAS",s) then 
        print StatusDFeas
    else if match("dUNBD|pINF_dFEAS",s)  then(
        print StatusDInfeas;
        y=null;Z=null; )
    else if match("pUNBD|pFEAS_dINF",s) then(
        print StatusPInfeas;
        X=null; )
    else if match("noINFO|pdINF",s) then(
        print StatusFailed;
        X=null;y=null;Z=null; )
    else
        print("Warning: Solver returns unknown message!!! " |s);
    return (X,y,Z);
    )

--###################################
-- Interface to MOSEK
--###################################

solveMOSEK = method( Options => {Verbose => false} )
solveMOSEK(Matrix,Sequence,Matrix) := o -> (C,A,b) -> (
    if mosekexec===null then error "mosek executable not found";
    n := numColumns C;
    fin := getFileName ".cbf";
    fout := replace(".cbf",".sol",fin);
    fout2 := getFileName "";
    writeMOSEK(fin,C,A,b);
    print("Executing MOSEK");
    print("Input file: " | fin);
    runcmd(mosekexec | " " | fin | ">" | fout2);
    print("Output file: " | fout);
    (X,y,Z) := readMOSEK(fout,fout2,n,o.Verbose);
    return (X,y,Z);
    )

-- write mosek input file (CBF format)
writeMOSEK = (fin,C,A,b) -> (
    version := 1;
    digits := 16;
    formatD := format_digits;
    m := length A;
    n := numColumns C;
    f := openOut fin;
    smat2str := (a,pref) -> (
        s := "";
        for i to n-1 do
            for j from 0 to i do
                if a_(i,j)!=0 then
                    s = s | pref | i | " " | j | " " | formatD a_(i,j) | "\n";
        return s;
        );
    nlines := str -> #select("^",str)-1;
    -- header
    f << "# CBF file generated by SOSm2" << endl;    
    f << "VER" << endl << version << endl;
    f << "OBJSENSE" << endl << "MIN" << endl;
    f << "PSDVAR" << endl << 1 << endl << n << endl;
    f << "CON" << endl;
    f << m << " " << 1 << endl << "L= " << m << endl;
    -- objective function
    f << "OBJFCOORD" << endl;
    Cstr := smat2str(C, "0 ");
    f << nlines Cstr << endl << Cstr;
    -- constraints
    f << "FCOORD" << endl;
    Astr := concatenate for i to m-1 list 
        smat2str(A#i, i|" 0 ");
    f << nlines Astr << endl << Astr;
    -- constants
    f << "BCOORD" << endl << m << endl;
    for i to m-1 do
        f << i << " " << formatD(-b_(i,0)) << endl;
    f << close;
    )

readMOSEK = (fout,fout2,n,Verbose) -> (
    splitline := l -> separate(" ", replace(" +"," ",l));
    if Verbose then print get fout2;
    text := get fout;
    L := lines replace("\\+","",text);
    -- READ SOL y
    Ly := select(L, match_" EQ ");
    y := matrix(RR, for l in Ly list(
        l = splitline(l);
        {value l#6 - value l#7} ));
    -- READ SOL X,Z
    Lpsd := select(L, match_"BARX1");
    k := #Lpsd;
    Xh := new MutableList from 2*k:null;
    Zh := new MutableList from 2*k:null;
    for s to k-1 do(
        l := splitline(Lpsd#s);
        i := value l#2;
        j := value l#3;
        Xij := value l#4;
        Zij := value l#5;
        Xh#(2*s) = (i,j)=>Xij;
        Xh#(2*s+1) = (j,i)=>Xij;
        Zh#(2*s) = (i,j)=>Zij;
        Zh#(2*s+1) = (j,i)=>Zij;
        (i,j,Xij,Zij)
        );
    X := map(RR^n,RR^n,toList Xh);
    Z := map(RR^n,RR^n,toList Zh);
    -- READ STATUS
    s := select(L, l -> match("PROBLEM STATUS",l));
    if #s==0 then( print StatusFailed; return (,,) );
    s = first s;
    if match("PRIMAL_AND_DUAL_FEASIBLE",s) then 
        print StatusFeas
    else if match("PRIMAL_FEASIBLE",s) then
        print StatusPFeas
    else if match("DUAL_FEASIBLE",s) then
        print StatusDFeas
    else if match("UNKNOWN|ILL_POSED",s) then(
        print StatusFailed;
        X=null; y=null;Z=null; )
    else if match("PRIMAL_INFEASIBLE",s) then(
        print StatusPInfeas;
        X=null; )
    else if match("DUAL_INFEASIBLE",s) then (
        print StatusDInfeas;
        y=null;Z=null; )
    else 
        print("Warning: Solver returns unknown message!!! " |s);
    return (X,y,Z);
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
        return checkMethod#fun(solver,i->true);
    if fun != "AllMethods" then
        error "No test implemented for this function";
    T := for f in keys checkMethod list(
        print "################################";
        print("checking method "|f);
        print "################################";
        t := checkMethod#f(solver,i->true);
        informAboutTests t;
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

-- In the following methods the argument applyTest 
-- allows to specify which tests should be performed

--checkSolveSDP
checkSolveSDP = (solver,applyTest) -> (
    tol := .001;
    equal := (y0,y) -> y=!=null and norm(y0-y)<tol*(1+norm(y0));
    checkZ := (C,A,y,Z) -> if y===null then false
        else ( yA := sum for i to #A-1 list y_(i,0)*A_i; norm(Z-C+yA)<MedPrecision );
    local C; local b; local A; local A1; local A2; local A3; 
    local y0; local y; local X; local Z; local yopt;

    t0:= if applyTest(0) then(
        C = matrix{{0,2,0,0,0,0},{2,0,0,0,0,0},
         {0,0,10,0,0,0},{0,0,0,10,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
        A1 = matrix{{-1,0,0,0,0,0},{0,0,0,0,0,0},
         {0,0,1,0,0,0},{0,0,0,0,0,0},{0,0,0,0,-1,0},{0,0,0,0,0,0}};
        A2 = matrix{{0,0,0,0,0,0},{0,-1,0,0,0,0},
         {0,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1}};
        A = (A1,A2);
        y0 = matrix{{7},{9}};
        b = matrix{{-1},{-1}};
        (X,y,Z) = solveSDP(C,A,b,y0,Solver=>solver);
        yopt = matrix{{2.},{2.}};
        equal(yopt,y)
        );

    t1:= if applyTest(1) then(
        C = matrix {{2,1,-1},{1,0,0},{-1,0,5}};
        A1 = matrix {{0,0,1/2},{0,-1,0},{1/2,0,0}};
        A2 = matrix {{1,0,0},{0,1,0},{0,0,1}};
        A = (A1,A2);
        b = matrix {{0},{1}};
        y0 = matrix {{0},{-.486952}};
        (X,y,Z) = solveSDP(C,A,b,y0,Solver=>solver);
        yopt = matrix{{1.97619},{.466049}};
        equal(yopt,y)
        );

    t2:= if applyTest(2) then(
        C = matrix{{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
        A1 = matrix{{-1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
        A2 = matrix{{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
        A = (A1,A2);
        b = matrix{{-1},{0}};
        (X,y,Z) = solveSDP(C,A,b,Solver=>solver);
        yopt = matrix{{0.},{4.}};
        equal(yopt,y)
        );

    t3:= if applyTest(3) then( -- not strictly feasible
        C = matrix {{2,2,-1,3},{2,0,0,2},{-1,0,1,0},{3,2,0,1}};
        A1 = matrix {{0,0,0,1/2},{0,-1,0,0},{0,0,0,0},{1/2,0,0,0}};
        A = sequence A1;
        b = matrix {{-1}};
        (X,y,Z) = solveSDP(C,A,b,Solver=>solver);
        yopt = 4.;
        equal(yopt,y)
        );

    t4:= if applyTest(4) then( -- zero objective
        C = matrix(RR, {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
        A1 = matrix(RR, {{1, 3/2, 3/2}, {3/2, 0, 1/2}, {3/2, 1/2, 0}});
        A2 = matrix(RR, {{0, 1/2, 3/2}, {1/2, 0, 3/2}, {3/2, 3/2, 1}});
        A3 = matrix(RR, {{0, 0, 1/2}, {0, -1, 0}, {1/2, 0, 0}});
        A = (A1,A2,A3);
        b = matrix(RR, {{0}, {0}, {0}});
        (X,y,Z) = solveSDP(C, A, b, Solver=>solver);
        checkZ(C,A,y,Z)
        );

    return {t0,t1,t2,t3,t4};
    )

--checkSolveSOS
checkSolveSOS = (solver,applyTest) -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local w; w= symbol w;
    local t; t= symbol t;
    isGram := (f,mon,Q) -> (
        if Q===null then return false;
        e := eigenvalues(Q,Hermitian=>true);
        tol := MedPrecision;
        if min e < -tol then return false;
        S := ring mon;
        return isZero(tol, sub(f,S) - transpose(mon)*Q*mon);
        );
    isGramParam := (f,mon,Q,tval) ->
        if tval===null then false else isGram(sub(f,t=>tval_(0,0)),mon,Q);

    ---------------GOOD CASES1---------------
    t0:= if applyTest(0) then(
        R := QQ[x,y];
        f := 4*x^4+y^4;
        (mon,Q,X,tval) := readSdpResult solveSOS(f,Solver=>solver);
        isGram(f,mon,Q)
        );

    t1:= if applyTest(1) then(
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
        isGram(f,mon,Q)
        );

    t2:= if applyTest(2) then(
        R = QQ[x,y,z];
        f = x^4+y^4+z^4-4*x*y*z+x+y+z+3;
        (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
        isGram(f,mon,Q)
        );
    
    t3:= if applyTest(3) then(
        R = QQ[x,y,z,w];
        f = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
        (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver);
        isGram(f,mon,Q)
        );

    ---------------PARAMETRIC1---------------
    t4:= if applyTest(4) then(
        R = QQ[x][t];
        f = (t-1)*x^4+1/2*t*x+1;
        (mon,Q,X,tval) = readSdpResult solveSOS (f,Solver=>solver);
        isGramParam(f,mon,Q,tval)
        );

    ---------------QUOTIENT1---------------
    t5:= if applyTest(5) then(
        R = QQ[x,y];
        S := R/ideal(x^2 + y^2 - 1);
        f = sub(10-x^2-y,S);
        (mon,Q,X,tval) = readSdpResult solveSOS (f, 2, TraceObj=>true);
        isGram(f,mon,Q) and rank Q == 2
        );

    ---------------BAD CASES1---------------
    t6:= if applyTest(6) then(
        R = QQ[x,y][t];
        f = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1; --Motzkin
        (mon,Q,X,tval) = readSdpResult solveSOS(f,Solver=>solver); 
        ( Q === null )
        );

    t7:= if applyTest(7) then(
        (mon,Q,X,tval) = readSdpResult solveSOS(f-t,-t, Solver=>solver); 
        ( Q === null )
        );

    results := {t0,t1,t2,t3,t4,t5,t6,t7};
    return results;
    )

-- check sosdecTernary
checkSosdecTernary = (solver,applyTest) -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;

    cmp := (f,p,q) -> (
        if p===null then return false;
        d := product(sumSOS\p) - f*product(sumSOS\q);
        return isZero(LowPrecision, d);
        );

    t0:= if applyTest(0) then(
        R:= QQ[x,y,z];
        f := x^6 + y^6 +z^6;
        (p,q) := sosdecTernary (f, Solver=>solver);
        cmp(f,p,q)
        );

    t1:= if applyTest(1) then(
        R = QQ[x,y,z];
        f = x^4*y^2 + x^2*y^4 + z^6 - 4*x^2 *y^2 * z^2;
        (p,q) = sosdecTernary (f, Solver=>solver);
        (p===null)
        );

    t2:= if applyTest(2) then(
        R = RR[x,y,z];
        f = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2; --Motzkin
        (p,q) = sosdecTernary (f, Solver=>solver);
        cmp(f,p,q)
        );

    return {t0,t1,t2};
    )


-- check sosInIdeal
checkSosInIdeal = (solver,applyTest) -> (
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local sol; local s; local mult;
    cmp := (h,s,mult) -> (
        if s===null then return false;
        h = sub(h,ring s);
        d := (h*mult)_(0,0) - sumSOS s;
        return isZero(MedPrecision, d);
        );

    t0:= if applyTest(0) then(
        R:= QQ[x];
        h:= matrix {{x+1}};
        (sol,mult) = sosInIdeal (h,2, Solver=>solver);
        s = sosPoly sol;
        cmp(h,s,mult)
        );
    
    t1:= if applyTest(1) then( --similar to test0
        R= RR[x];
        h= matrix {{x+1}};
        (sol,mult) = sosInIdeal (h,4, Solver=>solver);
        s = sosPoly sol;
        cmp(h,s,mult)
        );

    t2:= if applyTest(2) then(
        R = RR[x,y,z];
        h = matrix {{x-y, x+z}};
        (sol,mult) = sosInIdeal (h,2, Solver=>solver);
        s = sosPoly sol;
        cmp(h,s,mult)
        );

    t3:= if applyTest(3) then( --similar to test 2
        R = RR[x,y,z];
        h = matrix {{x-y, x+z}};
        (sol,mult) = sosInIdeal (h,6, Solver=>solver);
        s = sosPoly sol;
        cmp(h,s,mult)
        );

    -----------------QUOTIENT1-----------------
    t4:= if applyTest(4) then(
        R = QQ[x,y,z]/ideal {x^2+y^2+y, y-z^2};
        s = sosPoly sosInIdeal (R,2,Solver=>solver);
        s=!=null and sumSOS s==0
        );
    
    return {t0,t1,t2,t3,t4};
    )


-- check lowerBound
checkLowerBound = (solver,applyTest) -> (
    tol := 0.001;
    local x; x= symbol x;
    local y; y= symbol y;
    local z; z= symbol z;
    local mult;
    equal := (a,b) -> (
        if a===null then return false;
        d := if abs(b)<1 then abs(a-b) else abs(a-b)/abs(b);
        return d < tol;
        );
    cmp := (f,h,bound,mon,Q,mult) -> (
        if Q===null then return false;
        d := f - bound + (h*mult - transpose mon * Q * mon)_(0,0);
        return isZero(MedPrecision, d);
        );

    --------------UNCONSTRAINED1--------------
    t0:= if applyTest(0) then(
        R := QQ[x];
        f := (x-1)^2 + (x+3)^2;
        (bound,sol) := lowerBound(f, Solver=>solver);
        equal(bound,8)
        );

    t1:= if applyTest(1) then(
        R = RR[x,y];
        f = (x-pi*y)^2 + x^2 + (y-4)^2;
        (bound,sol) = lowerBound(f, Solver=>solver);
        equal(bound,16*pi^2/(2+pi^2))
        );

    t2:= if applyTest(2) then(
        R = QQ[x,z];
        f = x^4+x^2+z^6-3*x^2*z^2;
        (bound,sol) = lowerBound (f,Solver=>solver,RoundTol=>infinity);
        equal(bound,-.17798)
        );

    t3:= if applyTest(3) then( --rational function
        R = QQ[x];
        f = (x^2-x)/(x^2+1);
        (bound,sol) = lowerBound(f, Solver=>solver, RoundTol=>infinity);
        equal(bound,1/2-1/sqrt(2))
        );

    ---------------CONSTRAINED1---------------
    t4:= if applyTest(4) then(
        R = RR[x,y];
        f = y;
        h := matrix {{y-pi*x^2}};
        (bound,sol,mult) = lowerBound (f, h, 2, Solver=>solver);
        (mon,Q,X,tval) := readSdpResult sol;
        equal(bound,0) and cmp(f,h,bound,mon,Q,mult)
        );

    t5:= if applyTest(5) then(
        R = QQ[x,y,z];
        f = z;
        h = matrix {{x^2 + y^2 + z^2 - 1}};
        (bound,sol,mult) = lowerBound (f, h, 4, Solver=>solver);
        (mon,Q,X,tval) = readSdpResult sol;
        equal(bound,-1) and cmp(f,h,bound,mon,Q,mult)
        );

    -----------------QUOTIENT1-----------------
    t6:= if applyTest(6) then(
        R = QQ[x,y];
        I := ideal (x^2 - x);
        S := R/I;
        f = sub(x-y,S);
        h = matrix {{sub(y^2 - y,S)}};
        (bound,sol,mult) = lowerBound(f, h, 2, Solver=>solver);
        (mon,Q,X,tval) = readSdpResult sol;
        equal(bound,-1) and cmp(f,h,bound,mon,Q,mult)
        );
    
    return {t0,t1,t2,t3,t4,t5,t6};
    )

--##########################################################################--
-- Documentation and Tests
--##########################################################################--

beginDocumentation()

load "./SOS/SOSdoc.m2"

--0
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

--1
TEST /// --SOSmult
    debug needsPackage "SOS"
    R = QQ[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1,2,3})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3,1/2,1/4})
    assert(sumSOS(p1*p2)==sumSOS(p1)*sumSOS(p2))
    assert(sumSOS(p1^4)==sumSOS(p1)^4)

    equal = (f1,f2) -> norm(f1-f2) < HighPrecision;
    R = RR[x,y,z,w]
    p1=sosPoly(R,{x^2-x*y,y^2+1,x},{1.32,1.47,12./7})
    p2=sosPoly(R,{y^3,x*w*z,y*z^2},{3.1,1.31,2.0})
    assert( equal(sumSOS(p1*p2),sumSOS(p1)*sumSOS(p2)) )
    assert( equal(sumSOS(p1^4),sumSOS(p1)^4) )
///

--2
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

--3
TEST ///--substitute SOSPoly
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3})
    S = QQ[x,y,z]
    t1 = sosPoly(S, {x+1,y}, {2,3})
    t2 = sub (s, S)
    assert (t1 == t2)
///

--4
TEST ///--toRing
    debug needsPackage "SOS"
    R = QQ[x,y];
    s = sosPoly(R, {x+1,y}, {2,3});
    S = RR[x,y];
    s2 = toRing_S s;
    assert instance(coefficientRing ring s2, RealField)
    s3 = toRing_R s2;
    assert (s==s3)
    
    tol := HighPrecision;
    f = 0.1*x_S^2 + y^2
    g = 1/10*(symbol x)_R^2 + (symbol y)_R^2
    -- comparison in rationals is complicated:
    residual = sum \\ abs \ (x -> lift (x,QQ)) \ flatten entries last coefficients (toRing_R f - g)
    assert (residual < tol)
    -- comparison in reals:
    assert (norm (toRing_S g - f) < tol)
///

--5
TEST /// --sosdec
    R=QQ[x,y,z]
    Q=matrix{{1,-1/2,1},{-1/2,1,-1/2},{1,-1/2,1}}
    Q=promote(Q,QQ)
    mon=matrix{{x^3},{x^2*z},{y*z^2}}
    f=sosPoly(mon,Q)
    assert(f=!=null and sumSOS f==transpose mon * Q *mon)
///

--6
TEST /// --chooseMons
    debug needsPackage "SOS"
    R = QQ[x,y];
    f = x^4+2*x*y-x+y^4
    lmsos = chooseMons(f)
    assert( lmsos === null )

    R = QQ[x,y]
    f = (x+2*y)^2 + (x-y)^4
    lmsos = chooseMons f
    assert( lmsos=!=null and numRows lmsos == 5 )

    R = QQ[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = chooseMons(f-t)
    assert( lmsos=!=null and ring lmsos===R and numRows lmsos == 6 )
    
    R = RR[x,y][t];
    f = x^4+2*x*y-x+y^4
    lmsos = chooseMons(f-t)
    assert( lmsos=!=null and ring lmsos===R and numRows lmsos == 6 )

    R = QQ[x,y][l,t];
    f = y + l *(y-x^2) - t
    mon = chooseMons f
    assert(mon == matrix{{1_R},{x_R}})
///

--7
TEST /// --createSOSModel
    debug needsPackage "SOS"
    eval = (Q,v) -> (transpose v * Q * v)_(0,0)
    
    R = QQ[x][t];
    f = x^4 - 2*x + t;
    mon = matrix{{1},{x},{x^2}}
    (C,Ai,p0,V,A,B,b) = createSOSModel(f,mon)
    assert( eval(C,mon) == x^4 - 2*x )
    assert( #Ai==2 and all({0,1}, j -> eval(Ai#j,mon) == V_(0,j)) )
    
    equal = (f1,f2) -> norm(f1-f2) < HighPrecision;
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

--8
TEST /// --LDLdecomposition
    debug needsPackage "SOS"
    A = matrix(QQ, {{5,3,5},{3,2,4},{5,4,10}})
    (L,D,P,err) = LDLdecomposition A
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    (L,D,P,err) = LDLdecomposition promote(A,RR)
    assert(err==0 and L*D*transpose L == transpose P * A * P)
    
    V = random(QQ^12,QQ^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and L*D*transpose L == transpose P * A * P)

    equal = (f1,f2) -> norm(f1-f2) < MedPrecision;
    V = random(RR^12,RR^8)
    A = V * transpose V 
    (L,D,P,err) = LDLdecomposition(A)
    assert(err==0 and equal(L*D*transpose L, transpose P * A * P))

    -- this matrix is not psd, but its principal minors are zero
    A = matrix(QQ,{{1,-1,1},{-1,1,1},{1,1,1}})
    (L,D,P,err) = LDLdecomposition A
    assert(err>0)
///

--9
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

--10
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

--11
TEST ///--recoverSolution
    R = RR[x,y];
    mon = matrix {{1},{x},{y}};
    X = matrix(RR, {{1,0,1},{0,0,0},{1,0,1}} );
    sol = recoverSolution(mon,X);
    assert(sol#0#1==0 and sol#1#1==1)
///

--12
TEST /// --solveSDP
    debug needsPackage "SOS"

    -- trivial cases (solved in preprocessing)
    (X,y,Z) = solveSDP (matrix{{1,0},{0,-1}},(),zeros(QQ,0,1),Solver=>"M2");
    assert(y===null and X===null);
    (X,y,Z) = solveSDP (matrix{{1,0},{0,1}},(),zeros(QQ,0,1),Solver=>"M2");
    assert(y==0);

    tests := {0,1,2,4};
    results := checkSolveSDP("M2",i->member(i,tests));
    assert all(results,t->t=!=false);
///

--13
TEST /// --solveSOS
    debug needsPackage "SOS"
    results := checkSolveSOS("M2",i->true)
    assert all(results,t->t=!=false);
///

--14
TEST /// --lowerBound
    debug needsPackage "SOS"
    tests := set{0,1,2,4,5,6};
    results := checkLowerBound("M2",i->member(i,tests))
    assert all(results,t->t=!=false);
///
