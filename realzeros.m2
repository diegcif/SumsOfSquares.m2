needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp"})
needsPackage ("NumericalAlgebraicGeometry")

sortBy = (L, fun) -> 
    last \ sort for l in L list (fun l, l);

zeroSolve = f -> (
    G := squareUp polySystem f;
    sol := solveSystem G;
    pts := for p in sol list
        if not isRealPoint p then continue
        else new Point from {Coordinates => realPart\coordinates p};
    F := matrix {f};
    for p in pts do p.Residual = norm( infinity, evaluate(F,p) );
    pts = sortBy(pts, p -> p.Residual);
    return pts;
    )

roundPoints = (d,pts) -> (
    -- d is the precision in decimal digits
    for p in pts do 
        p.Coordinates = round_d \ p.Coordinates;
    n := #pts;
    recur := {};
    for i to n-1 do
        for j from i+1 to n-1 do
            if coordinates pts#i == coordinates pts#j then
                recur = append(recur,j);
    recur = set recur;
    pts = for i to n-1 list 
        if member(i,recur) then continue else pts#i;
    return pts;
    )

realZeros = method(
     Options => {RndTol => -3, Solver=>"CSDP", Verbose => false, CleanTol => 1e-3, ResTol => 1e-2} )

realZeros(List,ZZ) := o -> (h,d) -> (
    if #h==0 then error "list of polynomials is empty";
    R := ring h#0;
    if coefficientRing R === QQ then (
        R = changeRingField(RR,R);
        h = toRing_R \ h; 
        );
    (s,mult) := sosInIdeal(h,d,RndTol=>o.RndTol,Solver=>o.Solver,Verbose=>o.Verbose);
    if s===null then
        error "SOS polynomial not found. Try increasing the degree.";
    s' := clean(o.CleanTol, s);
    h' := h | gens s';
    pts := zeroSolve(h');
    pts' := select(pts, p -> p.Residual<o.ResTol);
    return pts';
    )

R = RR[x1,x2,x3]
x4 = 1

mon6 = {x4^4, x3*x4^3, x3^2*x4^2, x3^3*x4, x3^4, x2*x4^3, x2*x3*x4^2, x2*x3^2*x4, x2*x3^3, x2^2*x4^2, x2^2*x3*x4, x2^2*x3^2, x2^3*x4, x2^3*x3, x2^4, x1*x4^3, x1*x3*x4^2, x1*x3^2*x4, x1*x3^3, x1*x2*x4^2, x1*x2*x3*x4, x1*x2*x3^2, x1*x2^2*x4, x1*x2^2*x3, x1*x2^3, x1^2*x4^2, x1^2*x3*x4, x1^2*x3^2, x1^2*x2*x4, x1^2*x2*x3, x1^2*x2^2, x1^3*x4, x1^3*x3, x1^3*x2, x1^4}

c = { 0.0178596043395922,0.0233257692779691,0.116691455690316,0.0313008384002122,0.0213029214009932,0.077794303793544,0.0939025152006359,0.255635056811919,0.0403901170088572,0.127817528405959,0.121170351026572,0.140774282093463,0.0938495213956418,0.0507647493041006,0.0259727378574786,0.0313008384002122,0.255635056811919,0.121170351026572,0.0938495213956418,0.121170351026572,0.563097128373852,0.152294247912303,0.152294247912303,0.311672854289742,0.0626218091300803,0.140774282093463,0.152294247912303,0.15583642714487,0.311672854289742,0.187865427390241,0.173320429336353,0.0626218091300803,0.115546952890901,0.0761882337453432,0.0322651411496233 };

p = sum apply(mon6,c, (i,j)->i*j);

(Q,mon,X) = solveSOS(p,Solver=>"CSDP");
s = sosdec(Q,mon)
s' = clean(1e-4, s);
pts = solveSystem gens s'
--pts := zeroSolve(h');
--pts' := select(pts, p -> p.Residual<o.ResTol);
end


-- One polynomial
R = QQ[x,y]
h1 = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1 --Motzkin
--h1 = (x^6 + y^6 + 1 ) + 3*x^2*y^2 - (x^4*y^2 + x^4 + x^2*y^4 + x^2 + y^4 + y^2) -- robinson
h = {h1}
pts = realZeros(h,8, CleanTol=>1e-3, ResTol=>1e-2)
print roundPoints(4,pts)


-- Simple example (one zero)
--R = QQ[x,y,z]
--h1 = (x-1)^2 + (y-2)^2
--h2 = (y-2)^2 + z^2
--h3 = z*(z-1)
--h = {h1,h2,h3}
--pts = realZeros(h,4, CleanTol=>1e-3, ResTol=>1e-2)
--print roundPoints(4,pts)


-- More complicated (eight zeros)
--R = QQ[x,y,z]
--h1 = 5*x^9 - 6*x^5*y + x*y^4 + 2*x*z
--h2 = -2*x^6*y + 2*x^2*y^3 + 2*y*z
--h3 = x^2 + y^2 - 17/64
--h = {h1,h2,h3}
--pts = realZeros(h,10, CleanTol=>1e-3, ResTol=>1e-2)
--print roundPoints(4,pts)
