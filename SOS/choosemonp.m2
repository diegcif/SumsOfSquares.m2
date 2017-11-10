-- (lmf, lmsos) = choosemonp(f,[p])
-- 
-- CHOOSEMON   Creates a list of monomials for a SOS decomposition of a
--    polynomial with Newton Polytope
--
-- Input:  (f,[p]): f: polynomial,
--                  p: (optional) list of linear parameters
-- Output: (lmf,lmsos): List of f's monomials, List of monomials for SOS decomp.
--
-- Authors: Helfried Peyrl, Pablo Parrilo
-- $Id: choosemonp.m2,v 1.3 2013-01-19 14:36:09 hpeyrl Exp $

needsPackage "FourierMotzkin";

choosemonp = args -> (
             
     args = sequence args;
     f := args#0;
     p := {}; if #args == 2 then p = args#1;

     -- Get rid of parameters in polynomial:
     genpos := positions(gens ring f,i->not any(p,j->j==i));
     ringf := QQ[(gens ring f)_genpos];
     n := #genpos;
     p1 := apply(p,i->i=>1);
     lmf := unique(apply(flatten entries (coefficients f)_0,i->(substitute(i,p1))));
     falt := sum lmf;
     
     -- Get exponent-points for newton polytope:
     points := substitute(matrix (transpose exponents falt)_genpos,QQ);
     maxdeg := first degree falt;
     mindeg := floor first min entries (transpose points*matrix map(ZZ^n,ZZ^1,i->1));
     maxdegs := apply(entries points, i-> max i);
     mindegs := apply(entries points, i-> min i);
     
     -- Regard exponent-points in a possible subspace
     numpoints := #entries transpose points;
     shift := first entries transpose points;
     V := matrix transpose apply(entries transpose points, i -> i - shift);
     basV := mingens image V;
     basVdim := numgens image basV;
     if basVdim != n then T := id_(QQ^n)//basV else T = id_(QQ^n);
     basVtrans := mingens kernel transpose basV;
     
     -- Get candidate points from basis of f:
     mon := flatten apply( toList(mindeg..maxdeg), k -> flatten entries basis(k, ringf));
     cp := apply (mon, i -> flatten exponents (i));
     stdio << "#candidate points: " << #cp << endl;

     -- Filter candidate points:
     -- Only the even ones within the box of degrees[mindegs:maxdegs]:
     cpf := select(cp,i-> all(i,even) and all(i-mindegs,j->j>=0) and all(maxdegs-i,j->j>=0)); 
     stdio << "#points (even and within box of polynomial degrees): " << #cpf << endl;
     -- Drop points that do not live on the subspace: 
     cpf2 := select(cpf,i-> matrix{i-shift}*basVtrans==0);
     stdio << "#points in subspace of exponent-points: " << #cpf2 << endl;
     
     -- Compute convex hull:
     liftedpts := T*V || map (QQ^1,QQ^(size falt),i->1);
     polytope := transpose substitute(first fourierMotzkin(liftedpts),QQ);

     -- Find points within the polytope:
     lexponents := cpf2_(positions(cpf2, i-> 
	       max flatten entries (polytope * ((T * transpose matrix {i-shift})||1)) <=0))/2;
     lmSOS := apply(lexponents, i-> product(n,j->(
        assert (denominator i#j==1);
         (ring f)_(genpos#j)^(numerator i#j)
         )));
     stdio << "#points inside Newton polytope: " << #lmSOS << endl;
     
     -- Workaround: cast generators back from ringf to ring f!!!!
     -- dummyring = ring f;
     
     (lmf,lmSOS)
     )
