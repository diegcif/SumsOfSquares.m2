-- (ok,Q,mon) = findSOS(f)
-- 
-- findSOS   Compute a rational SOS decomposition for a polynomial f
--    (ok,Q,mon) = findSOS (f) computes a rational SOS decomposition
--    of the multivariate polynomial f such that f = mon'*Q*mon, where
--    Q is a positive semidefinite matrix and mon denotes the vector of
--    monomials. ok is of boolean 
-- input:  f:          polynomial for SOS decomp.
-- output: (ok,Q,mon): ok: boolean: successful decomposition possible?
--                     Q: dual variable
--                     mon: vector of moniomials
--
-- Author: Helfried Peyrl
-- $Id: findSOS.m2,v 1.6 2013-01-19 14:36:09 hpeyrl Exp $

load "./simpleSDP.m2";
load "./LDL.m2";
load "./BlkDiag.m2";

opts = {rndTol => -3}

findSOS = opts >> o -> args -> (
     needs "./createSOSModel.m2";
     needs "./project2linspace.m2";
     needs "./getRationalSOS.m2";
     needs "./norm2.m2"; 
     
     stdio << "findSOS by H. Peyrl and P. A. Parrilo 2007-2013" << endl;
     
     args = sequence args;
     -- f := args#0;
     p := {}; if #args >= 2 then p = args#1;
     if #args >=3 then (
	  objFcn := args#2;
	  if first degree objFcn > 1 then error("Only linear objective function allowed.")
	  );
     if #args >=4 then (
     	  lB := promote(args#3,QQ);
	  uB := promote(args#4,QQ);
	  parBounded := true;
	  ) else parBounded = false;
	      
     -- build SOS model --	 
     (C,Ai,Bi,A,B,b,mon,GramIndex,LinSpaceIndex) := createSOSModel args;

     ndim := #entries C;
     mdim := #Ai;
     
     if #p!=0 and #args>2 then (
	  -- compute an optimal solution --
	  stdio << "Solving SOS optimization problem..." << endl;
	  objFcnCF := coefficients objFcn;
	  objFcnHT := hashTable transpose {flatten entries objFcnCF_0, flatten entries objFcnCF_1};
	  obj := transpose matrix {apply(p,i->substitute(-objFcnHT#i,QQ))} || map(QQ^#Ai,QQ^1,i->0);
	  
	  if parBounded then (
	       C = blkDiag(C,diagonalMatrix(-lB),diagonalMatrix(uB));
	       Ai = apply(0..#Ai-1,i -> blkDiag(Ai_i, map(QQ^(2*#p), QQ^(2*#p), j -> 0)));
     	       Bi = apply(0..#Bi-1,i -> blkDiag(Bi_i, 
		    	 map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then 1_QQ else 0_QQ),
	       	    	 map(QQ^#p,QQ^#p, (j,k) -> if j==k and j==i then -1_QQ else 0_QQ)));
	       );
          y := - solveSDP(C, Bi | Ai, obj);
     ) else (
          -- compute a feasible solution --
     	  stdio << "Solving SOS feasibility problem..." << endl;
	  lambda := min eigenvalues (promote (C,RR), Hermitian=>true);
	  if lambda >=0 then (
	       stdio << "SDP solving maybe not necessary. Checking...." << endl;
	       (L,D,P,CnotPSD) := ldl(C);
	       if CnotPSD == 0 then return (true,C,mon);
	       );
	  obj = map(RR^(#Ai+#Bi),RR^1,i->0) || matrix{{-1_RR}};
	  y0 := map(RR^(#Ai+#Bi),RR^1,i->0) || matrix{{lambda*1.1}};
	  y = - solveSDP(C, append (Bi | Ai, id_(QQ^ndim)), obj, y0);
     );

     -- round and project --
     ynum :=  apply(0..#entries y-1, i-> round(y_(i,0)*2^52)/2^52);
     Qnum := C + sum toList apply(0..#Bi-1, i-> matrix(ynum#i*entries Bi_i)) + (
	  sum toList apply(0..#Ai-1, i-> matrix(ynum#(i+#Bi)*entries Ai_i)));
    
     if parBounded then Qnum = Qnum^{0..ndim-1}_{0..ndim-1};
     
     dhi := 52;
     d := o.rndTol;
     
     while (d < dhi) do (
	  stdio << "rounding step #" << d << endl;
  	  if #p!=0 then (
	       pVec := map(QQ^#p,QQ^1,(i,j) -> round(y_(i,0)*2^d)/2^d);
	       bPar := b - B*pVec;
	       ) else bPar= b;

     	  (ok,Qp) := getRationalSOS(Qnum,A,bPar,d,GramIndex,LinSpaceIndex);
	  if ok then break else d = d + 1;
	  );     	       
     if #p!=0 then return (ok,Qp,mon,flatten entries pVec) else
     	  return (ok,Qp,mon)
     )

