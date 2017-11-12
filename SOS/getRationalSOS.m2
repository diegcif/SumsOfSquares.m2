-- (Qp,ok) = getRationalSOS(Q,A,b,d,GramIndex,LinSpaceIndex)
-- 
-- GETRATIONALSOS compute rational SOS decomposition for given precision.
--    (Qp,ok) = getRationalSOS(Q,A,b,d) returns the projection of the
--    rounded matrix Q onto the affine subspace A*q=b. ok is true if Qp is
--    positive semidefinite. d denotes the rounding precision
--
--    GramIndex and LinSpaceIndex are hash tables for the correspondence
--    between the columns of A and the entries of Q.
--
--    Author: Helfried Peyrl
--    $Id: getRationalSOS.m2,v 1.2 2013-01-19 14:36:09 hpeyrl Exp $

needs "./LDL.m2";
     
getRationalSOS = (Q,A,b,d,GramIndex,LinSpaceIndex) -> (

     needs "./project2linspace.m2";

     ndim := #entries Q;
     
     stdio << "Rounding precision: " << d << endl;
     Q0 := matrix pack (apply(flatten entries Q, i -> round(i*2^d)/2^d),ndim);
     x0 := transpose matrix {{apply(0..numgens source A-1, i -> Q0_(toSequence (GramIndex#i-{1,1})))}};
     t := timing (xp := project2linspace(A,b,x0););
     stdio << "Time needed for projection: " << t#0 << endl;
     Q = map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then xp_(LinSpaceIndex#{i+1,j+1},0) 
     	  else xp_(LinSpaceIndex#{j+1,i+1},0));

     t = timing((L,D,P,Qpsd) := ldl(Q););
     stdio << "Time needed for LDL decomposition: " << t#0 << endl;
     if Qpsd == 0 then (true, Q) else (false,Q)
     )
