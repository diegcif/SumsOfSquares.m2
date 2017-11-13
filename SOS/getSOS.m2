-- (g,d) = getSOS(p)
-- 
-- getSOS   Get rational SOS decomposition of a polynomial f
--    (g,d) = getSOS(f) returns the rational SOS decomposition
--    of the multivariate polynomial f:
--
--    f = sum_i d_i * g_i^2 
--
--    input: f:      the polynomial
--    output: (g,d): g: sequence of polynomials
--                   d: sequence of weights     	       	    
--   
-- Author: Helfried Peyrl
-- $Id: getSOS.m2,v 1.4 2013-01-19 14:36:09 hpeyrl Exp $

needs "./LDL.m2";

getSOS = args -> (
     needs "./findSOS.m2";
          
     if #args!=1 then (ok,Q,mon,pVec) := findSOS args else (ok,Q,mon) = findSOS args;
          
     (L,D,P,err) := LDLdecomposition(Q);
     if err != 0 then error ("Gram Matrix is not positive semidefinite");
     n := #entries Q;
     g := toList flatten entries (transpose mon * transpose inverse P * L);
     d := apply(toList(0..n-1),i->D_(i,i));
     idx := positions (d, i->i!=0);
     d = d_idx;
     g = g_idx;
          
     if #args==1 then return (g,d) else return (g,d,pVec)
     )     

