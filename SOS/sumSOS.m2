-- p = sumSOS (g,d)
-- 
-- sumSOS   compute sum of weighted squares
--    p = sumSOS(g,d) returns the sum of weighted
--    sqaures such that:
--
--    p = sum_i d_i * g_i^2 
--
--    input: (g,d):  g: sequence of polynomials
--                   d: sequence of weights     	       	    
--    output: p:     the polynomial
--   
-- Author: Helfried Peyrl
-- $Id: sumSOS.m2,v 1.2 2013-01-19 14:36:09 hpeyrl Exp $


sumSOS = (g,d) -> (
     return sum apply(toList (0..#g-1), i -> g_i^2 * d_i)
     )     

