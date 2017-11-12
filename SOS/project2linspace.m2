-- xp = project2linspace(A,b,x0)
-- 
-- PROJECT2LINSPACE   projects a rational point x0 onto the 
--                    affine subspace given by Ax=b
--
-- input:  (A,b,x0)
-- output: xp: projection of x0
--
-- Authors: Helfried Peyrl, Pablo Parrilo
-- $Id: project2linspace.m2,v 1.1 2007/03/12 14:08:13 hpeyrl Exp $

project2linspace = (A,b,x0) -> (     
     -- cast into QQ (necessary class to compute inverse)
     A2 := promote (A,QQ);
     -- ugly hack to convert b into a matrix if it is a scalar in QQ/ZZ:
     b2 := promote (matrix{{b}},QQ);
     x02 := promote (x0,QQ);
          
     -- compute projection: 
     xp := x02 - transpose(A2)*((A2*x02-b2)//(A2*transpose(A2)))
     )



