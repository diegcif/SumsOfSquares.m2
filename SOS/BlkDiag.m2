-- B = blkdiag(A1,A2,...) 
--
-- BLKDIAG   blkdiag creates a block diagonal matrix 
--
-- Author: Helfried Peyrl
-- $Id$

newPackage(
	"BlkDiag",
    	Version => "1.0", 
    	Date => "April 27, 2007",
    	Authors => {
	     {Name => "Helfried Peyrl", Email => "peyrl@control.ee.ethz.ch"}
	     },
    	Headline => "Creating block diagonal matrices",
    	DebuggingMode => false 
    	)

export {"blkDiag"}

blkDiag = args -> (
     args = sequence args;
     if #args<2 then error ("expected at least 2 input arguments.");
     
     r := ring args#0;
     B := args#0;
     
     for i from 2 to #args do (
	  n1 := numgens source B;
     	  m1 := numgens source B;
     	  n2 := numgens source args#(i-1);
     	  m2 := numgens source args#(i-1);
	  B = matrix{{B,map(r^m1,n2,(i,j)->0_r)},{map(r^m2,r^n1,(i,j)->0),args#(i-1)}};
	  );
     return B;
     )
