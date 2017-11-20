-- B = blkdiag(A1,A2,...) 
--
-- BLKDIAG   blkdiag creates a block diagonal matrix 
--
-- Author: Helfried Peyrl
-- $Id$

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

beginDocumentation()
document {
     Key => {blkDiag},
     Headline => "construct a block diagonal matrix",
     "This method returns the block diagonal matrix with blocks ",
     TT "A1,A2,...,An.",
     Usage => "D = blkDiag(A1,A2,...,An)",
     Inputs => { "Ai" => {"square matrices"} },
     Outputs => { "D" => {"block diagonal matrix"} },
     
     EXAMPLE lines ///
          A1 = matrix {{0,1},{1,0}};
          A2 = matrix {{1,2},{2,2}};
          A3 = matrix {{3}};
          blkDiag(A1,A2,A3)
          ///
     }
