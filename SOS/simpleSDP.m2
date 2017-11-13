-- (y) = simpleSDP(C,A,b,y0)
-- 
-- SIMPLESDP   SIMPLESDP is a very simple and primitive SDP solver using
--    a pure dual interior point method based on a damped Newton steps.
--    y = simplesdp(C,A,b,y0) solves the minimization problem
--        
--             min sum( y_i * b_i )
--       s.t.  C - sum( y_i * A_i ) >=0,
--
--    where the y_i denote the decision variables with corresponding 
--    weight b_i. C and A_i are symmetric matrices. y0 is an optional
--    argument defining a starting point.
--
--    y0, y, b are vectors, C is a matrix, and A is a list of matrices.
-- 
--
--    References: 
--    Boyd, Vandenberghe: Convex Optimization, Cambridge University Press,
--    2004, pp. 618-619, pp. 463-466
--  
--    Authors: Helfried Peyrl, Johan Loefberg
--    $Id: simpleSDP.m2,v 1.5 2013-01-19 14:31:23 hpeyrl Exp $


export {"solveSDP", "untilObjNegative", "workingPrecision"}

solveSDP = method(
     TypicalValue => Matrix,
     Options => {untilObjNegative => false, workingPrecision => 53}
     )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C, sequence A, b)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
     prec := o.workingPrecision;
     R := RR_prec;
     C = promote(C,R);
     n := numgens target C;
     
     -- try to find strictly feasible starting point --
     lambda := min eigenvalues (C, Hermitian=>true);
     if lambda > 0 then (
	  y := map(R^#A,R^1,(i,j)-> 0);
	  ) else (
	  stdio << "Computing strictly feasible solution..." << endl; 
	  y =  map(R^#A,R^1,i->0) || matrix{{lambda*1.1}};
	  obj :=  map(R^#A,R^1,i->0) || matrix{{-1_R}};
	  y = solveSDP(C,append(A,id_(R^n)), obj, y, untilObjNegative=>true);
	  y = transpose matrix {take (flatten entries y,numgens target y - 1)};   
	  );
     stdio << "Computing an optimal solution..." << endl;
     solveSDP(C, A, b, y)
     )

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C, sequence A,b,y)

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y) -> (
     prec := o.workingPrecision;
     R := RR_prec;
     C = promote(C,R);
     A = apply(0..#A-1, i -> promote(A_i,R));
     b = promote(b,R);
     n := numgens target C;	       	    	      	   
     y = promote(y,R);

     m := numgens target y;
     mu := 1_R;
     theta := 10_R;
     iter := 1;
     NewtonIterMAX := 40;

     stdio << endl << "simpleSDP by H. Peyrl and J. Loefberg 2007" << endl;    
     stdio << endl << "#It:	   b'y	  dy'Hdy   mu   alpha" << endl;	    

     while mu > 0.000001 do (
	  mu = mu/theta;
     	  while true do (
	       S := C - sum toList apply(0..m-1, i-> y_(i,0) * A_i);
	       Sinv :=  solve(S, id_(target S));
	       -- compute Hessian:
	       H := map(R^m,R^m,(i,j) -> trace(Sinv*A_i*Sinv*A_j));
	       -- compute gradient:
	       g := map(R^m,R^1,(i,j) -> b_(i,0)/mu + trace(Sinv*A_i));
	       
	       -- compute damped Newton step:
	       dy := -g//H;
	       alpha := backtrack(S, -sum toList apply(0..m-1, i -> matrix(dy_(i,0) * entries A_i)));
	       y = y + transpose matrix {alpha* (flatten entries dy)};
	       lambda := (transpose dy*H*dy)_(0,0);
	       obj := transpose b * y;
	       
	       -- print some information:
	       stdio << iter << ":  " << obj << "    " << lambda << "    " << mu;
	       stdio << "    " << alpha << endl;

               iter = iter + 1;
	       if iter > NewtonIterMAX then (
		    stdio << "Warning: exceeded maximum number of iterations" << endl;
		    return y);
	       if (o.untilObjNegative == true) and (obj_(0,0) < 0)  then return y;
	       if lambda < 0.4 then break;
	       ); 
	  );
     y	        	
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
	  if cnt > BacktrackIterMAX then error ("line search did not converge.");
	  );
     alpha          
     )


beginDocumentation()
document {
     Key => {solveSDP},
     Headline => "A simple numerical SDP solver",
     EM "solveSDP", " is a pure dual interior point solver for semidefinite ",
     "programs of the form", BR{}, BR{}, TT "min sum b", SUB "i", TT " y", SUB "i", BR{}, 
     TT "s.t. C - sum y", SUB "i", TT " A", SUB "i", " > 0,", BR{}, BR{},
     "where ", TT "y", " denotes the decision variables and ", TT "C", " and ",
     TT "A", SUB "i", " are symmetric n by n matrices. A strictly feasible ",
     "initial point ", TT "y0", " may be provided by the user. If no initial point ",
     "is given, simpleSDP tries to find one.",
     Usage => "y = solveSDP(C,A,b,y0)",
     Inputs => { "C" => Matrix => {"a symmetric n by n matrix, over ", TT "RR"},
	  "A" => Sequence => {"consisting of m symmetric n by n matrices over ", TT "RR"},
	  "b" => Matrix => {"an m by 1 matrix over ", TT "RR"},
	  "y0" => Matrix => {"an m by 1 matrix over ", TT "RR"}},
     Outputs => { "y" => Matrix => {"an m by 1 matrix over ", TT "RR"}},
     
     EXAMPLE lines ///
          C = matrix {{1,0},{0,2}};
          A = (matrix {{0,1},{1,0}});
	  b = matrix {{1}};
          y = solveSDP(C,A,b)
          ///
     }
 

TEST ///
     C = matrix{{0,2,0,0,0,0},{2,0,0,0,0,0},
	  {0,0,10,0,0,0},{0,0,0,10,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
     A1 = matrix{{-1,0,0,0,0,0},{0,0,0,0,0,0},
	  {0,0,1,0,0,0},{0,0,0,0,0,0},{0,0,0,0,-1,0},{0,0,0,0,0,0}};
     A2 = matrix{{0,0,0,0,0,0},{0,-1,0,0,0,0},
	  {0,0,0,0,0,0},{0,0,0,1,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1}};
     A = (A1,A2);
     y0 = matrix{{7},{9}};
     b = matrix{{1},{1}};
     y = solveSDP(C,A,b,y0);
     yopt = matrix{{2.},{2.}};
     
    assert ( sqrt( sum apply(toList (0..#entries y-1), i-> (y_(i,0)-yopt_(i,0))^2)) <
	 0.001 * (1. + sqrt(sum apply(flatten entries yopt, i->i^2))) )
///
