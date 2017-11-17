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

csdpexec=((options SOS).Configuration)#"CSDPexec"

solveSDP = method(
     Options => {untilObjNegative => false, workingPrecision => 53, Algorithm=>"M2"}
     )

solveSDP(Matrix, Matrix, Matrix) := o -> (C,A,b) -> solveSDP(C, sequence A, b)

solveSDP(Matrix, Matrix, Matrix, Matrix) := o -> (C,A,b,y) -> solveSDP(C, sequence A,b,y)

solveSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
    y:=null; X:=null; Z:= null;
    if o.Algorithm == "M2" then
        y = simpleSDP(C,A,b,untilObjNegative=>o.untilObjNegative)
    else if o.Algorithm == "CSDP" then
        (y,X,Z) = solveCSDP(C,A,b)
    else
        error "unknown algorithm";
    return (y,X);
)

solveSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y0) -> (
    if o.Algorithm != "M2" then
        return solveSDP(C,A,b,o);
    y := simpleSDP(C,A,b,y0,untilObjNegative=>o.untilObjNegative);
    return (y,null);
)

--##############################
--####### SIMPLE SDP ###########
--##############################

simpleSDP = method(
     TypicalValue => Matrix,
     Options => {untilObjNegative => false}
     )

simpleSDP(Matrix, Sequence, Matrix) := o -> (C,A,b) -> (
     R := RR;
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
	  y = simpleSDP(C,append(A,id_(R^n)), obj, y, untilObjNegative=>true);
	  y = transpose matrix {take (flatten entries y,numgens target y - 1)};   
	  );
     stdio << "Computing an optimal solution..." << endl;
     simpleSDP(C, A, b, y)
     )

simpleSDP(Matrix, Sequence, Matrix, Matrix) := o -> (C,A,b,y) -> (
     R := RR;
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
           if H==0 then error "Hessian is zero";
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


--##############################
--########### CSDP #############
--##############################

solveCSDP = (C,A,b) -> (
    n := numColumns C;
    fin := getFileName();
    fout := getFileName();
    fout2 := getFileName();
    writeSDPA(fin,C,A,b);
    r := run(csdpexec | " " | fin | " " | fout);
    if r == 32512 then error "csdp executable not found";
    r = run("cat " | fout | " | tr -d + > " | fout2);
    (y,X,Z) := readSDPA(fout2,n);
    return (y,X,Z);
)

getFileName = () -> (
     filename := temporaryFileName();
     while fileExists(filename) or fileExists(filename|".txt") do filename = temporaryFileName();
     return filename
)

writeSDPA = (fname,C,A,b) -> (
    m := length A;
    n := numColumns C;
    A = prepend(C,A);
    f := openOut fname;
    inputMatrix := l -> (
        a := -A_l;
        pref := toString l | " 1 ";
        for i to n-1 do
            for j from i to n-1 do
                if a_(i,j)!=0 then
                    f << pref | toString(i+1) | " " | toString(j+1) | " " | toString a_(i,j) << endl;
    );
    f << "*SDPA file generated by SOSm2" << endl;    
    f << toString m << " =mdim" << endl;
    f << "1 =nblocks" << endl;
    f << toString n << endl;
    f << demark(" ", toString\flatten entries b) << endl;
    for l to m do(
        inputMatrix(l);
    );
    f << close;
)

readSDPA = (fname,n) -> (
    sdpa2matrix := s -> (
        e := for i in s list (i_2-1,i_3-1) => i_4;
        e' := for i in s list (i_3-1,i_2-1) => i_4;
        return map(RR^n, RR^n, e|e');
    );
    readLine := l -> for s in separate(" ",l) list if s=="" then continue else value s;
    L := lines get openIn fname;
    y := transpose matrix{readLine L_0};
    S := readLine \ drop(L,1);
    S1 := select(S, l -> l_0==1);
    S2 := select(S, l -> l_0==2);
    Z := sdpa2matrix(S1); -- slack matrix
    X := sdpa2matrix(S2); -- dual solution
    return (y,X,Z);
)


beginDocumentation()
document {
     Key => {solveSDP,(solveSDP,Matrix,Matrix,Matrix),(solveSDP,Matrix,Matrix,Matrix,Matrix),(solveSDP,Matrix,Sequence,Matrix),(solveSDP,Matrix,Sequence,Matrix,Matrix)},
     Headline => "solve a semidefinite program",
     "This method solves a semidefinite program of the form ",
     BR{}, BR{}, TT "min sum b", SUB "i", TT " y", SUB "i", BR{}, 
     TT "s.t. C - sum y", SUB "i", TT " A", SUB "i", " > 0,", BR{}, BR{},
     "where ", TT "y", " denotes the decision variables and ", TT "C", " and ",
     TT "A", SUB "i", " are symmetric n by n matrices. A strictly feasible ",
     "initial point ", TT "y0", " may be provided by the user. ",
     "The default algorithm is a dual interior point method implemented in M2, but an interface to ",
     TO2 {[solveSDP,Algorithm],"CSDP"},
     " is also available. ",
     Usage => "(y,X) = solveSDP(C,A,b),\n (y,X) = solveSDP(C,A,b,y0),",
     Inputs => { "C" => Matrix => {"a symmetric n by n matrix, over ", TT "RR"},
	  "A" => Sequence => {"consisting of m symmetric n by n matrices over ", TT "RR"},
	  "b" => Matrix => {"an m by 1 matrix over ", TT "RR"},
	  "y0" => Matrix => {"an m by 1 matrix over ", TT "RR", " (optional)"}},
     Outputs => { "y" => {"an m by 1 ",TO matrix,", primal solution"},
                  "X" => {"an n by n ",TO matrix,", dual solution (not available if ", TT "Algorithm=>\"M2\"", " )"} },
     
     EXAMPLE lines ///
          C = matrix {{1,0},{0,2}};
          A = (matrix {{0,1},{1,0}});
	  b = matrix {{1}};
          (y,X) = solveSDP(C,A,b);
          y
          ///
     }
document {
     Key => {[solveSDP,Algorithm]},
     Headline => "semidefinite programming solver",
     "The following SDP solvers are available:",
     UL{
	   {"\"M2\"", " -- use a simple dual interior point method implemented in Macaulay2"},
	   {"\"CSDP\"", " -- use the CSDP solver, available at ", TT "https://projects.coin-or.org/Csdp/" },
	   },
     "The CSDP executable can be specified when loading the package, as follows ",BR{},
     TT "loadPackage(SOS,Configuration=>{\"CSDPexec\"=>\"csdp\"})",
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
     (y,X) = solveSDP(C,A,b,y0);
     yopt = matrix{{2.},{2.}};
     
    assert ( sqrt( sum apply(toList (0..#entries y-1), i-> (y_(i,0)-yopt_(i,0))^2)) <
	 0.001 * (1. + sqrt(sum apply(flatten entries yopt, i->i^2))) )
///
