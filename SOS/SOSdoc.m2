
document { 
     Key => SOS,
     Headline => "An SOS package",
     EM "SOS", " is a package for solving sum of squares (SOS) problems."
     }
document {
     Key => {findSOS},
     Headline => "Computation of a SOS decomposition of a polynomial",
     EM "findSOS", " uses ", TO{"solveSDP"}, " to compute an SOS ",
     "decomposition of a polynomial. It tries to obtain an exact solution by rounding ",
     "the numerical result and checking positive definiteness of the Gram matrix ",
     "afterwards. If successful the polynomial ", TT"f", " can be written as", BR{}, BR{},    
     TT "f = mon' * Q * mon", ",", BR{}, BR{}, "where ", TT"Q", " is a rational, positive ",
     "semidefinite matrix, and ", TT"mon", " is a vector of monomials.",     
     Usage => "(ok,Q,mon) = findSOS f",
     Inputs => { "f" => PolynomialRing => {"a polynomial with coefficients in ", TT "QQ"}},
     Outputs => { "ok" => Boolean => {"indicates whether a rational SOS decomposition was found"},
      "Q" => Matrix => {"the rational n by n Gram matrix of the polynomial ", TT "f"},
      "mon" => Matrix => {"a n by 1 matrix of monomials"}},
     SeeAlso => {getSOS,Solver},
     BR{},
     "Find a SOS decomposition of a given polynomial:",
     EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (ok,Q,mon) = findSOS f
     transpose(mon)*Q*mon - f
          ///
     }
document {
     Key => {getSOS},
     Headline => "SOS decomposition of a polynomial",
     EM "getSOS", " uses ", TO{"findSOS"}, " to compute a rational SOS decomposition ",
     "of a polynomial: ", BR{}, BR{},
     TT "f = sum d", SUB "i", TT " g", SUB "i", SUP "2", ",", BR{},BR{},
     "where the g", SUB "i", " are polynomials in ", TT "QQ[x]", " and the w", SUB "i", 
     " are weights in ", TT "QQ", ". The function yields an error if such a decomposition ",
     "could not be obtained.",
     Usage => "(g,d) = getSOS f",
     Inputs => { "f" => PolynomialRing => {"a polynomial with coefficients in ", TT "QQ"}},
     Outputs => { "g" => Sequence => {"of polynomials with coefficients in ", TT "QQ"},
           "d" => Sequence => {"of scalar weights in ", TT "QQ"}},
     SeeAlso => {findSOS,Solver},
     EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (g,d) = getSOS f
     sumSOS(g,d) - f
     ///,
     BR{},
     EM "getSOS", " can also solve parametric SOS problems that depend affinely of some decision variables. 
     For instance, we can find an SOS lower bound for the dehomogenized Motzkin polynomial:",
     EXAMPLE lines ///      
     R = QQ[x,z,gam];
     f = x^4+x^2+z^6-3*x^2*z^2-gam;
     (g,d,tval) = getSOS (f,{gam},-gam,rndTol=>12)  
                ///
     }
document {
     Key => {sumSOS},
     Headline => "Computes the expansion of a weighted SOS",
     EM "sumSOS", " expands a weighted SOS decomposition: ", BR{}, BR{},
     TT "f = sum d", SUB "i", TT " g", SUB "i", SUP "2", ",", BR{},BR{},
     "where g", SUB "i", " are polynomials and d", SUB "i", 
     " are scalars.",
     Usage => "f = getSOS (g,d)",
     Outputs => { "f" => RingElement => {"a polynomial"}},
     Inputs => { "g" => Sequence => {"of polynomials"},
           "d" => Sequence => {"of scalar weights"}},
     EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (g,d) = getSOS f
     sumSOS(g,d) - f
     ///
     }

document {
     Key => {rndTol},
     Headline => "rounding tolerance",
     "Minimal rounding precision in x binary digits.",
     }

document {
        Key => {LDLdecomposition},
        Headline => "LDL' factorization of a positive semidefinite matrix",
    "If ", TT "A", " is a positive semidefinite matrix, ", EM "LDLdecomposition", " returns a lower 
    triangular matrix ", TT "L", " with ones in the diagonal, a diagonal matrix ",
    TT "D", " and a permutation matrix ", TT "P", " such that ", TT "L'*D*L = P'*A*P.",
        Usage => "(L,D,P,err) = LDLdecomposition A",
        Inputs => { "A" => Matrix => {"over ", TT "QQ", " or ", TT "ZZ." } },
        Outputs => { "L" => Matrix => {"a lower triangular matrix over ", TT "QQ."},
    "D" => Matrix => {"a diagonal matrix over ", TT "QQ."},
    "P" => Matrix => {"a permutation matrix over ", TT "QQ."},
    "err" => ZZ => {"which is 0 when the factorization was successful, i.e., if ", TT "A", 
        " is positive semidefinite."}},
        -- SourceCode => {LDLdecomposition},
        EXAMPLE lines ///
          A = matrix {{5,3,5},{3,2,4},{5,4,10}}
      (L,D,P,err) = LDLdecomposition(A)
      L*D*transpose(L) == transpose(P)*A*P
      ///
        }

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
     TO2 {[solveSDP,Solver],"CSDP"},
     " is also available. ",
     Usage => "(y,X) = solveSDP(C,A,b),\n (y,X) = solveSDP(C,A,b,y0),",
     Inputs => { "C" => Matrix => {"a symmetric n by n matrix, over ", TT "RR"},
      "A" => Sequence => {"consisting of m symmetric n by n matrices over ", TT "RR"},
      "b" => Matrix => {"an m by 1 matrix over ", TT "RR"},
      "y0" => Matrix => {"an m by 1 matrix over ", TT "RR", " (optional)"}},
     Outputs => { "y" => {"an m by 1 ",TO matrix,", primal solution"},
                  "X" => {"an n by n ",TO matrix,", dual solution (not available if ", TT "Solver=>\"M2\"", " )"} },
     
     EXAMPLE lines ///
          C = matrix {{1,0},{0,2}};
          A = (matrix {{0,1},{1,0}});
      b = matrix {{1}};
          (y,X) = solveSDP(C,A,b);
          y
          ///
     }
document {
     Key => {Solver,[solveSDP,Solver]},
     Headline => "semidefinite programming solver",
     "The following SDP solvers are available:",
     UL{
       {"\"M2\"", " -- use a simple dual interior point method implemented in Macaulay2"},
       {"\"CSDP\"", " -- use the CSDP solver, available at ", TT "https://projects.coin-or.org/Csdp/" },
       },
     "The CSDP executable can be specified when loading the package, as follows ",BR{},
     TT "loadPackage(SOS,Configuration=>{\"CSDPexec\"=>\"csdp\"})",
     }
