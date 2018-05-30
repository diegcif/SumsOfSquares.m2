
document { 
    Key => SOS,
    Headline => "An SOS package",
    EM "SOS", " is a package for solving sum-of-squares (SOS) problems.",
    EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (ok,Q,mon) = solveSOS f;
     (g,d) = sosdec(Q,mon)
     sumSOS(g,d) - f
    ///,
    }

document {
     Key => {sumSOS},
     Headline => "Computes the expansion of a weighted SOS",
     EM "sumSOS", " expands a weighted SOS decomposition: ", BR{}, BR{},
     TT "f = sum d", SUB "i", TT " g", SUB "i", SUP "2", ",", BR{},BR{},
     "where g", SUB "i", " are polynomials and d", SUB "i", 
     " are scalars.",
     Usage => "sumSOS(g,d)",
     Outputs => { "f" => RingElement => {"a polynomial"}},
     Inputs => { "g" => Sequence => {"of polynomials"},
           "d" => Sequence => {"of scalar weights"}},
     EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (ok,Q,mon) = solveSOS f
     (g,d) = sosdec(Q,mon)
     sumSOS(g,d) - f
     ///
     }

document {
     Key => {rndTol},
     Headline => "rounding tolerance",
     "Minimal rounding precision in x binary digits.",
     }

--  References: 
--  Gene Golub and Charles van Loan: Matrix Computations, Johns Hopkins
--  series in the Mathematical Science, 2 ed., pp. 133-148,
--  Baltimore Maryland, 1989.
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

--  References: 
--  Boyd, Vandenberghe: Convex Optimization, Cambridge University Press,
--  2004, pp. 618-619, pp. 463-466
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

doc /// --sosdec
    Key
        sosdec
    Headline
        SOS decomposition of a polynomial
    Usage
        (g,d) = sosdec(Q,mon)
    Inputs
        Q:Matrix
          the rational $n \times n$ Gram matrix of the polynomial f
        mon:Matrix
          a $n x 1$ matrix of monomials
    Outputs
        g:Sequence
          of polynomials with coefficients in $\QQ$
        d:Sequence
          of scalar weights in $\QQ$
        tval:List
          of parameter values
    Consequences
    Description
      Text
        This method computes a rational SOS decomposition of a polynomial:
        $$f = \sum_i d_i g_i^2$$
        where the $g_i$ are polynomials in $\QQ[x]$ and the $d_i$ are weights in $\QQ$.
        The input is a Gram matrix $Q$ and a vector of monomials $mon$, as produced by the method @TO solveSOS@.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        (ok,Q,mon) = solveSOS f;
        (g,d) = sosdec(Q,mon)
        sumSOS(g,d) - f
      Code
      Pre
    SeeAlso
        solveSOS
        Solver
///

doc /// --solveSOS
    Key
        solveSOS
    Headline
        solve a sum-of-squares problem
    Usage
        (ok,Q,mon) = solveSOS f
        (ok,Q,mon,tval) = solveSOS(f,p,objFun)
    Inputs
        f:RingElement
          a polynomial with coefficients in $\QQ$
        p:List
          of parameters (optional)
        objFun:RingElement
          a polynomial with coefficients in $\QQ$ (optional)
    Outputs
        ok:Boolean
          indicates whether a rational SOS decomposition was found
        Q:Matrix
          the rational $n x n$ Gram matrix of the polynomial f
        mon:Matrix
          a $n x 1$ matrix of monomials
        tval:List
          of parameter values
    Consequences
    Description
      Text
        This method solves SOS problems.
        Given a rational polynomial $f$, it attempts to find a rational positive semidefinite matrix $Q$ and a vector of monomials $mon$ such that
        $$f = mon' Q mon.$$ 
        The algorithm first computes a floating point solution, 
        and tries to obtain an exact solution by rounding the numerical result and checking positive definiteness of the Gram matrix afterwards. 
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        (ok,Q,mon) = solveSOS f
        transpose(mon)*Q*mon - f
      Text
        The method can also solve parametric SOS problems that depend affinely of some decision variables. 
        For instance, we can find an SOS lower bound for the dehomogenized Motzkin polynomial:
      Example
        R = QQ[x,z,t];
        f = x^4+x^2+z^6-3*x^2*z^2-t;
        (ok,Q,mon,tval) = solveSOS (f,{t},-t,rndTol=>12);
        tval
      Code
      Pre
    SeeAlso
        sosdec
        Solver
///

doc /// --getRationalSOS
    Key
        getRationalSOS
    Headline
        compute rational SOS decomposition for given precision
    Usage
        (Qp,ok) = getRationalSOS(Q,A,b,d)
        (Qp,ok) = getRationalSOS(Q,A,b,d,GramIndex,LinSpaceIndex)
    Inputs
        Q:Matrix
          Gram matrix to be rounded
        A:Matrix
        b:Matrix
          a vector
        d:RR
          the rounding precision
    Outputs
        Qp:Matrix
          the rounded matrix
        ok:Boolean
          true if Qp is positive semidefinite
    Consequences
    Description
      Text
        Returns the projection of the rounded matrix Q onto the affine subspace $Aq=b$.

        GramIndex and LinSpaceIndex are hash tables for the correspondence between the columns of $A$ and the entries of $Q$.
      Code
      Pre
    SeeAlso
        createSOSModel
        project2linspace
///

doc /// --choosemonp
    Key
        choosemonp
    Headline
        create list of monomials based on the Newton polytope
    Usage
        (lmf, lmsos) = choosemonp(f,p)
    Inputs
        f:RingElement
          a polynomial
        p:List
          of parameters
    Outputs
        lmf:List
          of monomials of f
        lmsos:List
          of monomials for the SOS factors
    Consequences
    Description
      Text
        Creates a list of monomials for an SOS decomposition.
        The monomials are chosen based on the Newton polytope.
      Code
      Pre
    SeeAlso
///

doc /// --project2linspace
    Key
        project2linspace
    Headline
        project a rational point onto affine subspace
    Usage
        xp = project2linspace(A,b,x0)
    Inputs
        A:Matrix
        b:Matrix
          a vector
        x0:Matrix
          a rational vector
    Outputs
        xp:Matrix
          the projection of x0
    Consequences
    Description
      Text
        Projects a rational point $x_0$ onto the affine subspace given by $Ax=b$
      Code
      Pre
    SeeAlso
///

doc /// --createSOSModel
    Key
        createSOSModel
    Headline
        model of the Gram matrix representations of a polynomial
    Usage
        (C,Ai,mon,A,b,GramIndex) = createSOSModel(f)
    Inputs
        f:RingElement
          a polynomial
    Outputs
        C:Matrix
    Consequences
    Description
      Text
        This method creates the kernel and image model of the Gram matrices of a polynomial f.

        A Gram matrix representation of f is a symmetric matrix X such that
        $f = mon' X mon$,
        where $mon$ is a vector of monomials.
        The set of all Gram matrices $X$ is an affine subspace.
        This affine subspace can be described in image form as
        $X = C - \sum_i y_i A_i$
        where $y_i$ are free parameters,
        or in kernel form as
        $A x = b$
        where $x_i = X_{GramIndex#i}$.
      Code
      Pre
    SeeAlso
///
