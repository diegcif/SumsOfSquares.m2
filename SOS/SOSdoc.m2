
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

--###################################
-- Methods
--###################################

doc /// --sumSOS
    Key
        sumSOS
    Headline
        expansion of a weighted SOS decomposition
    Usage
        sumSOS(g,d)
    Inputs
        g:Sequence
          of polynomials
        d:Sequence
          of numbers
    Outputs
        :RingElement
          a polynomial
    Consequences
    Description
      Text
        Given polynomials $g_i$ and scalars $d_i$,
        this method computes
        $f = \sum_i d_i g_i^2$.
      Example
        R = QQ[x,y];
        sumSOS( (x+1,y), (2,3) )
      Code
      Pre
    SeeAlso
        sosdec
///

doc /// --sosdec
    Key
        sosdec
    Headline
        SOS decomposition of a polynomial
    Usage
        (g,d) = sosdec(Q,mon)
    Inputs
        Q:Matrix
          the rational $n\times n$ Gram matrix of the polynomial f
        mon:Matrix
          a $n\times 1$ matrix of monomials
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
        (solveSOS,RingElement)
        (solveSOS,RingElement,List)
        (solveSOS,RingElement,List,RingElement)
        (solveSOS,RingElement,List,RingElement,List)
    Headline
        solve a sum-of-squares problem
    Usage
        (ok,Q,mon) = solveSOS f
        (ok,Q,mon,tval) = solveSOS(f,p,objFun)
        (ok,Q,mon,tval) = solveSOS(f,p,objFun,bounds)
    Inputs
        f:RingElement
          a polynomial with coefficients in $\QQ$
        p:List
          of parameters (optional)
        objFun:RingElement
          a polynomial with coefficients in $\QQ$ (optional)
        bounds:List
          a lower and upper bound for the parameters (optional)
    Outputs
        ok:Boolean
          indicates whether a rational SOS decomposition was found
        Q:Matrix
          the rational $n\times n$ Gram matrix of the polynomial f
        mon:Matrix
          a $n\times 1$ matrix of monomials
        tval:List
          of parameter values
    Consequences
    Description
      Text
        This method solves SOS problems.
        Given a rational polynomial $f$, it attempts to find a rational positive semidefinite matrix $Q$ and a vector of monomials $mon$ such that
        $$f = mon' Q mon.$$ 
        The algorithm first computes a floating point solution, 
        and then tries to obtain an exact solution by rounding the numerical result. 
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
        Returns the projection of the rounded matrix Q onto the affine subspace $A q = b$.

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
        Projects a rational point $x_0$ onto the affine subspace given by $A x = b$
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
        (C,Ai,Bi,A,B,b,mon,GramIndex,LinSpaceIndex) = createSOSModel(f,p)
    Inputs
        f:RingElement
          a polynomial
        p:List
          of parameters
    Outputs
        C:Matrix
        Ai:Sequence
        Bi:Sequence
        A:Matrix
        B:Sequence
        b:Matrix
        mon:Matrix
        GramIndex:HashTable
        LinSpaceIndex:HashTable
    Consequences
    Description
      Text
        This method creates the kernel and image model of the Gram matrices of a polynomial $f$.

        A Gram matrix representation of $f$ is a symmetric matrix X such that
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

doc /// --LDLdecomposition
    Key
        LDLdecomposition
    Headline
        LDL factorization of a positive semidefinite matrix
    Usage
        (L,D,P,err) = LDLdecomposition A
    Inputs
        A:Matrix
          over $\QQ$ or $\ZZ$
    Outputs
        L:Matrix
          lower triangular
        D:Matrix
          diagonal
        L:Matrix
          lower triangular
        P:Matrix
          permutation matrix
        err:ZZ
          which is 0 when the factorization was successful, i.e., if A is positive semidefinite.
    Consequences
    Description
      Text
        Given a positive semidefinite matrix $A$, this method returns a lower triangular matrix $L$ with ones in the diagonal, a diagonal matrix $D$ and a permutation matrix $P$ such that $L' D L = P' A P.$
      Example
        A = matrix {{5,3,5},{3,2,4},{5,4,10}}
        (L,D,P,err) = LDLdecomposition(A)
        L*D*transpose(L) == transpose(P)*A*P
      Text
        {\bf References:}
        Gene Golub and Charles van Loan: Matrix Computations, Johns Hopkins
        series in the Mathematical Science, 2 ed., pp. 133-148,
        Baltimore Maryland, 1989.
      Code
      Pre
    SeeAlso
///

doc /// --blkDiag
    Key
        blkDiag
    Headline
        construct a block diagonal matrix
    Usage
        D = blkDiag(A1,A2,...,An)
    Inputs
        Ai:
          square matrices
    Outputs
        D:
          block diagonal matrix
    Consequences
    Description
      Text
        This method returns the block diagonal matrix with blocks 
        $A1,A2,...,An.$
      Example
        A1 = matrix {{0,1},{1,0}};
        A2 = matrix {{1,2},{2,2}};
        A3 = matrix {{3}};
        blkDiag(A1,A2,A3)
      Code
      Pre
    SeeAlso
///

doc /// --solveSDP
    Key
        solveSDP
        (solveSDP,Matrix,Matrix,Matrix)
        (solveSDP,Matrix,Matrix,Matrix,Matrix)
        (solveSDP,Matrix,Sequence,Matrix)
        (solveSDP,Matrix,Sequence,Matrix,Matrix)
    Headline
        solve a semidefinite program
    Usage
        (y,X,Q) = solveSDP(C,A,b)
        (y,X,Q) = solveSDP(C,A,b,y0)
    Inputs
        C:Matrix
          a symmetric $n\times n$ matrix over $\RR$
        A:Sequence
          consisting of $m$ symmetric $n\times n$ matrices over $\RR$
        y0:Matrix
          an $m\times 1$ matrix over $\RR$ (optional)
    Outputs
        y:
          an $m\times 1$ matrix, the primal solution
        X:
          an $n\times n$ matrix, the dual solution (not available if Solver=>"M2")
    Consequences
    Description
      Text
        This method solves a semidefinite program of the form 

        $$min_{y,Q} \, \sum_i b_i y_i \,\,\, s.t. \,\,\, Q = C - \sum_i y_i A_i \, and \, Q \geq 0$$

        where $y,Q$ are the decision variables and $C, A_i$ are symmetric $n\times n$ matrices. 
        A strictly feasible initial point $y0$ may be provided by the user. 
        The default algorithm is a dual interior point method implemented in M2. 
        Alternatively, there is an interface to the @TO2 {[solveSDP,Solver],"solvers"}@ CSDP and SDPA.
      Example
        C = matrix {{1,0},{0,2}};
        A = matrix {{0,1},{1,0}};
        b = matrix {{1}};
        (y,X,Q) = solveSDP(C,A,b);
        y
      Text
        {\bf References:}
        Boyd, Vandenberghe: Convex Optimization, Cambridge University Press,
        2004, pp. 618-619, pp. 463-466
      Code
      Pre
    Caveat
        Then "M2" solver might fail to compute the solution if the problem is not strictly feasible.
    SeeAlso
///

doc /// --checkSolver
    Key
        checkSolver
    Headline
        tests method "solveSDP" (for developers)
    Usage
        checkSolver solver
    Inputs
        solver:String
          either "M2" or "CSDP"
    Consequences
    Description
      Text
        This function tests that @TO solveSDP@ works properly.
      Code
      Pre
    SeeAlso
        Solver
///

--###################################
-- Symbols
--###################################

doc /// --rndTol
    Key
        rndTol
        [solveSOS,rndTol]
    Headline
        construct a block diagonal matrix
    Consequences
    Description
      Text
        Minimal rounding precision in $x$ binary digits.
      Code
      Pre
    SeeAlso
///

document { --Solver
    Key => {Solver,[solveSDP,Solver],[solveSOS,Solver]},
    Headline => "semidefinite programming solver",
    "The following SDP solvers are available:",
    UL{
      {"\"M2\"", " -- use a simple dual interior point method implemented in Macaulay2"},
       {"\"CSDP\"", " -- use the CSDP solver, available at ", TT "https://projects.coin-or.org/Csdp/" },
       {"\"SDPA\"", " -- use the SDPA solver, available at ", TT "http://sdpa.sourceforge.net/" },
      },
    "The CSDP and SDPA executables can be specified when loading the package, as follows ",BR{},
    TT "loadPackage(SOS,Configuration=>{\"CSDPexec\"=>\"csdp\",\"SDPAexec\"=>\"sdpa\"})",BR{},
    }

