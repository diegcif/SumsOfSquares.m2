
document { 
    Key => SOS,
    Headline => "An SOS package",
    TT "SOS", " is a package for solving sum-of-squares (SOS) problems.",

    HEADER5 "Computing SOS decompositions",
    "The following example demonstrates how to compute an SOS decomposition of a polynomial.",
    EXAMPLE lines ///
      R = QQ[x,y];
      f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y
      (mon,Q,X,tval) = solveSOS f;
      s = sosPoly(mon,Q)
    ///,
    "We can check with the command ", TT "sumSOS", " whether the found decomposition matches the original polynomial",
    EXAMPLE lines ///
      sumSOS(s)
    ///,

    HEADER5 "SOS with parameters",
    "If the coefficients of the polynomial are linearly parameterized, we can also search for parameters which render a polynomial to be a SOS. In the following example, the variable ", TT "t", " will be treated as a free parameter.",
    EXAMPLE lines ///
      R = QQ[x][t];
      f = (t-1)*x^4+1/2*t*x+1;
      (mon,Q,X,tval) = solveSOS (f);
      sosPoly(mon,Q)
      tval
    ///,

    HEADER5 "SOS with parameter optimization",
    "Semidefinite programming also allows to optimize a linear functional of the decision variables. As an example we compute a lower bound of a polynomial by minimizing its constant term. Since there is a tradeoff between rounding and optimality, we specify the required rounding precision as an optional input argument.",
    EXAMPLE lines ///
      R = QQ[x][t];
      f = x^4 - 2*x + t;
      (mon,Q,X,tval) = solveSOS (f,t);
      sosPoly(mon,Q)
      tval
    ///,
    }

--###################################
-- SOSPoly
--###################################

doc /// --SOSPoly
    Key
        SOSPoly
        (ring, SOSPoly)
        (gens, SOSPoly)
        (coefficients, SOSPoly)
        (length, SOSPoly)
        (net, SOSPoly)
        (substitute, SOSPoly, Ring)
        (symbol +, SOSPoly, SOSPoly)
        (symbol *, SOSPoly, SOSPoly)
        (symbol *, Number, SOSPoly)
        (symbol ^, SOSPoly, ZZ)
        (symbol ==, SOSPoly, SOSPoly)
        (symbol ==, SOSPoly, RingElement)
        (symbol ==, RingElement, SOSPoly)
    Headline
        A type to store SOS decompositions of polynomials
    Description
      Text
        A polynomial $f\in K[x]$ is a sum-of-squares (SOS) if it can be written as
        $$f = \sum_i d_i g_i^2,$$
        where the $g_i$ are polynomials in $K[x]$ and the $d_i$ are weights in $K$.
        This data type stores sums of squares in terms of the summands.  
        The type is a hash table consisting of the polynomials to be 
        squared and summed (the 'generators'), corresponding coefficients,
        and the base ring.
      Example
        R = QQ[x,y];
        s = sosPoly(R, {x+1,y}, {2,3} )
        peek s
      Text
        The ingredients of a SOS can be recovered using the expected commands:
      Example
        gens s
        ring s
        coefficients s
      Text
        The length of an SOS is the number of summands:
      Example
        length s
      Text
        Sums of squares support many common operations with polynomials:
      Example
        2 * s
        s + s
        s * s
        s == s
      Text
        The actual polynomial can be recovered using @TO sumSOS@:
      Example
        sumSOS s
      Text
        @TO SOSPoly@ supports the @TO substitute@ command.  This
        cannot be used to change the coefficient field, use @TO toRing@ for
        this purpuse.
      Example
        S = QQ[x,y,z];
        sub (s, S)
    SeeAlso
        sosPoly
///


--###################################
-- Methods
--###################################

doc /// --cleanSOS
    Key
        (clean,RR,SOSPoly)
    Headline
        Remove terms with very small coefficients from a sum of squares.
    Usage
        clean (tol, s)
    Inputs
	  tol:RR
	    the tolerance for the coefficients.
      s:SOSPoly
    Outputs
        :SOSPoly
          a cleaned up @TO SOSPoly@
    Consequences
    Description
      Text
        Given an @TO SOSPoly@ with coefficients in the reals, 
	this method removes terms with 
	coefficients smaller than the given tolerance.  It does nothing 
	on inputs with rational coefficients.
      Example
        R = RR[x,y];
        s = sosPoly(R, {x^2+.0001*x+1, y}, {2, .0001})
        clean( .001, s )
      Code
      Pre
    SeeAlso
        SOSPoly
///

doc /// --sumSOS
    Key
        sumSOS
        (sumSOS,SOSPoly)
	(sumSOS, List, List)
    Headline
        expansion of a weighted SOS decomposition
    Usage
        sumSOS(s)
	sumSOS(g,d)
    Inputs
        s:SOSPoly
	g:List
	  a list of polynomials
	d:List
	  a list of coefficients
    Outputs
        :RingElement
          a polynomial
    Consequences
    Description
      Text
        Given polynomials $g_i$ and coefficients $d_i$,
        this method computes $f = \sum_i d_i g_i^2$.
	The polynomials and coefficients can be given as lists or
	encapsulated in an object of type @TO SOSPoly@.
      Example
        R = QQ[x,y];
	sumSOS( {x+1,y}, {2,3}  )
        s = sosPoly(R, {x+1,y}, {2,3} )
        sumSOS( s )
      Code    
      Pre
    SeeAlso
        (sosPoly,Matrix,Matrix)
///

doc /// --sosPoly
    Key
        sosPoly
        (sosPoly,Ring,List,List)
        (sosPoly,List,List)
        (sosPoly,Matrix,Matrix)
    Headline
        make an SOS polynomial
    Usage
        s = sosPoly(R,polys,coeffs)
        s = sosPoly(mon,Q)
    Inputs
        R:Ring
        polys:List
          of polynomials
        coeffs:List
          of scalars
        Q:Matrix
          the Gram matrix of the polynomial f
        mon:Matrix
          a vector of monomials
    Outputs
        s:SOSPoly
    Consequences
    Description
      Text
        This method creats an object of type @TO SOSPoly@.
        An SOS polynomial can be created from a list of generators and weights.
      Example
        R = QQ[x,y];
        s = sosPoly(R, {x+1,y}, {2,3} )
      Text
        Alternatively, one can input a Gram matrix $Q$ and a vector of monomials $mon$, as produced by the method @TO solveSOS@.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        (mon,Q,X,tval) = solveSOS f;
        s = sosPoly(mon,Q)
        sumSOS(s)
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
        (solveSOS,RingElement,RingElement)
        (solveSOS,RingElement,Matrix)
        (solveSOS,RingElement,RingElement,Matrix)
    Headline
        solve a sum-of-squares problem
    Usage
        (mon,Q,X,tval) = solveSOS(f)
        (mon,Q,X,tval) = solveSOS(f,objFun)
        (Q,X,tval) = solveSOS(f,mon)
        (Q,X,tval) = solveSOS(f,objFun,mon)
    Inputs
        f:RingElement
          a polynomial
        objFun:RingElement
          a linear function of the parameters (optional)
        mon:Matrix
          a vector of monomials (optional)
    Outputs
        Q:Matrix
          the $n\times n$ Gram matrix of the polynomial f
        mon:Matrix
          a $n\times 1$ matrix of monomials
        X:Matrix
          the $n\times n$ moment matrix (dual to Q)
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
        If the rounding fails, the numerical solution is returned.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        (mon,Q,X,tval) = solveSOS f
        transpose(mon)*Q*mon - f
      Text
        The method can also solve parametric SOS problems that depend affinely of some decision variables. 
        For instance, we can find an SOS lower bound for the dehomogenized Motzkin polynomial:
      Example
        R = QQ[x,z][t];
        f = x^4+x^2+z^6-3*x^2*z^2-t;
        (mon,Q,X,tval) = solveSOS (f,-t,RndTol=>12);
        tval
      Code
      Pre
    SeeAlso
        (sosPoly,Matrix,Matrix)
        Solver
///

doc ///
    Key
    	toRing
	(toRing, Ring, RingElement)
	(toRing, Ring, SOSPoly)
    Headline
        Move a polynomial to a ring with different coefficients
    Usage
        f = toRing (R, f)
        s = toRing (R, s)
    Inputs
        R:Ring
          a polynomial ring with rational or real coefficients
        f:RingElement
          a polynomial or
        s:SOSPoly
          an @TO SOSPoly@
    Outputs
    	f:RingElement
	  the input polynomial in the new ring, or
	s:SOSPoly
	  the input SOS polynomial in the new ring
    Description
    	Text
    	    This method can be used to change a polynomial or SOSPoly with
	    rational coefficients into one with real coefficients and vice
	    versa.
    	Example
	    R = QQ[x,y];
    	    s = sosPoly(R, {x+1,y}, {2,3});
    	    S = RR[x,y];
    	    s2 = toRing_S s
    	    s3 = toRing_R s2
    Caveat
    	The function is designed to switch from real coefficients to rational
	coefficients and back.  It's behaviour is undefined if both rings
	have the same coefficient ring.  The obvious rounding issues apply.
///

doc /// --roundPSDmatrix
    Key
        roundPSDmatrix
    Headline
        rational rounding of a PSD matrix
    Usage
        (Qp,ispsd) = roundPSDmatrix(Q,A,b,d)
    Inputs
        Q:Matrix
          a symmetric matrix (real)
        A:Matrix
          (rational)
        b:Matrix
          (rational)
        d:ZZ
          the rounding precision
    Outputs
        Qp:Matrix
          the rounded matrix (rational)
        ispsd:Boolean
          true if Qp is positive semidefinite
    Consequences
    Description
      Text
        Let $S^n$ be the space of symmetric $n\times n$ matrices,
        and let $L \subset S^n$ be a rational affine subspace.
        By @TO2 {smat2vec,"vectorization"}@ we may describe this subspace in the form  $A q = b$ for some matrix $A$ with $n(n+1)/2$ columns.
        Given a real matrix $Q\in S^n$, this method finds a nearby rational matrix $Q_p$ on $L$.
      Code
      Pre
    SeeAlso
        createSOSModel
        project2linspace
        smat2vec
///

doc /// --smat2vec
    Key
        smat2vec
        (smat2vec,Matrix)
        (smat2vec,List)
        vec2smat
        (vec2smat,Matrix)
        (vec2smat,List)
    Headline
        vectorization of a symmetric matrix
    Usage
        v = smat2vec A
        A = vec2smat v
    Inputs
        A:Matrix
          symmetric
        v:Matrix
          a vector
    Outputs
        v:Matrix
        A:Matrix
    Consequences
    Description
      Text
        The method {\tt smat2vec} obtains the vectorization of a symmetric matrix.
        The method {\tt vec2smat} performs the reverse operation.
      Example
        A = matrix(QQ, {{1,2,3,4},{2,5,6,7},{3,6,8,9},{4,7,9,10}})
        v = smat2vec A
        vec2smat v
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
    Outputs
        L:Matrix
          lower triangular
        D:Matrix
          diagonal
        P:Matrix
          permutation matrix
        err:ZZ
          which is 0 when the factorization was successful, i.e., if A is positive semidefinite.
    Consequences
    Description
      Text
        Given a positive semidefinite matrix $A$, this method returns a lower triangular matrix $L$ with ones in the diagonal, a diagonal matrix $D$ and a permutation matrix $P$ such that $L D L' = P' A P.$
      Example
        A = matrix(QQ, {{5,3,5},{3,2,4},{5,4,10}})
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
          a symmetric $n\times n$ matrix
        A:Sequence
          consisting of $m$ symmetric $n\times n$ matrices
        y0:Matrix
          an $m\times 1$ matrix (optional)
    Outputs
        y:
          an $m\times 1$ matrix, primal variable
        X:
          an $n\times n$ matrix, dual variable (not available if Solver=>"M2")
        Q:
          an $n\times n$ matrix, primal variable
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

doc /// --sosdecTernary
    Key
        sosdecTernary
    Headline
       Sum of squares decomposition for ternary forms.
    Usage
        (p,q) = sosdecTernary(f, Solver=>"CSDP")
    Inputs
        f:RingElement
          a homogeneous polynomial in 3 variables
    Outputs
        p:List
          of sum of squares
        q:List
          of sum of squares
    Consequences
    Description
      Text
        Given a non-negative ternary form $f$, this method uses Hilbert's algorithm to compute a decomposition of $f$ as sum of squares of rational polynomials:
        $$f=\frac{\prod_ip_i}{\prod_iq_i}$$
        The method returns null if $f$ is not non-negative.

        {\bf References:}
        de Klerk, Etienne and Pasechnik, Dmitrii V.: Products of positive forms, linear matrix inequalities, and Hilbert 17th problem for ternary forms, European J. Oper. Res., 39-45 (2004).
    Caveat
        This implementation only works with the solver CSDP.
///

doc /// --sosInIdeal
    Key
        sosInIdeal
    Headline
       Sum of squares polynomial in ideal
    Usage
        (p,mult) = sosInIdeal(f,d,Solver=>"CSDP")
    Inputs
        f:List
          a list of polynomials
        d:ZZ
          positive integer greater than the degrees of polynomials in the list f.
    Outputs
        p:SOSPoly
        mult:List
          of polynomial multipliers
    Consequences
    Description
      Text
        This method takes a list of polynomials as an input and finds a sum of squares polynomial upto degree d in the ideal generated by the polynomials in the list f. 
        Output is null if it doesn't find a sum of squares polynomial upto that degree.
    Caveat
        This implementation only works with the solver CSDP.
///

doc /// --lowerBound
    Key
        lowerBound
        (lowerBound, RingElement)
        (lowerBound, RingElement, ZZ)
        (lowerBound, RingElement, Matrix, ZZ)
    Headline
        finds a lower bound for a polynomial
    Usage
        (bound,mon,Q,X) = lowerBound(f)
        (bound,mon,Q,X) = lowerBound(f,D)
        (bound,mon,Q,X,mult) = lowerBound(f, h, D)
    Inputs
        f:RingElement
          a polynomial
        h:List
          a polynomials (optional)
        D:ZZ
          a degree bound for the SDP relaxation (optional)
    Outputs
        bound:
          a lower bound on f
        mon:
          the minimizer of f (if it can be recovered)
    Consequences
    Description
      Text
        This method finds a lower bound for a function $f$, which is either a polynomial or a rational function.
        In some cases the minimizer can be recovered by using the method @TO recoverSolution@.
      Example
        R=QQ[x];
        f = (x-1)^2 + (x+3)^2;
        (bound,mon,Q,X) = lowerBound(f);
        bound
      Text
        More generally, the method computes lower bounds for a polynomial optimization problem of the form
        $$min_x \, f(x) \,\,\, s.t. \,\,\, h_i(x) = 0, \, i=1..m$$
      Example
        R = QQ[x,y];
        f = y;
        h = matrix{{y-x^2}};
        (bound,mon,Q,X,mult) = lowerBound (f, h, 4);
        bound
      Text
        The method also works in quotient rings.
      Example
        R = QQ[x,y]/ideal(x^2 - x, y^2 - y);
        f = x - y;
        (bound,mon,Q,X) = lowerBound(f,2);
        bound
    SeeAlso
        recoverSolution
///

doc /// --recoverSolution
    Key
        recoverSolution
    Headline
        recover the solution of an SOS problem
    Usage
        sol = recoverSolution(mon,X)
    Inputs
        mon:Matrix
          of monomials
        X:Matrix
          the moment matrix
    Outputs
        sol:List
          the solution
    Consequences
    Description
      Text
        This method attempts to find the solution of an SOS problem
        by checking if the moment matrix is rank one.
      Example
        R = RR[x,y];
        mon = matrix {{1},{x},{y}};
        X = matrix(RR, {{1,0,1},{0,0,0},{1,0,1}} );
        sol = recoverSolution(mon,X)
      Code
      Pre
    SeeAlso
        solveSOS
        lowerBound
///

doc /// --checkSolver
    Key
        checkSolver
        (checkSolver,String)
        (checkSolver,String,String)
        (checkSolver,String,Function)
    Headline
        tests an SDP solver
    Usage
        checkSolver(solver)
        checkSolver(solver,fun)
    Inputs
        solver:String
          either "M2" or "CSDP" or "SDPA"
        fun:Function
          (optional)
    Consequences
    Description
      Text
        This method tests that a function works works properly using a specified solver.
      Code
      Pre
    SeeAlso
        Solver
///

--###################################
-- Unexported methods (for developers)
--###################################

doc /// --createSOSModel
    Key
        createSOSModel
        (createSOSModel,RingElement,Matrix)
    Headline
        space of Gram matrices of a polynomial (for developers)
    Usage
        (C,Ai,Bi,A,B,b) = createSOSModel(f,mon)
    Inputs
        f:RingElement
          a polynomial
        mon:Matrix
          a vector of monomials
    Outputs
        C:Matrix
        Ai:Sequence
        Bi:Sequence
        A:Matrix
        B:Matrix
        b:Matrix
    Consequences
    Description
      Text
        This method creates the kernel and image model of the Gram matrices of a polynomial $f$.

        A Gram matrix representation of $f$ is a symmetric matrix $Q$ such that
        $f = mon' Q mon$,
        where $mon$ is a vector of monomials.
        The set of all Gram matrices $Q$ is an affine subspace.
        This subspace can be described in image form as
        $Q = C - \sum_i y_i A_i$,
        or in kernel form as
        $A q = b$
        where $q$ is the @TO2 {smat2vec,"vectorization"}@ of $Q$.

        For parametric SOS problems the image form is
        $Q = C - \sum_i y_i A_i - \sum_j p_j B_j$,
        where $p_j$ are the parameters,
        and the kernel form is
        $A q + B p = b$.
      Code
      Pre
    SeeAlso
///

doc /// --choosemonp
    Key
        choosemonp
    Headline
        create list of monomials based on the Newton polytope
    Usage
        (lmf, lmsos) = choosemonp(f)
    Inputs
        f:RingElement
          a polynomial
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

doc /// --makeMultiples
    Key	   
        makeMultiples
    Headline
        Multiply a list of polynomials by monomials
    Usage
        (H,m) = makeMultiples (h, D, homog)
    Inputs
        h:List
          a list of polynomials
        D:ZZ
          degree bound
	homog:Boolean
	  whether the whole list should be homogenous
    Outputs
        f:RingElement
	  the generic combination
	H:List
	  consisting of multiples of h
	m:List
	  consisting of monomials
    Description
      Text
        This method takes a list of polynomials as an input and multiplies them by all the monomials up to a certain degree bound.
      Example
        R = QQ[x,y,z];
    	f1 = x + y;
    	f2 = x + z;
    	(H,m) = makeMultiples ({f1,f2},2, false);
        H
///

--###################################
-- Symbols
--###################################

doc /// --RndTol
    Key
        RndTol
        [solveSOS,RndTol]
    Headline
        rounding precision
    Consequences
    Description
      Text
        Minimal rounding precision in $x$ binary digits.
      Code
      Pre
    SeeAlso
///

document { --Solver
    Key => {
        Solver,
        [solveSDP,Solver],
        [solveSOS,Solver],
        [sosInIdeal,Solver],
        [sosdecTernary,Solver],
        [lowerBound,Solver],
        },
    Headline => "semidefinite programming solver",
    "The following SDP solvers are available:",
    UL{
      {"\"M2\"", " -- use a simple dual interior point method implemented in Macaulay2"},
       {"\"CSDP\"", " -- use the CSDP solver, available at ", TT "https://projects.coin-or.org/Csdp/" },
       {"\"SDPA\"", " -- use the SDPA solver, available at ", TT "http://sdpa.sourceforge.net/" },
      },
    "The CSDP and SDPA executables can be specified when loading the package, as follows ",BR{},
    TT "loadPackage(\"SOS\",Configuration=>{\"CSDPexec\"=>\"csdp\",\"SDPAexec\"=>\"sdpa\"})",BR{},
    }

