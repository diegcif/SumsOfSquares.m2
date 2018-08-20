
document { 
    Key => SOS,
    Headline => "A package for sums-of-squares problems",
    TT "SOS", " is a package to solve sum-of-squares (SOS) problems.",
    
    HEADER4 "Introduction",
    "Writing a polynomial as a sum-of-squares proves its non-negativity for all arguments,
    but not all non-negative polynomials are sum-of-squares.  
    While non-negativity of a polynomial is hard to check, there are efficient methods to find
    sums-of-squares decompositions and this package makes some of them available in Macaulay2.  
    These methods rely on semi-definite-programming solvers from 
    mathematical optimization.  While there is a built in solver in the package, 
    it is highly recommended to configure an external ", TO Solver, ".",
    
    HEADER4 "Usage examples",
    "The most basic application is to (try to) decompose a polynomial as a sum-of-squares using the function ", TO "solveSOS",
    EXAMPLE lines ///
      R = QQ[x,y];
      f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
      sol = solveSOS f;
    ///,
    "The return value is an object of type ", TO "SDPResult", " which, in the case of success, contains in particular the SOS decomposition.",
    " It can be extracted with ", TO "sosPoly", ".  This returns an object of type ", TO "SOSPoly", 
    " which supports many operations that polynomials support.",
    EXAMPLE lines ///
      s = sosPoly sol
    ///,
    "The command ", TO "sumSOS", " can be used to check that the found decomposition matches the original polynomial:",
    EXAMPLE lines ///
      sumSOS(s)
    ///,
    
    HEADER4 "Sums of squares modulo equality constraints",
    "The package supports SOS decompositions in ",  
    TO2 {"solveSOS(RingElement,Matrix)","quotient rings"},
    ". This can be useful to prove non-negativity of a polynomial on a variety.  The following example is 
    taken from [P05].  Consider the problem
    of proving that the polynomial ", ITALIC TEX "f = 10-x^2-y", " is non-negative on the circle defined by ", ITALIC TEX "g = x^2 + y^2 - 1", ". ",
    "To do this we check if ", ITALIC TEX "f", " is a sum-of-squares in the quotient ring modulo ", ITALIC TEX "g", ". ",
    "For such a computation, a degree bound must be given by the user",
    EXAMPLE lines ///
        R = QQ[x,y];
        S = R/ideal(x^2 + y^2 - 1);
        f = 10-x^2-y;
        sol = solveSOS (f, 2);
        sosPoly sol
    ///,
    "See ", TO "TraceObj", " for how to reduce the number of summands to 2.",
    
    HEADER4 "Other cool stuff",
    UL {
	LI {"The package implements Hilbert's algorithm to decompose a non-negative ternary form into a sum-of-squares of rational functions: ", TO "sosdecTernary"},
	LI {"Sums of squares problems can be solved parametrically: ", TO "solveSOS"},
	LI {"Optimization over varieties can run using ", TO "lowerBound"},
	},
    
    HEADER4 "On the role of coefficient fields",  
    "The ", TT "SOS", " package interfaces tries to hide 
    some of the difficulties that arise from using these numerical procedures. ", 
    "The SOS package works with two coefficient rings: the rational numbers ", TO "QQ", " and the real numbers ", TO "RR",". ",
    "Almost all operations in this package rely on a numerical SDP ", TO "Solver", ".  When calling such a solver, even if the 
	input was a polynomial with rational coefficients, the result is numerical.  The package makes some 
	effort to round and return a rational result, but this can fail, independent of whether a 
	rational SOS decomposition exists or not.  In this case of failure, a real result is returned. 
	The following example of Scheiderer is SOS, but does not admit any
	rational SOS decomposition.  Consequently the package must return a real solution:",
    EXAMPLE lines ///
        f = nonnegativeForm("Scheiderer", QQ[x,y,z])
        sol = solveSOS (f, Solver=>"CSDP");
        sosPoly sol
    ///,
    "Once the rational world has been left, there is usually now way back.
    The package offers the function ", TO "clean(RR,SOSPoly)", " which from an ", TO "sosPoly", " removes
    all summands whose coefficient is smaller than a given tolerance.  This can be useful sometimes and here is how to invoke it:",
    EXAMPLE lines ///
        clean (0.001, sosPoly sol)
    ///,

    HEADER4 "Literature",
    UL {
	LI {"[BPT12] ", EM "Semidefinite Optimization and Convex Algebraic Geometry", " SIAM Textbook, edited by G. Blekherman, P. Parrilo, and R. Thomas, (2012)"},
	LI {"[P05] ", EM "Exploiting Algebraic Structure in Sum of Squares Programs", " P. Parrilo in ", EM"Positive polynomials in control", " (2005)"},
	LI {"[PP] ", EM "Computing sum-of-squares decompositions with rational coefficients", " H. Peyrl and P. Parrilo ", " in Theoretical Computer Science 409 (2008) p. 269–281"},
	}
    }

doc /// 
    Key 
        "Coefficients"
    Headline
    	Coefficient rings in the SOS package
    Description
      Text
///

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
        This data type stores SOS polynomials in terms of the summands.  
        The type is a hash table consisting of the polynomials to be 
        squared and summed (the 'generators'), corresponding coefficients,
        and the base ring.  The most common way an SOSPoly comes to life is
	as the result of an SOS decomposition.  See @TO SDPResult@ for more on this.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        sol = solveSOS f;
        sosPoly sol
      Text
        Of course one can also construct SOSPolys by specifying all their ingredients.
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
        cannot be used to change the coefficient field, though.  See @TO "Coefficients"@ for some
	of the limitations.
      Example
        S = QQ[x,y,z];
        sub (s, S)
    SeeAlso
        sosPoly
///


doc /// -- SDPResult
    Key
       SDPResult
       (net, SDPResult)
       GramMatrix
       MomentMatrix
       Parameters
    Headline
       result of an SDP computation
    Description
      Text
        This type encapsulates the result of an SDP computation.
      Example
        R = QQ[x][t];
        f = x^2 - 3*x - t;
        sol = solveSOS (f, -t, RoundTol=>12)
        peek sol
      Text
        The fields can be extracted with the operator "#"
      Example
        sol#GramMatrix
      Text
        If the Gram matrix is different from null, then the SOS polynomial can be recovered with @TO sosPoly@.
      Example
        sosPoly sol
///

--###################################
-- Methods
--###################################

doc /// --cleanSOS
    Key
        (clean,RR,SOSPoly)
    Headline
        Remove terms with very small coefficients from a sum-of-squares.
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
        (sosPoly,SDPResult)
    Headline
        make an SOS polynomial
    Usage
        s = sosPoly (SDPR) 
        s = sosPoly (mon,Q) 
        s = sosPoly (R,polys,coeffs) 
    Inputs
        SDPR:SDPResult
          the result of an SDP computation
        Q:Matrix
          positive semidefinite (Gram matrix)
        mon:Matrix
          a vector of monomials
        R:Ring
        polys:List
          of polynomials
        coeffs:List
          of scalars
    Outputs
        s:SOSPoly
    Consequences
    Description
      Text
        This method creates an object of type @TO SOSPoly@.  Very often this is applied to an
	object of type @TO SDPResult@, the result of an SDP computation.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        sosPoly solveSOS f
      Text
        One can also input a Gram matrix $Q$ and a vector of monomials $mon$.
      Example
        Q = matrix(QQ,{{1,1,1},{1,1,1},{1,1,1}});
        mon = matrix{{1},{x},{y}};
        sosPoly(mon,Q)
      Text
        Alternatively, an SOS polynomial can be created from a list of generators and weights.
      Example
        s = sosPoly(R, {x+1,y}, {2,3} )
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
    Headline
        solve a sum-of-squares problem
    Usage
        solveSOS(f) 
        solveSOS(f,objFun) 
    Inputs
        f:RingElement
          a polynomial
        objFun:RingElement
          a linear function of the parameters (optional)
    Outputs
        :SDPResult
    Consequences
    Description
      Text
        This method solves SOS problems.
        Given a rational polynomial $f(x)$, it attempts to find a rational positive semidefinite matrix $Q$ and a vector of monomials $mon$ such that
        $$f(x) = mon' Q mon.$$ 
        The algorithm first computes a floating point solution, 
        and then tries to obtain an exact solution by rounding the numerical result. 
        If the rounding fails, the numerical solution is returned.
      Example
        R = QQ[x,y];
        f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
        sol = solveSOS f;
        Q = sol#GramMatrix
        mon = sol#Monomials
        transpose(mon)*Q*mon - f
      Text
        {\bf SOS with parameters:}
        If the coefficients of the polynomial are linearly parametrized, we can search for parameters which render a polynomial to be a SOS.
        In the following example, the variable $t$ will be treated as a free parameter.
      Example
        R = QQ[x][t];
        f = (t-1)*x^4+1/2*t*x+1;
        sol = solveSOS (f);
        sosPoly(sol)
        sol#Parameters
      Text 
        {\bf SOS with parameter optimization:}
        The method also allows to optimize a linear function of the parameters.
        More precisely, given a polynomial $f(x;p)$ that depends affinely on some parameters $p$, we can solve the problem

        $$min_{p} \, objFun(p) \,\,\, s.t. \,\,\, f(x; p) \, is SOS $$

        In the following example we minimize $-t$ in order to find a lower bound for the polynomial $x^2-3x$:
      Example
        R = QQ[x][t];
        f = x^2 - 3*x - t;
        sol = solveSOS (f, -t, RoundTol=>12);
        sol#Parameters
      Text
        By default the method tries to obtain rational values of the parameters.
        Since there is a trade-off between rounding and optimality, we specify the @TO2 {RoundTol,"rounding precision"}@ as an optional input argument.
      Code
      Pre
    SeeAlso
        (solveSOS,RingElement,Matrix)
        sosPoly
        SDPResult
        Solver
///

doc /// --solveSOS (quotient ring)
    Key
        (solveSOS,RingElement,RingElement,Matrix)
        (solveSOS,RingElement,Matrix)
        (solveSOS,RingElement,ZZ)
        (solveSOS,RingElement,RingElement,ZZ)
    Headline
        sum-of-squares problem in a quotient ring
    Usage
        solveSOS(f,mon) 
        solveSOS(f,objFun,mon) 
        solveSOS(f,D) 
        solveSOS(f,objFun,D) 
    Inputs
        f:RingElement
          a polynomial
        objFun:RingElement
          a linear function of the parameters (optional)
        mon:Matrix
          a vector of monomials (alternatively a degree bound D)
    Outputs
        :SDPResult
    Consequences
    Description
      Text
        This method allows to compute SOS decompositions in quotient rings.
        A vector of monomials must be provided in this case.
      Example
        R = QQ[x,y]/ideal(x^2 + y^2 - 1);
        f = 10-x^2-y;
        mon = matrix {{1}, {x}, {y}};
        solveSOS (f, mon)
      Text
        If a degree bound $D$ is given, the method will use the vector of monomials of degree at most $D/2$.
      Example
        solveSOS (f, 2)
      Text
        Parametrized SOS problems can also be solved in quotient rings.
      Example
        S = R[t];
        solveSOS(f-t,-t,mon,RoundTol=>12)
      Code
      Pre
    SeeAlso
        solveSOS
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
        smat2vec
        RoundTol
///

doc /// --smat2vec
    Key
        smat2vec
        (smat2vec,Matrix)
        (smat2vec,List)
        vec2smat
        (vec2smat,Matrix)
        (vec2smat,List)
        Scaling
     	[smat2vec,Scaling]
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
      Text
        The scaling of the off-diagonal entries can be controlled with the optional parameter.
      Example
        smat2vec(A,Scaling=>2)
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
        {\it Matrix Computations}, Gene Golub and Charles van Loan. Johns Hopkins
        series in the Mathematical Science (1989), 2 ed., pp. 133-148.
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
        (X,y,Q) = solveSDP(C,A,b) 
        (X,y,Q) = solveSDP(C,A,b,y0) 
    Inputs
        C:Matrix
          a symmetric $n\times n$ matrix
        A:Sequence
          consisting of $m$ symmetric $n\times n$ matrices
        y0:Matrix
          a dual feasible $m\times 1$ matrix (optional)
    Outputs
        X:
          an $n\times n$ matrix, primal variable (not available if Solver=>"M2")
        y:
          an $m\times 1$ matrix, dual variable
        Q:
          an $n\times n$ matrix, dual variable
    Consequences
    Description
      Text
        This method solves a primal/dual pair of semidefinite programs.
        Given symmetric matrices $C, A_i$ and a vector $b$, the primal problem is

        $$min_{X} \, C \bullet X \,\,\, s.t. \,\,\, A_i \bullet X = b_i \, and \, X \geq 0$$

        and the dual problem is

        $$max_{y,Q} \, \sum_i b_i y_i \,\,\, s.t. \,\,\, Q = C - \sum_i y_i A_i \, and \, Q \geq 0$$

        The default algorithm is a dual interior point method implemented in M2. 
        A strictly feasible initial point $y_0$ may be provided by the user. 
        Alternatively, there is an interface to the @TO2 {[solveSDP,Solver],"solvers"}@ CSDP and SDPA.
      Example
        C = matrix {{1,0},{0,2}};
        A = matrix {{0,1},{1,0}};
        b = matrix {{-1}};
        (X,y,Q) = solveSDP(C,A,b);
        y
      Text
        {\bf References:}
        {\it Convex Optimization}, Boyd, Vandenberghe, Cambridge University Press (2004), pp. 618-619, pp. 463-466
      Code
      Pre
    Caveat
        The "M2" solver does not return the primal solution.
        It might also fail if the dual problem is not strictly feasible.
    SeeAlso
///

doc /// --sosdecTernary
    Key
        sosdecTernary
        (sosdecTernary, RingElement)
    Headline
       Sum of squares decomposition for ternary forms.
    Usage
        (p,q) = sosdecTernary(f, Solver=>"CSDP") 
    Inputs
        f:RingElement
          a homogeneous polynomial in 3 variables
    Outputs
        p:List
          of @TO SOSPoly@s
        q:List
          of @TO SOSPoly@s
    Consequences
    Description
      Text
        Given a non-negative ternary form $f$, this method uses Hilbert's algorithm to compute a decomposition of 
	$f$ as sum-of-squares of rational functions: $f=\frac{\prod_ip_i}{\prod_iq_i}$. 
        The method returns null if $f$ is not non-negative.  
	As an example, consider the homogeneous Motzkin polynomial:
      Example
         R = RR[x,y,z];
    	 f = nonnegativeForm ("Motzkin", {x,y,z});
    	 (p,q) = sosdecTernary (f, Solver=>"CSDP");
      Text
        The result, in this case, is just a quotient of two sums of squares.
      Example
         (#p, #q)
         f * sumSOS q#0 == sumSOS p#0
      Text
        {\bf References:}
        {\it Products of positive forms, linear matrix inequalities, and Hilbert 17th problem for ternary forms}, E. de Klerk, and D.V. Pasechnik, European J. Oper. Res. (2004), pp. 39-45.
    Caveat
        This implementation only works with the @TO Solver@ CSDP.
///

doc /// --sosInIdeal
    Key
        sosInIdeal
        (sosInIdeal, Matrix, ZZ)
        (sosInIdeal, Ring, ZZ)
    Headline
        Sum of squares polynomial in ideal
    Usage
        (sol,mult) = sosInIdeal(h,D,Solver=>"CSDP") 
        sol = sosInIdeal(R,D,Solver=>"CSDP") 
    Inputs
        h:Matrix
          row vector with polynomial entries
        D:ZZ
          bound on the degree of the SOS polynomial
        R:QuotientRing
          a quotient of a polynomial ring
    Outputs
        sol:SDPResult
        mult:Matrix
          column vector with polynomial multipliers
    Consequences
    Description
      Text
    	This methods finds sum-of-squares in ideals.
        It accepts two types of inputs that are useful for different purposes.
        The first invocation is to give a one row matrix with polynomial entries and a degree bound.
        The method then tries to find an SOS polynomial in the generated ideal.
        More precisely, given equations $h_1(x),...,h_m(x)$, the method looks for polynomial multipliers $h_i(x)$ such that $\sum_i l_i(x) h_i(x)$ is SOS.
     Example
        R = QQ[x,y,z];
        h = matrix {{x^2+y^2+y, y-z^2}};
        (sol,mult) = sosInIdeal (h, 2, Solver=>"CSDP");
        sosPoly sol
        h * mult == sumSOS sosPoly sol
     Text
        The second invocation is on a quotient ring, also with a degree bound.
        This tries to decompose the zero of the quotient ring as a sum-of-squares.
     Example
        S = R/ideal h;
        sol = sosInIdeal (S, 2, Solver=>"CSDP");
        sosPoly sol
        sumSOS sosPoly sol
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
        (bound,sol) = lowerBound(f) 
        (bound,sol) = lowerBound(f,D) 
        (bound,sol,mult) = lowerBound(f,h,D) 
    Inputs
        f:RingElement
          a polynomial or a rational function
        D:ZZ
          degree bound for the SDP relaxation (optional)
        h:Matrix
          row vector with polynomial entries (optional)
    Outputs
        bound:
          a lower bound on f
        sol:SDPResult
        mult:Matrix
          column vector with polynomial multipliers
    Consequences
    Description
      Text
        This method finds a lower bound for a polynomial function $f(x)$.
        More precisely, this method solves the following relaxation

        $$max_{t} \, t \,\,\, s.t. \,\,\, f(x) - t \, is SOS $$

        Lower bounds for rational functions can be obtained in a similar way.
        In some cases the minimizer can be extracted with the method @TO recoverSolution@.
      Example
        R=QQ[x];
        f = (x-1)^2 + (x+3)^2;
        (bound,sol) = lowerBound(f);
        bound
        f - bound == sosPoly sol
      Text
        By default the method tries to obtain a rational lower bound.
        Since there is a trade-off between rounding and optimality, we specify the @TO2 {RoundTol,"rounding precision"}@ as an optional input argument.
      Text
        {\bf Quotient rings:}
        Given an ideal $I$, we can also find a lower bound for $f$ on the variety of $I$.
        This can be done by constructing the associated quotient ring.
        A degree bound must be provided.
      Example
        R = QQ[x,y]/ideal(x^2 - x, y^2 - y);
        f = x - y;
        (bound,sol) = lowerBound(f,2);
        bound
        f - bound == sosPoly sol
      Text
        {\bf Avoiding quotient rings:}
        Constructing the quotient ring is sometimes too expensive since it requires Gröbner bases.
        There is an alternative (though weaker) relaxation that avoids Gröbner bases computation.
        Given equations $h_1(x),...h_m(x)$, we can look for multipliers $l_i(x)$ such that $f(x) - t + \sum_i l_i(x) h_i(x)$ is SOS.
      Example
        R = QQ[x,y];
        f = x - y;
        h = matrix{{x^2 - x, y^2 - y}};
        (bound,sol,mult) = lowerBound (f, h, 2);
        bound
        f - bound + h*mult == sumSOS sosPoly sol
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
        This method tests that a function works properly using a specified solver.
      Code
      Pre
    SeeAlso
        Solver
///

doc /// --nonnegativeForm
    Key
        nonnegativeForm
        (nonnegativeForm,String,List)
        (nonnegativeForm,String,Ring)
    Headline
        dictionary of interesting nonnegative forms
    Usage
        nonnegativeForm(name,R) 
        nonnegativeForm(name,var) 
    Inputs
        name:String
          either "Motzkin", "Robinson", "Schmuedgen", "Lax-Lax", "Choi-Lam", "Scheiderer", "Harris"
        R:Ring
          a polynomial ring
        var:List
          of variables
    Outputs
        :RingElement
          a nonnegative form
    Consequences
    Description
      Text
        This method contains a dictionary of some interesting nonnegative forms.
      Text
        The Motzkin polynomial is a ternary sextic that is not SOS.
        It was the first example of a nonnegative polynomial that is not SOS.
      Example
        R = QQ[x,y,z];
        nonnegativeForm("Motzkin", R)
      Text
        The Robinson and Schmüdgen polynomials are also ternary sextics that are not SOS.
      Example
        nonnegativeForm("Robinson", R)
        nonnegativeForm("Schmuedgen", R)
      Text
        The Lax-Lax and Choi-Lam polynomials are quaternary quartics that are not SOS.
      Example
        R = QQ[x,y,z,w];
        nonnegativeForm("Lax-Lax", R)
        nonnegativeForm("Choi-Lam", R)
      Text
        The Scheiderer polynomial is SOS over the reals, but not over the rationals.
      Example
        R = QQ[x,y,z];
        nonnegativeForm("Scheiderer", R)
      Text
        The Harris polynomial is a ternary form of degree 10 with 30 projective zeros (the largest number known in August 2018).
      Example
        nonnegativeForm("Harris", R)
      Text
        {\bf References:}
        {\it Some concrete aspects of Hilbert's 17th problem}. B. Reznick. Contemporary mathematics (2000), 253, pp. 251-272
      Code
      Pre
    SeeAlso
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

doc /// --chooseMons
    Key
        chooseMons
    Headline
        create list of monomials based on the Newton polytope
    Usage
        lmsos = chooseMons(f) 
    Inputs
        f:RingElement
          a polynomial
    Outputs
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

--###################################
-- Symbols
--###################################

doc /// -- RoundTol
    Key
        RoundTol
        [solveSOS,RoundTol]
        [lowerBound,RoundTol]
        [sosdecTernary,RoundTol]
        [sosInIdeal,RoundTol]
    Headline
        tolerance for rational rounding
    Consequences
    Description
      Text
        The optional argument {\tt RoundTol} specifies the minimal rounding precision in $d$ binary digits.
        
        SOS problems are solved numerically using an SDP solver, and afterwards the package attempts to round the floating point solution to rational numbers.
        The rounding strategy is guaranteed to work whenever the space of Gram matrices is full dimensional.
        For SOS optimization problems the rounding may cause a loss in optimality.
        The argument {\tt RoundTol} allows to control the trade-off between optimality and simplicity.
        Higher values of {\tt RoundTol} lead to better solutions.
      Example
        R = QQ[x,z];
        f = x^4+x^2+z^6-3*x^2*z^2;
        (bound,sol) = lowerBound (f,RoundTol=>4);
        bound
        (bound,sol) = lowerBound (f,RoundTol=>12);
        bound
      Text
        One can also skip the rounding by setting {\tt RoundTol => infinity}.
      Example
        (bound,sol) = lowerBound (f,RoundTol=>infinity);
        bound
      Text
        {\bf References:}
        {\it Computing sum-of-squares decompositions with rational coefficients}, H. Peyrl and P. Parrilo, in Theoretical Computer Science 409 (2008) p. 269–281
      Code
      Pre
    SeeAlso
        roundPSDmatrix
///

doc /// -- Verbose
     Key
     	[solveSDP,Verbose]
        [sosdecTernary,Verbose]
        [roundPSDmatrix,Verbose]
        [solveSOS,Verbose]
        [sosInIdeal,Verbose]
        [lowerBound,Verbose]
     Headline 
        non-essential but informative output
     Description
        Text
           Setting this option to true enables additional informative output. The default is @TO false@.
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
    Headline => "picking a semidefinite programming solver",
    "Many important computations in this package rely on an efficient SDP solver.  There
    is a very rudimentary implementation of such a solver in the Macaulay2 language.  It is 
    called the M2 solver but for most applications it will be insufficient. For this reason it is
    almost mandatory to install another solver.  The package supports csdp and sdpa which are open
    source.  This option can take the following values:" ,
    UL{
      {"\"M2\"", " -- use a simple dual interior point method implemented in Macaulay2"},
       {"\"CSDP\"", " -- use the CSDP solver, available at ", TT "https://projects.coin-or.org/Csdp/" },
       {"\"SDPA\"", " -- use the SDPA solver, available at ", TT "http://sdpa.sourceforge.net/" },
      },
    "Before any serious computation the user should install CSDP or SDPA.  In our experience CSDP gives
    the best results.  An easy way to a solver available to Macaulay2  
    is to add the executable to the PATH environment variable.  Another way is
    to explicitly specify the location of the executable when loading the package:",
    EXAMPLE lines ///
       needsPackage ("SOS", Configuration=>{"CSDPexec"=>"/some/path/csdp", "SDPAexec"=>"/some/path/sdpa"})
    ///,
    "The method ", TO "checkSolver", " can be used to check if a solver works.",
    }

doc /// --TraceObj
    Key
        TraceObj
        [solveSOS,TraceObj]
    Headline
        whether to use trace as the objective function
    Consequences
    Description
      Text
        Using the trace as the objective function is a heuristic for obtaining SOS decompositions with small number of summands.  
        Here we repeat Example 5 from [P05] and recover the shorter solution from that paper:
      Example
        R = QQ[x,y]/ideal(x^2 + y^2 - 1);
        f = 10-x^2-y;
        sosPoly solveSOS (f, 2)
        sosPoly solveSOS (f, 2, TraceObj=>true)
      Code
      Pre
    SeeAlso
///

