export {"findSOS","getSOS","rndTol","sumSOS","createSOSModel"}

load "./findSOS.m2"
load "./createSOSModel.m2";
load "./project2linspace.m2";
load "./getRationalSOS.m2";
load "./choosemonp.m2";
load "./getSOS.m2";
load "./sumSOS.m2";


beginDocumentation()
document { 
     Key => SOS,
     Headline => "An SOS package",
     EM "SOS", " is a package for solving sum of squares (SOS) problems."
     }
document {
     Key => {findSOS},
     Headline => "Computation of a SOS decomposition of a polynomial",
     EM "findSOS", " uses the SDP solver ", TO{"solveSDP"}, " to compute an SOS ",
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
     "where g", SUB "i", " are polynomials in ", TT "QQ[x]", " and w", SUB "i", 
     " are weights in ", TT "QQ", ".",
     Usage => "f = getSOS (g,d)",
     Outputs => { "f" => PolynomialRing => {"a polynomial with coefficients in ", TT "QQ"}},
     Inputs => { "g" => Sequence => {"of polynomials with coefficients in ", TT "QQ"},
	       "d" => Sequence => {"of scalar weights in ", TT "QQ"}},
     EXAMPLE lines ///
     R = QQ[x,y];
     f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
     (g,d) = getSOS f
     sumSOS(g,d) - f
     ///
     }


