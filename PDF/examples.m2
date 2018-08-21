-- This file contains code for the examples in the paper:
needsPackage "SOS"

R = QQ[x,y];
f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
sosPoly solveSOS f

R = QQ[x,y,z]
f = nonnegativeForm ("Motzkin", {x,y,z})
r = sosdecTernary

R = QQ[x,y]/ideal(x^2 + y^2 - 1);
f = 10-x^2-y;
sosPoly solveSOS (f, 2)
sosPoly solveSOS (f, 2, TraceObj=>true)
