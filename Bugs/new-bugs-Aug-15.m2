needsPackage ("SOS" , Configuration => { "CSDPexec" => "SDPsolvers/csdp"})

-- the following two are not solved because dual infeasible
R = QQ[x,y,z]
f = (x+y+z)^2 -- dual infeasible without calling csdp?
solveSOS (f, Solver=>"CSDP")

R = QQ[x,y,z,w]
f = 1/4 * (x+y+z)^2 + 3*(x-w^3+2*z^2)^2
solveSOS (f, Solver=>"CSDP")

-- This one is not solved but prints "Success: Dual infeasible"
R = QQ[x,y,z]
(a,b,c,d,e) = (16,-36,20,57,-38)
f4 = a*( x^10 + y^10 + z^10)+
    b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
    c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) +
    d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
    e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4)
h = x^2 + y^2 + z^2
(mon4,Q4,X4,tval) = solveSOS (f4*h, Solver=>"CSDP")
g4 = sosPoly (mon4,Q4)


-- Prints "Newton polytope has odd vertices", although sum of squares:
f = (x+2*y)^2 + (x-y)^4
solveSOS f





