needsPackage( "SOS", Configuration=>{"CSDPexec"=>"SDPsolvers/csdp","SDPAexec"=>"SDPsolvers/sdpa"} )

-- Motzkin polynomial
R = QQ[x,y,z]
h = x^2 + y^2 + z^2
f1 = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2 --Motzkin
sol1 = solveSOS (f1*h, Solver=>"CSDP")
g1 = sosPoly sol1

--Robinson
f2 = x^6 + y^6 + z^6 - (x^4*y^2 + x^2*y^4 + x^4*z^2 + x^2*z^4 + y^4*z^2 + y^2*z^4) + 3*x^2*y^2*z^2 --Robinson
sol2 = solveSOS (f2*h, Solver=>"CSDP")
g2 = sosPoly sol2

--Lax-Lax
R = QQ[a,b,c,d]
h = a^2 + b^2 + c^2 +d^2
f3=(a-b)*(a-c)*(a-d)*a+(b-a)*(b-c)*(b-d)*b+(c-a)*(c-b)*(c-d)*c+(d-a)*(d-b)*(d-c)*d+a*b*c*d
sol3 = solveSOS (f3*h, Solver=>"CSDP")
g3 = sosPoly sol3

--Harris polynomial
R = QQ[x,y,z]
(a,b,c,d,e) = (16,-36,20,57,-38)
f4 = a*( x^10 + y^10 + z^10)+ 
    b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
    c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) + 
    d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
    e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4)
h = x^2 + y^2 + z^2
sol4 = solveSOS (f4*h^3, Solver=>"CSDP")
g4 = sosPoly sol4

-- Scheiderer polynomial
-- h is SOS, but no rational decomposition exists
R = QQ[x,y,z]
f = (x^4 + x*y^3 + y^4 - 3*x^2*y*z - 4*x*y^2*z + 2*x^2*z^2 + x*z^3 + y*z^3 + z^4)
sol5 = solveSOS (f, Solver=>"CSDP")
g5 = sosPoly sol5


