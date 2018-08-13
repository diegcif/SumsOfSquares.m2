needsPackage( "SOS", Configuration=>{"CSDPexec"=>"CSDP/csdp","SDPAexec"=>"CSDP/sdpa"} )

-- Motzkin polynomial
R = QQ[x,y,z]
h = x^2 + y^2 + z^2
f1 = x^4*y^2 + x^2*y^4 + z^6 - 3*x^2 *y^2 * z^2 --Motzkin
(mon1,Q1,X1,tval) = solveSOS (f1*h, Solver=>"CSDP")
g1 = sosPoly (mon1,Q1)

--Robinson
f2 = x^6 + y^6 + z^6 - (x^4*y^2 + x^2*y^4 + x^4*z^2 + x^2*z^4 + y^4*z^2 + y^2*z^4) + 3*x^2*y^2*z^2 --Robinson
(mon2,Q2,X2,tval) = solveSOS (f2*h, Solver=>"CSDP")
g2 = sosPoly (mon2,Q2)

--Lax-Lax
R = QQ[a,b,c,d]
h = a^2 + b^2 + c^2 +d^2
f3=(a-b)*(a-c)*(a-d)*a+(b-a)*(b-c)*(b-d)*b+(c-a)*(c-b)*(c-d)*c+(d-a)*(d-b)*(d-c)*d+a*b*c*d
(mon3,Q3,X3,tval) = solveSOS (f3*h, Solver=>"CSDP")
g3 = sosPoly (mon3,Q3)

--Harris polynomial
R = QQ[x,y,z]
(a,b,c,d,e) = (16,-36,20,57,-38)
f4 = a*( x^10 + y^10 + z^10)+ 
    b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
    c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) + 
    d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
    e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4)
h = x^2 + y^2 + z^2
(mon4,Q4,X4,tval) = solveSOS (f4*h^3, Solver=>"CSDP")
g4 = sosPoly (mon4,Q4)
