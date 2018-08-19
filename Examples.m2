needsPackage( "SOS", Configuration=>{"CSDPexec"=>"SDPsolvers/csdp","SDPAexec"=>"SDPsolvers/sdpa"} )

-- Motzkin polynomial
R = QQ[x,y,z]
h = x^2 + y^2 + z^2
f1 = nonnegativeForm("Motzkin",R)
sol1 = solveSOS (f1*h, Solver=>"CSDP")
g1 = sosPoly sol1

--Robinson
f2 = nonnegativeForm("Robinson",R)
sol2 = solveSOS (f2*h, Solver=>"CSDP")
g2 = sosPoly sol2

--Lax-Lax
R = QQ[a,b,c,d]
h = a^2 + b^2 + c^2 +d^2
f3 = nonnegativeForm("Lax-Lax",R)
sol3 = solveSOS (f3*h, Solver=>"CSDP")
g3 = sosPoly sol3

--Harris polynomial
R = QQ[x,y,z]
f4 = nonnegativeForm("Harris",R)
h = x^2 + y^2 + z^2
sol4 = solveSOS (f4*h^3, Solver=>"CSDP")
g4 = sosPoly sol4

-- Scheiderer polynomial
-- h is SOS, but no rational decomposition exists
R = QQ[x,y,z]
f5 = nonnegativeForm("Scheiderer",R)
sol5 = solveSOS (f5, Solver=>"CSDP")
g5 = sosPoly sol5


