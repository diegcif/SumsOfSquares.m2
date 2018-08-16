needsPackage( "SOS", Configuration=>{"CSDPexec"=>"SDPsolvers/csdp","SDPAexec"=>"SDPsolvers/sdpa"} )

-- Bug 1: sdp is feasible, but csdp fails
-- Solution: set perturbobj=0
C = matrix {{8, 0, -4}, {0, 2, -3.14159265358979}, {-4, -3.14159265358979, 10.8696044010894}}
A1 = matrix {{.707106781186547, 0, 0}, {0, 0, 0}, {0, 0, 0}}
b = matrix {{-.707107}}
(y,X,Z) = solveSDP(C,Ai,b,Solver=>"CSDP")

-- Bug 2: this fails, but works if we round the entries
C = matrix{{8,0,-4},{0,2,-1},{-4,-1,2}}
Ai = matrix{{.707106781186547,0,0},{0,0,0},{0,0,0}}
b = matrix{{-.707107}}
(y,X,Z) = solveSDP(C,Ai,b,Solver=>"CSDP")

