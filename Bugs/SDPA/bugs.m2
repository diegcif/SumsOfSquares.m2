needsPackage( "SOS", Configuration=>{"CSDPexec"=>"CSDP/csdp","SDPAexec"=>"CSDP/sdpa"} )

-- Bug 1: sdp is feasible, bug sdpa doesn't believe it
C = matrix(RR, {{1, 0, 0}, {0, 0, 0}, {0, 0, -1}} )
A1 = matrix(RR, {{0, 0, 1/2}, {0, -1, 0}, {1/2, 0, 0}} )
A2 = matrix(RR, {{0, -1/4, 0}, {-1/4, 0, 0}, {0, 0, -1}} )
A = (A1,A2)
b = matrix {{0}, {0}}
(y,X,Z) = solveSDP(C,A,b,Solver=>"SDPA")

-- Bug 1a: optimal value should be negative
A0 = matrix(RR, {{1,0,0},{0,1,0},{0,0,1}} )
b' = matrix(RR, {{-1}, {0}, {0}} )
A' = (A0,A1,A2)
(y',X',Z') = solveSDP(C,A',b',Solver=>"SDPA")
opt' = transpose b'*y'
