needsPackage ("SOS" , Configuration => { "CSDPexec" => "CSDP/csdp"})

R = QQ[x,y,z]
(a,b,c,d,e) = (16,-36,20,57,-38);
f = a*( x^10 + y^10 + z^10)+ 
    b*( x^8* y^2 + x^2* y^8 + x^8* z^2 + x^2* z^8 + y^8* z^2 + y^2* z^8 ) +
    c*( x^6* y^4 + x^4* y^6 + x^6* z^4 + x^4* z^6 + y^6* z^4 + y^4* z^6 ) + 
    d*( x^6* y^2* z^2 + x^2* y^6* z^2 + x^2* y^2* z^6) +
    e*( x^4* y^4* z^2 + x^4* y^2* z^4 + x^2* y^4* z^4);
(nums, dens) = sosdecTernary(f, Solver=>"CSDP");
-- nums lie in different rings
product nums
