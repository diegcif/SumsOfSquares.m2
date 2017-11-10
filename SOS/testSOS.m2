loadPackage "SOS"
R = QQ[x,y]

-- Works well
f = 2*x^4+5*y^4-2*x^2*y^2+2*x^3*y;
(g,d) = getSOS f
sumSOS(g,d)-f

-- Does not work
f = x^4+2*x*y-x+y^4
(g,d) = getSOS f
sumSOS(g,d)-f
