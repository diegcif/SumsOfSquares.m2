newPackage(
        "SOS",
        Version => "1.5", 
        Date => "January 10, 2007",
        Authors => {{Name => "Helfried Peyrl", 
                  Email => "peyrl@control.ee.ethz.ch"},
	     {Name => "Pablo A. Parrilo",      
		  Email => "parrilo@mit.edu"}},
        Headline => "Sum of Squares Package",
        DebuggingMode => true,
	AuxiliaryFiles => true
        )

load "./SOS/SOSp.m2"



TEST ///
 	 R = QQ[x,y];
         p = 4*x^4+y^4;
	 (g,d) = getSOS(p)
	 assert( p  == sumSOS(g,d) )


         R = QQ[x,y,z,w];
	 p = 2*x^4 + x^2*y^2 + y^4 - 4*x^2*z - 4*x*y*z - 2*y^2*w + y^2 - 2*y*z + 8*z^2 - 2*z*w + 2*w^2;
	 (g,d) = getSOS p
	 assert( p == sumSOS (g,d) )

	 ///
	    