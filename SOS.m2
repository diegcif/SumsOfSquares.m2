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

load "./SOS/SOS.m2"



TEST ///

 	 R = QQ[x,y];
         p = 4*x^4+y^4;
	 (g,d) = getSOS(p)
	 assert( p  == sumSOS(g,d) )
	 ///
	    