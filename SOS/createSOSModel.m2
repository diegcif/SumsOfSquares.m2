-- (C,Ai,mon,A,b,GramIndex) = CreateSOSModel(f)
--
-- CreateSOSModel   creates kernel and image model of the Gram matrix 
--    representation of a polynomial f:
--
--       f = mon' X mon,     	       	    (*)
--    
--    where X is a symmetric matrix indexed by mon, a vector of 
--    monomials.
-- 
--    C and Ai are symmetric matrices such that X can we written as
--
--       X = C - sum_i (Ai_i * y_i),
--    
--    where y_i are free parameters.
-- 
--    The linear subspace of all X such that (*) holds is determined by
--   
--       A x = b,
--
--    where x_i = X_(GramIndex#i).
--
-- Authors: Pablo Parrilo, Helfried Peyrl
-- $Id: createSOSModel.m2,v 1.3 2013-01-19 14:36:09 hpeyrl Exp $

createSOSModel = args -> (
    -- needs "./choosemonp.m2";

     args = sequence args;
     f := args#0;
     p := {}; if #args >= 2 then p = args#1;
     
     -- Degree and number of variables
     n := numgens ring f;
     d := (first degree f)//2;
     -- Get list of monomials for SOS decomposition
     (lmf,lm) := choosemonp (f,p);
             
     Hm := hashTable apply(lm, toList(1..#lm), identity);
     HHm := combine(Hm,Hm,times,(j,k)-> if j>=k then (j,k) else () , join );
     HHHm := applyValues(HHm,k->pack(2,k));
     -- This is a hash table that maps monomials into pairs of indices
     -- Writes the matrix, in sparse form
     ndim := #lm; -- binomial(n+d,d);
     mdim := #HHHm; --binomial(n+2*d,2*d);

     -- A hash table with the coefficients (should be an easier way)
     cf := coefficients f ;
     Hf := hashTable transpose {flatten entries cf#0, flatten entries cf#1};
     K := keys HHHm;

     -- Linear constraints: b     
     b := transpose matrix{apply (K, k-> (if Hf#?k then substitute(Hf#k,QQ) else 0))};

     -- Linear constraints: Ai, Bi
     Ah := new MutableHashTable;
     Bh := new MutableHashTable;
     LinSpaceDim := floor(ndim^2/2+ndim/2);
     LinSpaceIndex := hashTable apply (flatten values HHHm, toList(0..LinSpaceDim-1),identity);
     GramIndex := hashTable apply(values LinSpaceIndex, keys LinSpaceIndex, (i,j) -> (i,j));
     for k from 0 to #K-1 do (
	  -- Set constraints for monomial K#k 
     	  PairsEntries := toList HHHm#(K#k) ;
	  scan(PairsEntries, p -> (
	    	    if p_0 == p_1 then Ah#(k,LinSpaceIndex#p)=1_QQ else Ah#(k,LinSpaceIndex#p)=2_QQ;)
       	       );
	   -- Consider search-parameters:
	  for i from 0 to #p-1 do (
		    mp := K#k*p_i;
		    if Hf#?mp then Bh#(k,i) = -leadCoefficient Hf#mp;
		    );
	  );
   
     A := map(QQ^#K,QQ^(LinSpaceDim),(i,j) -> if Ah#?(i,j) then Ah#(i,j) else 0);
     if #p!=0 then B := map(QQ^#K,QQ^#p,(i,j) -> if Bh#?(i,j) then Bh#(i,j) else 0)  else B = ();
	       	       
     -- compute the C matrix
     c := b//A;
     C := map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then c_(LinSpaceIndex#{i+1,j+1},0) 
	  else c_(LinSpaceIndex#{j+1,i+1},0));
     -- compute the B_i matrices
     if #p!=0 then (
     	  bi := -B//A;
     	  Bi := apply(0..#p-1, k->
     	       map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then bi_(LinSpaceIndex#{i+1,j+1},k)
	       	    else bi_(LinSpaceIndex#{j+1,i+1},k)));
	  ) else Bi = ();
     -- compute the A_i matrices     
     v := - generators kernel A;
     
     Ai := apply(0..(rank v) - 1,k ->
	   map(QQ^ndim,QQ^ndim, (i,j) -> if i>=j then v_(LinSpaceIndex#{i+1,j+1},k) 
	       else v_(LinSpaceIndex#{j+1,i+1},k))); 
          
     (C,Ai,Bi,A,B,b,transpose matrix {lm},GramIndex,LinSpaceIndex)
     )
