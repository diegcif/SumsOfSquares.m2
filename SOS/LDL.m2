-- (L,D,P,err) = ldl(A)
-- 
-- LDL   LDL' factorization of a positive semidefinte matrix.
--    (L,D,P,err) = ldl(A) returns a lower triangular matrix L with ones in
--    the diagonal, a diagonal matrix D, and a permutation matrix P such that
--    L'*D*L = P'*A*P. If A is positive semidefinite err=0, else err!=0.
--  
--    References: 
--    Gene Golub and Charles van Loan: Matrix Computations, Johns Hopkins
--    series in the Mathematical Science, 2 ed., pp. 133-148,
--    Baltimore Maryland, 1989.
--  
--    Author: Helfried Peyrl
--    $Id: LDL.m2,v 1.3 2013-01-19 14:39:21 hpeyrl Exp $


export {"ldl"}

ldl = args -> (
     args = sequence args;
     A := promote (args#0,QQ);
     if transpose A != A then error("Matrix must be symmetric.");
      
     n := #entries A;
     Ah := new MutableHashTable; map (QQ^n,QQ^n,(i,j)->Ah#(i,j) = A_(i,j));
     v := new MutableList from toList apply(0..n-1,i->0_QQ);
     d := new MutableList from toList apply(0..n-1,i->0_QQ);
     piv := new MutableList from toList(0..n-1);
     err := 0;
     
     for k from 0 to n-1 do (
	  q := maxPosition apply(k..n-1, i->Ah#(i,i)); q = q + k;
	  -- Symmetric Matrix Permutation:
	  tmp := piv#q; piv#q = piv#k; piv#k = tmp;
	  scan(0..n-1, i-> (tmp := Ah#(i,q); Ah#(i,q) = Ah#(i,k); Ah#(i,k) = tmp;));
	  scan(0..n-1, i-> (tmp := Ah#(q,i); Ah#(q,i) = Ah#(k,i); Ah#(k,i) = tmp;));
	       
	  --  positive semidefinite?
	  if Ah#(k,k) < 0 then (err = k+1; break;);
	  if (Ah#(k,k)==0) and (number(apply(0..n-1,i->Ah#(i,k)),f->f!=0)!=0) then (
	       err = k+1; break;);
	  
	  -- Perform LDL factorization step:
	  if Ah#(k,k) > 0 then (
      	       scan(0..k-1, i -> v#i = Ah#(k,i)*Ah#(i,i));	  
	       Ah#(k,k) = Ah#(k,k) - sum apply(toList(0..k-1), i -> Ah#(k,i)*v#i);
	       if Ah#(k,k) < 0 then (err = k+1; break;);
	       if Ah#(k,k) > 0 then (
		    scan(k+1..n-1, i ->
			 (Ah#(i,k) = (Ah#(i,k)-sum apply(toList(0..k-1),j->Ah#(i,j)*v#j)) 
			 / Ah#(k,k);))
	       	    );
	  );
     );

     A = map(QQ^n,QQ^n,(i,j)-> if i>j then Ah#(i,j) else if i==j then 1_QQ else 0_QQ);
     D := map(QQ^n,QQ^n,(i,j)->if i==j then Ah#(i,j) else 0_QQ);
     P := submatrix(id_(QQ^n),toList piv);
     (A,D,P,err)     
)


beginDocumentation()
document {
        Key => {ldl},
        Headline => "LDL' factorization of a positive semidefinite matrix",
	"If ", TT "A", " is a positive semidefinite matrix, ", EM "ldl", " returns a lower 
	triangular matrix ", TT "L", " with ones in the diagonal, a diagonal matrix ",
	TT "D", " and a permutation matrix ", TT "P", " such that ", TT "L'*D*L = P'*A*P.",
        Usage => "(L,D,P,err) = ldl A",
        Inputs => { "A" => Matrix => {"over ", TT "QQ", " or ", TT "ZZ." } },
        Outputs => { "L" => Matrix => {"a lower triangular matrix over ", TT "QQ."},
	"D" => Matrix => {"a diagonal matrix over ", TT "QQ."},
	"P" => Matrix => {"a permuation matrix over ", TT "QQ."},
	"err" => ZZ => {"which is 0 when the factorization was successful."}},
        --SourceCode => {ldl},
        EXAMPLE lines ///
          A = matrix {{5,3,5},{3,2,4},{5,4,10}}
	  (L,D,P,err) = ldl(A)
	  ///
        }
   
TEST ///

--  Simple example
    A = matrix {{5,3,5},{3,2,4},{5,4,10}}
    (L,D,P,err) = ldl(A)
    assert(L*D*transpose L == transpose P * A * P)
    
--  Random low-rank matrix
    V = random(QQ^12,QQ^8)
    A = V * transpose V 
    (L,D,P,err) = ldl(A)
    assert(L*D*transpose L == transpose P * A * P)

///

