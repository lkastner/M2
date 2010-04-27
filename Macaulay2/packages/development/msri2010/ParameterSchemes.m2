-- -*- coding: utf-8 -*-
newPackage(
	"ParameterSchemes",
    	Version => "0.02", 
    	Date => "April 26,2010",
    	Authors => {
	     {Name => "Mike Stillman", Email => "mike@math.cornell.edu", HomePage => "http://www.math.cornell.edu/~mike/"},
	     {Name => "Kristine Fisher", Email => "kejones@math.cornell.edu"}},
    	Headline => "a Macaulay2 package for local equations of Hilbert and other parameter schemes",
    	DebuggingMode => true
    	)

export { 
     standardMonomials,
     findWeightVector,
     findHeft,
     findEliminants,
     groebnerFamily,
     groebnerBasin,
     Minimalize,
     linearPart,
     AllStandard,
     smallerMonomials
     }

smallerMonomials = method()
smallerMonomials(Ideal) := M ->
    ( error "not implemented yet")

standardMonomials = method()
standardMonomials(ZZ,Ideal) :=
standardMonomials(List,Ideal) := (d,M) -> (
     (terms sum flatten entries basis(d,comodule M))/leadMonomial
     )

standardMonomials(Ideal) := (M) -> (
     apply(M_*, f -> standardMonomials(degree f, M))
     )
standardMonomials(List) := (L) -> standardMonomials ideal L

findHeft = (degs) -> (
     A := transpose matrix degs;
     needsPackage "FourierMotzkin";
     B := ((value getGlobalSymbol "fourierMotzkin") A)#0;
     r := rank source B;
     f := (matrix{toList(r:-1)} * transpose B);
     if f == 0 then return;
     heft := first entries f;
     g := gcd heft;
     if g > 1 then heft = apply(heft, h -> h // g);
     minheft := min heft;
     if minheft <= 0 then heft = apply(heft, a -> a - minheft + 1);
     heftvals := apply(degs, d -> sum apply(d, heft, times));
     if not all(heftvals, d -> d > 0) then return;
     (heft, heftvals)
     )

findWeightVector = method ()
findWeightVector(Ideal, List) := (M, L) -> (
     --input: a monomial ideal M and a list of lists of standard monomials L
     --output: a weight vector w (a list of integers of length the # of variables
       -- of the ring of M) 
     --w places the listed generators of M greater than the corresponding standard monomials)
     R := ring M;
     kk := coefficientRing R;
     nv := sum apply(L, s -> #s);
     Mlist := M_*;
     D := flatten apply(#Mlist, i -> (
	       m := Mlist_i;
	       apply(L#i, s -> (
			 first exponents m - first exponents s))));
     findHeft D
     )
     
findEliminants = (M,J) -> (
-- input: a monomial ideal M, and the parameter family J
-- output: an ordered pair of lists of indices; the first entry is the variables that
      --will be eliminated, the second entry is the variables that will not be eliminated    
     S := ring J;
     SP1 := (gens J) * sub(syz gens M, S);
     SP2 := SP1 % sub(M,S);
     L1 := trim ideal flatten last coefficients SP2;
     L1 = sub(L1,coefficientRing S);
     -- set1 will be the variables which are lead terms here:
     set1 := sort ((flatten entries leadTerm L1)/index);
     set2 := sort toList(set toList(0..numgens coefficientRing S - 1) - set set1);
     (set1, set2)
     )     

groebnerFamily = method (Options => {AllStandard => true, Variable => global t})
groebnerFamily (Ideal) := opts -> (M) -> ( 
     L := if opts.AllStandard then standardMonomials M else smallerMonomials M;
     groebnerFamily(M, L, opts)
     )
groebnerFamily (Ideal, List) := opts -> (M, L) -> (
     (w, wvals) := findWeightVector(M,L);
     R := ring M;
     kk := coefficientRing R;
     nv := #wvals;   
     R1 := kk[t_1..t_nv, Degrees => wvals];
     U := R1 (monoid[gens R, Join=>false, Degrees => w]);
     phi := map(U,R1);
     lastv := -1;
     Mlist := M_*;
     elems := apply(#Mlist, i -> (
	       m := Mlist_i;
	       substitute(m,U) + sum apply(L#i, p -> (
			 lastv = lastv + 1;
			 phi(R1_lastv) * substitute(p,U)))));
     J := ideal elems;
     --determining which parameters can be eliminated
     (indices1, indices2) := findEliminants(M,J);  
     gens1 := apply(indices1, i -> (wvals#i, i));
     indices1 = apply(rsort gens1, (h,i) -> i);
     -- Now make the actual coeff ring:
     R2 := kk[(gens R1)_indices1, (gens R1)_indices2, 
	  Degrees => join(wvals_indices1,wvals_indices2),
	  MonomialOrder => {Lex => (#indices1), #indices2}];
     U = R2 (monoid[gens R, Join=>false, Degrees => w]);
     f1 := map(R2,R1);
     f2 := map(U, ring J, vars U | f1.matrix);
     J = f2 J;
     J
     )

linearPart = (f) -> sum select (terms f, t->(
	  sum first exponents t === 1))

minimalizeFamily = method()
minimalizeFamily (Ideal, Ideal) := (J,H) -> (
     R := ring H;
     ringgens := flatten entries selectInSubring (1,vars ring H);
     newR := (coefficientRing R)[ringgens, Degrees => ringgens/degree];
     newringJ := newR( monoid[gens ring J, Join => false, Degrees => degrees(ring J)]);
     L := matrix {apply (H_*, f->linearPart f)};
     C := getChangeMatrix gb (L, ChangeMatrix => true);
     time elimH := trim ideal ((gens H)*C); 
     D := ideal leadTerm ((gens H)*C);
     assert( # unique D_* == numgens D);
     newH := compress ((gens H)% elimH);
     time newH = trim ideal substitute(newH, newR);
     elimH' := promote (elimH, ring J);
     newJ := ideal ((gens J) % elimH');
     (newJ, newH)
     )
     
groebnerBasin = method(Options => {AllStandard => true, Minimalize => true})
groebnerBasin Ideal := opts -> (M) -> (
    --input: a monomial ideal M
    -- output: (J, H) two ideals: J is the family and H is the
       -- ideal on the coefficients
    J := groebnerFamily (M, AllStandard => opts.AllStandard);
    R := ring J;
    G := forceGB gens J;
    syzM := substitute(syz gens M, R);
    eq := compress((gens J * syzM) % G);
    (mons,eqns) := coefficients(eq); 
    H := ideal lift(eqns,coefficientRing R);
    if not opts.Minimalize then return (J,H);
    minimalizeFamily (J,H)
    )   

---------------
--Documentation
----------------
beginDocumentation()
------------
--Front Page
-------------
doc///
  Key
    ParameterSchemes
  Headline
    a package for working with parameter schemes
  Description
    Text
      {\em ParameterSchemes} package is designed for computing homogeneous strata and parameter families for monomial ideals in polynomial rings.
      In certain instances homogeneous strata can be used to compute local coordinates on Hilbert Schemes.
///

----------
--Functions
-------------

doc ///
  Key
    standardMonomials
    (standardMonomials, Ideal)
    (standardMonomials, List)
    (standardMonomials, List, Ideal)
    (standardMonomials, ZZ, Ideal)
  Headline
     computes standard monomials  
  Usage
    L = standardMonomials M
      = standardMonomials G 
      = standardMonomials(d, M)
  Inputs
    M : Ideal
      M should be a monomial ideal
    G : List
      G is a list of monomials that generate M as an ideal
    d : List
      d should be an integer
  Outputs
    L : List
      L is a list of lists of standard monomials for the generators of M
  Description
    Text
      A monomial $$m$$ is standard with respect to a monomial ideal $$M$$ and a generator $$g$$ of $$M$$ if $$m$$ 
      is of the same degree as $$g$$ but is not an element of $$M$$. 
      
      Inputting an ideal $M$ returns the standard monomials of each of the cached generators of the ideal.
    Example
      R = ZZ/32003[a..d];
      M = ideal (a^2, a*b, b^3, a*c);
      L1 = standardMonomials M  
    Text
      Inputting a list of monomials $G$ returns the standard monomials for each member of the list, with respect to the ideal
      generated by $G$
    Example
      G = {a^2, a*b, b^3, a*c};
      L2 = standardMonomials G
    Text
      Inputting an integer $d$ and an ideal gives the standard monomials for the specified ideal in degree $d$.
    Example
      standardMonomials( 2, M)
  SeeAlso
    smallerMonomials
///

doc ///
  Key
    findWeightVector
    (findWeightVector, Ideal, List)
  Headline
    returns a weight vector
  Usage
    (w, h) = findWeightVector (M, L)
  Inputs
    M : Ideal
      M should be a monomial ideal
    L : List
      a list of lists of standard monomials
  Outputs 
    w : List
      w is a weight vector that places the specified generators of M greater than the corresponding
      standard monomials
    h : List
      h is a list of values for w dotted with the difference of the exponent of the each standard monomial for each 
      generator and the corresponding generator, in the order they are listed in L    
  Description
    Text
      In the first entry, this command returns a weight vector associated to a monomial order that places the generators 
      of a monomial ideal $M$ ahead of standard monomials of the same degree.  The second entry is a list of values for 
      the weight vector dotted with the difference of the exponent of each standard monomial for each generator and the 
      corresponding generator.
      
      Note that the desired term ordering, and hence weight vector, may not exist.
    Example
      R = ZZ/32003[a,b,c, d];
      M = ideal (a^2, a*b, b^2);
      L = standardMonomials M;
      findWeightVector(M,L)
    Text
      Note that the first generator listed for $M$ is $a^2$, and the first corresponding standard monomial is $a*c$.  The difference
      of these two monomials exponent vectors is $(1,0,-1,0)$.  This vector dotted with the weight vector $(2,2,1,1)$ gives the value $1$, 
      which is the first value in the second list.
      
      This command is used in the @TO groebnerBasin@ and @TO groebnerFamily@ routines.     
  SeeAlso
    findHeft
    groebnerBasin
    groebnerFamily
///

doc ///
  Key
    findHeft
  Headline
    computes heft vector
  Usage
    (w, h) = findHeft(D)
  Inputs
    D : List
      D is a list of the exponent vectors of standard monomials subtracted from the monomials with respect
      to which they are standard - one should have a set of generators for a monomial ideal in mind
  Outputs
    w : List
      w is a weight vector that places the specified monomials greater than the corresponding
      standard monomials
    h : List
      h is a list of values for w dotted with each element in D 
  Description
    Text
      This command returns a weight vector whose associated term ordering places desired monomials greater 
      than corresponding standard monomials.  It also returns a list of values for the weight vector dotted with each element
      in the input list.
      
      Note that such a weight vector may not exist for certain inputs.
    Example
      R = ZZ/32003[a,b]
      M = ideal a^2
      standardMonomials M
      D = {{1,-1}, {2,-2}}
      findHeft D
    Text
      Note that $(2,1)$ dotted with $(1,-1)$ gives $1$, and $(2,1)$ dotted with $(2,-2)$ gives $2$.  
      
      This command is primarily designed to enable the routines for @TO groebnerBasin@ and @TO groebnerFamily@. 
      The @TO findWeightVector@ command computes the same weight vector, but allows for input of an ideal and standard monmials
      instead of the list $D$.
  SeeAlso
    groebnerBasin
    groebnerFamily
///

doc ///
  Key
    findEliminants
  Headline
    determines a maximal set of eliminable parameters 
  Usage
    (A, B) = findEliminants(M, J)
  Inputs
    M : Ideal
      M should be a monomial ideal
    J : Ideal
      J is the parameter family of M
  Outputs
    A : List
      the indices of the parameters that will be eliminated
    B : List
      the indices of the parameters that will not be eliminated
  Description
    Text
      When computing the parameter family of a monomial ideal, a new parameter is introduced for each standard monomial for each
      generator of M.  Often when looking at the associated stratum many of these parameters are elimnable, and the stratum can be recognized 
      isomorphically in a smaller polynomial ring.  
    Example
      R = ZZ/32003[a,b,c]
      M = ideal (a^2, a*b)
      J = groebnerFamily(M)
      findEliminants (M, J)
      (groebnerBasin(M))_1
    Text
      In this case, the parameters t_4, t_8, t_3, and t_1 will be eliminated.  Note that the output of {\tt findEliminants} lists eliminants in terms 
      of their position in the list of variables in the polynomial ring associated to the parameter family.  The eliminants are NOT listed with respect
      to their subscripts.
  SeeAlso
    standardMonomials
    groebnerFamily
///
doc ///
  Key
    groebnerFamily
    (groebnerFamily, Ideal)
    (groebnerFamily, Ideal, List)
    [groebnerFamily, AllStandard]
  Headline
    computes families of ideals with a specified initial ideal 
  Usage
    J = groebnerFamily(M)
    J = groebnerFamily(M, L)
  Inputs
    M : Ideal
      M should be a monomial ideal
    L : List
      L is a list of lists of standard monomials or smaller standard monomials for the generators of M
      --Not yet functional for smaller standard monomials ie smallerMonomials!!!!!!
  Outputs
    J : Ideal
      J is the groebner family, an ideal in the polynomial ring over the original variables and the parameters
  Description
    Text
      Given a monomial ideal $M$ in a polynomial ring $R$, this computes the paramter families of homogeneous ideals where 
      $M$ could be their initial ideal. 
      These families are obtained from either the standard monomials to the generators of $M$, or the standard monomials smaller than
      the generators of $M$ but of the same degree as these generators. 
      In the former case we obtain a family of all ideals where $M$ could be their initial ideal.  
      In the latter case, we obtain such a family with respect 
      to a given term order.  This second functionality is not yet implemented.    
    Example
      R = ZZ/32003[a,b,c,d]
      M = ideal (a^2, a*b, b^2)
      J = groebnerFamily M
    Text
      Here, $J$ is the family of homogeneous ideals having $M$ as their initial ideal, under some term order.
      
      The optional argument @TO AllStandard@ is boolean, taking the value $true$ to compute the family of all homogeneous ideals
      with a given initial ideal and the value $false$ to compute the family with respect to a given order.  The default value for 
      this argument is $true$.
      
      If $L$, the list of standard monomials for $M$, has already been computed, this function will take it as an additional argument, so that it
      is not re-computed.      
    Example
      L = standardMonomials M
      J2 = groebnerFamily (M, L)
    Text  
      Note that $J$ and $J2$ are the same family.
  SeeAlso
    groebnerBasin
    smallerMonomials
    standardMonomials
///

doc ///
  Key
    groebnerBasin
    (groebnerBasin, Ideal)
    [groebnerBasin, Minimalize]
    [groebnerBasin, AllStandard]
  Headline
    returns homogeneous Groebner family and Groebner basin (stratum)
  Usage
    (J, H) = groebnerBasin(M)
  Inputs
    M : Ideal
      M should be a monomial ideal
  Outputs
    J : Ideal
    	 J is the Groebner family (parameter family)
    H : Ideal
    	 H is the Groebner basin (ideal on the coefficients, or homogeneous stratum)
  Description
    Text
      Given a monomial ideal $M$, this command returns a family of ideals having $M$ as an initial ideal, and conditions on the parameters so that the
      family is flat.  If the optional input @TO AllStandard@ is specified as $true$, the family of all homogeneous ideals having $M$ as an initial ideal is 
      computed, and if it is specified as $false$ the family of all homogenous ideals having $M$ as an initial ideal with respect to the given term order is 
      computed.  The second functionality is not yet implemented.  
    Example
      R = ZZ/32003[a,b,c]
      M = ideal (a^2, a*b, b^2)
      (J, H) = groebnerBasin M
      H
    Text
      In this example, J is the universal family, and H is the ideal giving the conditions on the parameters.  In general, several of the parameters are
      unnecesary. Note that H is an ideal in a ring with far fewer parameters.  This is because a maximal set of eliminable parameters from the original
      ideal of conditions on parameters have been eliminated.  If the full ideal in the polynomial ring over all the parameters is desired, set the optional 
      argument $Minimalize$ to false.
    Example
      (J1, H1) = groebnerBasin(M, Minimalize =>false)
      netList H1_*
    Text
      Notice that the parameters $t_3$, $t_6$, and $t_9$ are clearly eliminable.          
  SeeAlso
    findEliminants
    groebnerFamily
///

doc ///
  Key
    linearPart
  Headline
    returns linear part of a polynomial
  Usage
    l = linearPart(f)
  Inputs
    f : Thing
      f is a polynomial
  Outputs
    l : Thing
      l is the sum of the terms of f that have exponent equal to 1
  Description
    Text
      Use this command to return the sum of the terms of a polynomial with exponent 1.
    Example
      R = ZZ/32003[a,b,c];
      f = a^2 + 3*b + 5*a*b*c + 2*a +b^2;
      linearPart f
    Text
      If we have a polynomial ring over another polynomial ring, the command only looks at the exponents on the new set of variables.
    Example
      S = ZZ/32003[a,b,c][x,y,z];
      g = a*b*x + 3*y + a + b^2
      linearPart g
    Text 
      Note that in the previous example, {\tt linearPart g} does not include {\tt a} but does include {\tt a*b*x}
      
      This command is used in the @TO findEliminants@ routine to find eliminable parameters in the computation of 
      parameter families and strata.
  SeeAlso    
    findEliminants 
///

doc ///
  Key
    smallerMonomials
    (smallerMonomials, Ideal)
  Headline
    returns the standard monomials smaller than given monomials in given degrees
  Description
    Text
      This function is not yet implemented
  SeeAlso
    standardMonomials
///
---------
--Symbols/optional arguments
-------------
doc ///
  Key
    AllStandard
  Headline 
    boolean option for determining the use of all standard or smaller standard monomials
  Description
    Text
      This is an optional input for the @TO groebnerBasin@ and @TO groebnerFamily@ functions.  It takes values $true$ and $false$.  When 
      assigned the value $true$, the functions are computed with respect to all standard monomials for the specified generators.  When assigned the 
      value $false$, the functions are computed with respect to smaller monomials of the same degree as the specified generators but which do not lie in 
      the ideal.  
--the false value is not yet functional!!! 
  SeeAlso
    groebnerBasin
    groebnerFamily
    smallerMonomials
    standardMonomials
///

doc ///
  Key
    Minimalize
  Headline
    boolean option for determning whether excess parameters will be elimnated
  Description
    Text
      This an optional input for the @TO groebnerBasin@ function.  It takes values $true$ and $false$.  When assigned the value
      $true$, elimnable parameters are eliminated to obtain the groebner basin / stratum in a smaller ring.  When assigned the 
      value $false$, no parameters are eliminated, and the groebner basin / stratum is comupted in a ring containing all the parameters.
      The default input is $true$.
  SeeAlso
    findEliminants
    groebnerBasin
///

--------------------
--Tests
--------------------

--Test 0 standard monomials
TEST///
R = ZZ/32003[a..d]
M = ideal (a^2, a*b, b^3, a*c)
L = standardMonomials M
ans = {{b^2, b*c, c^2, a*d, b*d, c*d, d^2}, {b^2, b*c, c^2, a*d, b*d, c*d, d^2}, {b^2*c, b*c^2, c^3, b^2*d, b*c*d, c^2*d, a*d^2, b*d^2, c*d^2, d^3}, {b^2, b*c, c^2, a*d, b*d, c*d, d^2}}
assert (L == ans)
///
end














restart 
loadPackage "ParameterSchemes"

TEST ///
R = ZZ/32003[a..d]
M = ideal (a^2, a*b, b^3, a*c)
L = standardMonomials M
ans = {{b^2, b*c, c^2, a*d, b*d, c*d, d^2}, {b^2, b*c, c^2, a*d, b*d, c*d, d^2}, {b^2*c, b*c^2, c^3, b^2*d, b*c*d, c^2*d, a*d^2, b*d^2, c*d^2, d^3}, {b^2, b*c, c^2, a*d, b*d, c*d, d^2}}
assert (L == ans)
(w, hvals) = findWeightVector(M, L)
J = groebnerFamily(M,L)
netList J_*
(nJ, nH) = groebnerBasin M
netList (ideal gens gb nH)_*
netList nJ_*
/// 

TEST ///
R = ZZ/32003[a,b,c]
M = ideal (a^2, a*b, b^3)
J = groebnerFamily M
netList J_*
(J, H) = groebnerBasin M
netList (ideal gens gb H)_*
netList J_*
///

TEST ///
R = ZZ/32003[a,b,c, d]
M = ideal (a^2, a*b, b^2)
J = groebnerFamily M
netList J_*
(J, H) = groebnerBasin M
assert isHomogeneous J
assert isHomogeneous H
netList (ideal gens gb H)_*
///

TEST ///
R = ZZ/32003[a,b,c, d]
M = ideal (a,b^4,b^3*c)
M' = truncate (3, M)
J = groebnerFamily M'
(J, H) = groebnerBasin M'
netList J_*
///

TEST ///
--Hilbert polynomial 4z
R = ZZ/32003[a,b,c, d]
M = ideal (a^2, a*b, a*c, b^5, b^4*c)
M' = truncate (4, M)
J = groebnerFamily M'
netList J_*
///
export { 
     findWtVec,
     smallerMonomials, 
     standardMonomials,
     parameterFamily, 
     parameterIdeal, 
     parameterRing,
     pruneParameterScheme, 
     groebnerScheme
     }


smallerMonomials = method()
smallerMonomials(Ideal,RingElement) := (M,f) -> (
     -- input: a polynomial in a poly ring R
     -- output: an ordered list of monomials of R less than f, but of the same
     --   degree as (the leadterm of) f.
     R := ring M;
     d := degree f;
     m := flatten entries basis(d,coker gens M);
     m = f + sum m;
     b := apply(listForm m, t -> R_(first t));
     x := position(b, g -> g == f);
     drop(b,x+1))

smallerMonomials(Ideal,RingElement,ZZ) := (M,f,dummy) -> (
     -- input: a polynomial in a poly ring R
     -- output: an ordered list of monomials of R less than f, but of the same
     --   degree as (the leadterm of) f.
     d := degree f;
     m := flatten entries basis(d,coker gens M);
     select(m, m0 -> m0 < f))

smallerMonomials(Ideal) := (M) -> (
     Mlist := flatten entries gens M;
     apply(Mlist, m -> smallerMonomials(M,m)))

smallerMonomials(List) := (L) -> (
     M := ideal L;
     apply(L, m -> smallerMonomials(M,m)))

standardMonomials = method()
standardMonomials(ZZ,Ideal) :=
standardMonomials(List,Ideal) := (d,M) -> (
     (terms sum flatten entries basis(d,comodule M))/leadMonomial
     )

standardMonomials(Ideal) := (I) -> (
     L := flatten entries generators I;
     apply(L, f -> standardMonomials(degree f, I))
     )
standardMonomials(List) := (L) -> (
     I := ideal L;
     apply(L, f -> standardMonomials(degree f, I))
     )

findWtVec = method(Options => {Standard => false});
findWtVec Ideal := opts -> (I) -> (
     M := I_*;
     L := if opts.Standard 
           then standardMonomials I 
	   else smallerMonomials I;
     W := flatten apply(#M, i -> (
	       e := first exponents M#i;
	       apply(L#i, t -> e - first exponents t)));
     W = transpose matrix W;
     (C,H) := fourierMotzkin W;
     w := flatten entries sum(numColumns C, i -> C_{i});
     maxw := max w;
     result := apply(#w, i -> maxw+1-w#i);
     minw := min flatten entries(matrix {result} * W);
     if minw > 0 then (
	  << "minimal weight is " << minw << endl;
	  result
	  )
     )

findWeightCone = method()
findWeightCone Ideal := (I) -> (
     -- I should be a monomial ideal
     L := flatten entries gens I;
     M := transpose matrix flatten apply(L, f -> (
	       b := basis(degree f, ring I);
	       b = compress (b % I);
	       apply(first entries b, g -> first exponents f - first exponents g)
	       ));
     M0 = M;
     (C,H) := fourierMotzkin M;
     M = fourierMotzkin(C,H);
     C = -C;
     w := flatten entries sum(numColumns C, i -> C_{i});
     minw := min w;
     result := apply(#w, i -> -minw+1+w#i);
     (result,C,M)
     )

findEliminants = (M,J) -> (
     S := ring J;
     SP1 := (gens J) * sub(syz gens M, S);
     SP2 := SP1 % sub(M,S);
     L1 := trim ideal flatten last coefficients SP2;
     L1 = sub(L1,coefficientRing S);
     -- set1 will be the variables which are lead terms here:
     set1 := sort ((flatten entries leadTerm L1)/index);
     set2 := sort toList(set toList(0..numgens coefficientRing S - 1) - set set1);
     << netList J_* << endl;
     << set1 << endl;
     << set2 << endl;
     (set1, set2)
     )

parameterRing = method()
parameterRing(Ideal,List,Symbol) := (M,L,t) -> (
     -- M is a monomial ideal
     -- L is a list of lists of monomials, #L is the
     --  number of generators of M.
     R := ring M;
     kk := coefficientRing R;
     nv := sum apply(L, s -> #s);
     Mlist := M_*;
     D := flatten apply(#Mlist, i -> (
	       m := Mlist_i;
	       apply(L#i, s -> (
			 first exponents m - first exponents s))));
     Dxvars := entries id_(ZZ^(numgens R));
     Dall := join(Dxvars, D);
     s := local s;
     u := local u;
     R0 := kk (monoid[u_1..u_(numgens R), s_1..s_nv, Degrees => Dall]);
     w := heft R0;
     heftvals := flatten entries(matrix D * transpose matrix{w});
     R1 := kk[t_1..t_nv, Degrees => D, Heft => w];
     U := R1 (monoid([gens R, Join=>false, Degrees => Dxvars, Heft => w]));
     phi := map(U,R1);
     lastv := -1;
     elems := apply(#Mlist, i -> (
	       m := Mlist_i;
	       substitute(m,U) + sum apply(L#i, p -> (
			 lastv = lastv + 1;
			 phi(R1_lastv) * substitute(p,U)))));
     J := ideal elems;
     -- Next step: determine the variables to eliminate
     (indices1, indices2) := findEliminants(M,J); -- 
     gens1 := apply(indices1, i -> (heftvals#i, i));
     indices1 = apply(rsort gens1, (h,i) -> i);
     -- Now make the actual coeff ring:
     R2 := kk[(gens R1)_indices1, (gens R1)_indices2, Degrees => join(D_indices1,D_indices2), MonomialOrder => {Lex => (#indices1), #indices2}, Heft => w];
     U := R2 (monoid([gens R, Join=>false, Degrees => Dxvars, Heft => w]));
     f1 := map(R2,R1);
     f2 := map(U, ring J, vars U | f1.matrix);
     J = f2 J;
     -- last step: create the relations on these coeffs:
     J
     --(parameterIdeal(M, J), J)
     )

parameterFamily = method(Options=>{Local=>false, Weights=>null})
parameterFamily(Ideal,List,Symbol) := opts -> (M,L,t) -> (
     -- M is a monomial ideal
     -- L is a list of lists of monomials, #L is the
     --  number of generators of M.
     local R1;
     local U;
     local phi;
     R := ring M;
     kk := coefficientRing R;
     nv := sum apply(L, s -> #s);
     Mlist := flatten entries gens M;
     dot := (a,b) -> sum for i from 0 to #a-1 list (a#i * b#i);
     degs := if opts.Weights =!= null
       then (
	  flatten apply(#Mlist, i -> (
		    m := Mlist_i;
		    flatten apply(L#i, s -> (
			      dot(opts.Weights,first exponents m - first exponents s)
			      ))
		    ))
	  );
     if degs =!= null then if min degs <= 0 then error "expected positive weight values";
     if opts.Local then (
         R1 = kk{t_1..t_nv}; -- removed MonomialSize=>8
	 U = kk[gens R,t_1..t_nv,MonomialOrder=>{
		   Weights=>splice{numgens R:0,nv:-1},numgens R,nv},
	           Global=>false];
	 phi = map(U,R1,toList(U_(numgens R) .. U_(nv - 1 + numgens R)));
         )
     else (
	 R1 = if opts.Weights === null then 
	          kk[t_1..t_nv] -- removed MonomialSize=>8
	      else kk[t_1..t_nv,Degrees=>degs,MonomialOrder=>{Weights=>degs}];
     	 U = R1 (monoid R);
	 U = newRing(U, Join=>false);
	 phi = map(U,R1);
     );
     lastv := -1;
     elems := apply(#Mlist, i -> (
	       m := Mlist_i;
	       substitute(m,U) + sum apply(L#i, p -> (
			 lastv = lastv + 1;
			 phi(R1_lastv) * substitute(p,U)))));
     ideal elems
     )

parameterIdeal = method()
parameterIdeal(Ideal,Ideal) := (M,family) -> (
     -- M is a monomial ideal in a polynomial ring
     -- family is the result of a call to 'parameterFamily'
     R := ring M;
     time G = forceGB gens family;
     time syzM = substitute(syz gens M, ring family);
     time eq = compress((gens family * syzM) % G);
     --time (mons,eqns) := toSequence coefficients(toList(0..(numgens R)-1), eq);
     time (mons,eqns) = coefficients(eq); -- , Variables=>apply(gens R, x -> substitute(x,ring eq)));
     ideal lift(eqns,coefficientRing ring eqns))

pruneParameterScheme = method()
pruneParameterScheme(Ideal,Ideal) := (J,F) -> (
     R := ring F;
     A := coefficientRing R;
     if ring J =!= A then error "expected(ideal in coeffring A, family in A[x])";
     time J1 := minimalPresentation J; -- minPressy J;
     map1 := J.cache.minimalPresentationMap;
     map2 := J.cache.minimalPresentationMapInv;
     B := ring J1;
     phi := map1; -- map: A --> B
     -- want the induced map from A[x] -> B[x]
     S := B (monoid R);
     phi' := map(S,R,vars S | substitute(phi.matrix,S));
     (J1, phi' F)
     )

preprune = (J,F) -> (
     -- asumption: J is homogeneous
     R' := (coefficientRing ring J) [ gens ring J ];
     J' := sub(J,R');
     J0 := trim ideal apply(J'_*, f -> part_1 f);
     pos := set ((flatten entries leadTerm J0)/index//sort);
     wt := toList apply(0..numgens ring J - 1, i -> if member(i,pos) then 1 else 0);
     others := sort toList((set toList(0..numgens ring J-1)) - pos);
     pos = sort toList pos;
     degs := flatten degrees ring J;
     S := (coefficientRing ring J)[(gens ring J)_pos, (gens ring J)_others, 
	  MonomialOrder=>{Weights=>join(wt_pos,wt_others), #pos, #others},
     	  Degrees => join(degs_pos,degs_others), 
	  MonomialSize => 8];
     -- Now we need to change the ring of F
     T := S (monoid ring F);
     idS := map(S, ring J);
     Fnew := (map(T, ring F, vars T | sub(idS.matrix,T))) F;
     (sub(J,S), Fnew)
     )

preprune = (J,F) -> (
     -- asumption: J is homogeneous
     R' := (coefficientRing ring J) [ gens ring J ];
     J' := sub(J,R');
     J0 := trim ideal apply(J'_*, f -> part_1 f);
     pos := set ((flatten entries leadTerm J0)/index//sort);
     wt := toList apply(0..numgens ring J - 1, i -> if member(i,pos) then 1 else 0);
     others := sort toList((set toList(0..numgens ring J-1)) - pos);
     pos = sort toList pos;
     degs := flatten degrees ring J;
     S := (coefficientRing ring J)[(gens ring J)_pos, (gens ring J)_others, 
	  MonomialOrder=>{Weights=>join(wt_pos,wt_others), #pos, #others},
     	  Degrees => join(degs_pos,degs_others), 
	  MonomialSize => 8];
     -- Now we need to change the ring of F
     T := S (monoid ring F);
     idS := map(S, ring J);
     Fnew := (map(T, ring F, vars T | sub(idS.matrix,T))) F;
     (sub(J,S), Fnew)
     )

groebnerScheme = method(Options=>{Minimize=>true, Weights=>null})
groebnerScheme Ideal := opts -> (I) -> (
     L1 := smallerMonomials I;
     F0 := parameterFamily(I,L1,symbol t,Weights=>opts.Weights);
     J0 := parameterIdeal(I,F0);
     if false and isHomogeneous J0 then (
	  -- We change the ring so that the lead terms of the linear parts
	  -- will occur as lead terms of J0.
	  (J0,F0) = preprune(J0,F0);
	  );
     if opts.Minimize then
       (J0,F0) = pruneParameterScheme(J0,F0);
     (J0,F0)
     )

beginDocumentation()
document { 
     Key => ParameterSchemes,
     Headline => "a Macaulay2 package for local equations of Hilbert and other parameter schemes",
     EM "ParameterSchemes", " is a package containing tools to create parameter schemes, 
     especially centered about a monomial ideal.",
     PARA{},
     "An example of using the functions in this package:",
     EXAMPLE lines ///
         R = ZZ/101[a..e,MonomialOrder=>Lex];
         I = ideal"ab,bc,cd,ad";
         L1 = smallerMonomials I
	 ///,
     PARA{"We will construct the family of all ideals having I as its lexicographic initial ideal."},
     EXAMPLE lines ///
         F0 = parameterFamily(I,L1,symbol t);
	 netList F0_*
	 J0 = parameterIdeal(I,F0);
         ///,
     "At this point F0 is the universal family, and J0 is a (very non-minimally generated) ideal that
     the parameters must satisfy in order for the family to be flat.",
     PARA{},
     "We can minimalize the family, and the ideal.",
     EXAMPLE lines ///
         (J,F) = pruneParameterScheme(J0,F0);
	 J
	 netList first entries gens F
	 ///,
     "Notice that J is zero.  This means that the base is an affine space (in this case affine 8-space).
     Now let's find a random fiber over the base:",
     EXAMPLE lines ///
         B = ring J
	 S = ring F
	 rand = map(R,S,(vars R) | random(R^1, R^(numgens B)))
	 L = rand F
         ///,
     "In some sense, L is the 'generic' ideal having I as its lexicographic initial ideal.  
     Let's investigate
     L further:",
     EXAMPLE lines ///
	 leadTerm L
         betti res L
     	 primaryDecomposition L	 
     ///,
     "Note that this strongly indicates that every ideal with I as its 
     lexicographic initial ideal is not prime."
     }

document {
     Key => {(groebnerScheme,Ideal),groebnerScheme},
     Headline => "find the family of all ideals having a given monomial ideal as initial ideal",
     Usage => "(J,F) = groebnerScheme I\n(J,F) = groebnerScheme(I, Minimize=>false)",
     Inputs => { "I" => "a monomial ideal in a polynomial ring R",
	  Minimize => "set to false if minimalization of the ideal and family is 
	  not desired, or is too compute intensive" },
     Outputs => {
	  "J" => Ideal => "the ideal defining the base space",
	  "F" => Ideal => "the family"
	  },
     "The ideal J is in a ring A = kk[t_1, t_2, ....].  The scheme defined by J
     is the Groebner scheme of (I,>), where > is the monomial order in the ring of I.
     The ideal F is the ideal of the family, in the ring: A[gens R] (more precisely: A (monoid R)",
     PARA{},
     "As an example, we compute the groebner scheme of the following ideal.  The resulting parameter space
     is affine 8-space, and so is smooth, rational and irreducible.",	  
     EXAMPLE lines ///
         R = ZZ/101[a..e];
	 I = ideal"ab,bc,cd,ad";
	 (J,F) = groebnerScheme I;
	 J
	 netList first entries gens F
	 ///,
     SourceCode => {(groebnerScheme,Ideal)},
     SeeAlso => {parameterIdeal, parameterFamily, pruneParameterScheme, smallerMonomials}
     }


end
restart
loadPackage "ParameterSchemes"
installPackage ParameterSchemes
viewHelp ParameterSchemes

R = ZZ/101[x_1..x_4, MonomialOrder =>{Lex => 2, GroupLex => 2}, Global => false]
M = (monoid R).Options.MonomialOrder
w = {2, 4}
inducedMonomialOrder(M,w)


-- Example: triangle, giving twisted cubic --
kk = ZZ/101
R = kk[a..d]
I = ideal"ab,bc,ca"
time (J,F) = groebnerScheme(I);
A = ring J; B = ring F
-- Since J is 0, let's see what a random such fiber looks like
phi = map(R,B,(vars R)|random(R^1, R^(numgens A)))
L = phi F
leadTerm L
decompose L

-- Hi Amelia, here is a good example:
-- Amelia Amelia Amelia Amelia Amelia Amelia Amelia Amelia
kk = ZZ/101
R = kk[a..f]
I = ideal"ab,bc,cd,de,ea,ac"
time (J,F) = groebnerScheme(I);
time (J,F) = groebnerScheme(I, Minimize=>false);
time minimalPresentation J;
A = ring J
B = ring F
Alocal = kk{gens A, MonomialSize => 8} 
Jlocal = sub(J,Alocal);
gbTrace=3
Jlocal = ideal gens gb Jlocal
J1 = trim(ideal(t_40) + Jlocal)
J2 = Jlocal : t_40
P = prune J2;
-- Now: I want a point on V(J), on this component...
-- This should be a function?
A1 = ring P
phi2
psi = map(kk,A1,random(kk^1, kk^(numgens A1)))
g = psi * phi2
g Jlocal
g = map(kk,A,g.matrix)
g J
g' = map(R,B,(vars R) | sub(g.matrix,R))
L = g' F
leadTerm L
res L

use ring J
J' = ideal apply(flatten entries gens J, f -> f // t_40);
minimalPresentation J'
-- There are two smooth components through the ideal I.
-- Amelia Amelia Amelia Amelia Amelia Amelia Amelia Amelia

-- Example: in the local case --
kk = ZZ/101
R = kk[a,b,c]
I = ideal"ab,ac,bc"
Z = syz gens I
A = apply(degrees source gens I, d -> matrix basis(d,comodule I))
B = apply(degrees source Z, d -> matrix basis(d, coker Z))
ngens = sum apply(A, a -> numgens source a)
nsyz = sum apply(B, b -> numgens source b)
S = kk[gens R, s_1..s_nsyz, t_1..t_ngens, MonomialSize=>8]
firstvar = numgens R+nsyz-1;
G' = matrix {apply(numgens I, i -> (
	  substitute(I_i,S) + sum apply(first entries A_i, m -> (
		    firstvar=firstvar+1;substitute(m,S) * S_firstvar))
	  ))}
firstvar = numgens R - 1;
Z' = matrix{apply(numgens source Z, i -> (
	  substitute(Z_{i},S) + sum apply(numgens source B_i, j -> (
		    firstvar=firstvar+1;S_firstvar ** substitute((B_i)_{j}, S)))
	  ))}
coefficients(G' * Z', Variables => {S_0 .. S_(numgens R-1)})
J = ideal flatten entries oo_1
J0 = (minimalPresentation J)_0
S = kk{gens S, MonomialSize=>8}
J0 = substitute(J0,S)
gbTrace=3
gens gb J0

Istd = standardMonomials I
F = parameterFamily(I,Istd,symbol t,Local=>true)
parameterIdeal(I,F)
f = F_0
g = F_1
h = F_2
c*f-b*g + t_4*a*f -a*g*t_1 - b*t_2*h
(b*c^2) % G
b*c^2 - c*h + a*t_7*g
--------------------------------


R = ZZ/101[a..e,MonomialOrder=>Lex]
I = ideal"ab,bc,cd,ad"
time (J,F) = groebnerScheme I

L1 = smallerMonomials I
F0 = parameterFamily(I,L1,symbol t)
J0 = parameterIdeal(I,F0)
(J,F) = pruneParameterScheme(J0,F0)
B = ring J
rand = map(R,ring F,(vars R) | random(R^1, R^(numgens B)))
L = rand F
primaryDecomposition L
primaryDecomposition I

R = ZZ/101[a..f]
I = ideal"ab,bc,cd,ad"
L1 = smallerMonomials I
F0 = parameterFamily(I,L1,symbol t)
J0 = parameterIdeal(I,F0)
(J,F) = pruneParameterScheme(J0,F0)
B = ring J
rand = map(R,ring F,(vars R) | random(R^1, R^(numgens B)))
L = rand F
primaryDecomposition L
primaryDecomposition I

R = ZZ/101[a..f]
I = ideal"ab,bc,cd,de,ea,ad"
(J,F) = groebnerScheme I
Alocal = kk{gens ring J, MonomialSize=>8}
Jlocal = sub(J,Alocal)
gbTrace=3
time gens gb Jlocal;

L1 = smallerMonomials I
F0 = parameterFamily(I,L1,symbol t)
J0 = parameterIdeal(I,F0);
time (J,F) = pruneParameterScheme(J0,F0);
B = ring J
rand = map(R,ring F,(vars R) | random(R^1, R^(numgens B)))
L = rand F
decompose J
intersect oo == J -- yes
primaryDecomposition L
primaryDecomposition I

R = ZZ/101[a..e]
I = ideal"ab,bc,cd,ade"
time (J,F) = groebnerScheme(I, Minimize=>false);
time (J,F) = groebnerScheme(I);
B = ring J
rand = map(R,ring F,(vars R) | random(R^1, R^(numgens B)))
L = rand F
decompose L
intersect oo == L -- 
primaryDecomposition L
primaryDecomposition I


L1 = smallerMonomials I
F0 = parameterFamily(I,L1,symbol t)
J0 = parameterIdeal(I,F0);
time minimalPresentation J0;
time (J,F) = pruneParameterScheme(J0,F0);


R = ZZ/101[a..f,MonomialOrder=>Lex]
I = ideal"ab,bc,cd,ad,de"
L1 = smallerMonomials I
F0 = parameterFamily(I,L1,symbol t)
S = ring F1
J0 = parameterIdeal(I,F0)
(J,F) = pruneParameterScheme(J0,F0)

time minimalPresentation J
substitute(J.cache.minimalPresentationMap F1,S)
see oo

T = ZZ/101[gens ring J]
JT = substitute(J,T)
gens gb JT;
debug ParameterSchemes
B = flatten entries ((gens F1) * syzM)
B_0 % G
B_1 % G
H = flatten entries gens G
netList H
netList B
B_1 - e*t_18*H_1 +c*t_13*H_0 - d*t_13*t_19*H_0 + e*t_14*H_0 - e*t_13*t_20*H_0
B_1 % G

matrix{{B_1}} % G

A = (ZZ/101){t_1..t_30}
S = A[gens R]
smallerMonomials flatten entries gens I
L2 = standardMonomials I
F1 = parameterFamily(I,L1,symbol t)
S = ring F1
J = parameterIdeal(I,F1)
T = ZZ/101[gens ring J]
J = substitute(J,T)
gens gb J;
K = ideal select(flatten entries oo, f -> first degree f > 1)
minimalPresentation K
decompose oo
gens gb K
decompose K
(gens substitute(F1,T)) % J
ideal((gens F1) % (trim parameterIdeal(I,F1)))
compress gens oo
R = ZZ/101[a..e, MonomialOrder=>Lex]
I = ideal"ab,bc,cd,ad"
smallerMonomials I
smallerMonomials flatten entries gens I
standardMonomials I

-- AFTER THIS: TESTS FOR ROUTINES THAT HAVE BEEN RENAMED...
end
restart
path = prepend("/Users/mike/Macaulay2/code/",path)
load "localhilb.m2"
R = ZZ/32003[a..c]
I = ideal(a^2,a*b^2)
(J,fam) = localEquations(I,symbol t)
transpose gens gb J


R = ZZ/32003[a..d]
(J,fam) = localEquations(ideal(a*b, a*c, b^3), symbol t)
gbTrace = 3
gens gb J;
G = forceGB gens fam
syzI
eqs


F = J_(0,0)
F % G
g1 = fam_0 -- lead term a*b
g2 = fam_1 -- lead term a*c
g3 = fam_2 -- lead term b^3
F
use ring F
F + t_1 * a * g2 + t_5* d * g2 - t_9*a*g1 -t_10*g3 -t_13*d*g1 - t_2*t_9*b*g1
gens gb J


getStandardMonomials(I,3)
getStds I


family = familyIdeal(I,getStds I,t)
U = ring family
time G = forceGB gens family;
time syzI = substitute(syz gens I, ring family);
eqs = ((gens family) * syzI)
eqs1 = eqs % G
eqs2 = compress eqs1
coefficients({a,b,c},eqs2)
J = ideal oo_1
gbTrace = 3
gens gb J;
transpose gens gb J
R = ZZ/32003[s,t_1..t_11,MonomialOrder=>Eliminate 1]
J = ideal(t_2-t_1*t_6+t_6^2,
     t_4-t_3*t_6+t_2*t_7-t_1*t_8+2*t_6*t_8-t_6^2*t_7,
     -t_10+t_7*t_8+t_6*t_9-t_6*t_7^2,
     t_5+t_4*t_7-t_3*t_8+t_8^2+t_2*t_9-t_1*t_10+t_6*t_10-t_6*t_7*t_8,
     -t_11+t_8*t_9-t_6*t_7*t_9,
     t_5*t_7+t_4*t_9-t_3*t_10+t_8*t_10-t_1*t_11+t_6*t_11-t_6*t_7*t_10,
     t_5*t_9-t_3*t_11+t_8*t_11-t_6*t_7*t_11)
J = homogenize(J,s)
gens gb J
transpose gens gb J
mingens ideal substitute(leadTerm gens gb J, {s=>1})

-- example: initial ideal is an edge ideal
restart
path = prepend("/Users/mike/Macaulay2/code/",path)
load "localhilb.m2"
R = ZZ/32003[a..f]
I = ideal"ab,bc,cd,de,ea"
(J,fam) = localEquations(I,symbol t)
transpose gens gb J

--- Starting a possible test list for inducedMonomialOrder


restart
loadPackage"ParameterSchemes"

kk = ZZ/32003
R1 = kk[x_1..x_5, MonomialOrder => {GRevLex => 3, Weights => {2,2}, Lex => 2}]
R2 = kk[x_1..x_5, MonomialOrder => {GroupLex => 2, GroupRevLex => 1, Weights => {2,2}, Lex => 2}, Global => false]
R3 = kk[x_1..x_5, MonomialOrder => {Weights => {1,2,3}, Lex => 3, Weights => {2,2}, Lex => 2}]
M1 = (monoid R1).Options.MonomialOrder
M2 = (monoid R2).Options.MonomialOrder
M3 = (monoid R3).Options.MonomialOrder
l1 = {0, 1, 4}
l2 = {1, 2, 3}
l3 = {1, 4}
inducedMonomialOrder(M1,l1)
inducedMonomialOrder(M1,l2)
inducedMonomialOrder(M1,l3)
inducedMonomialOrder(M2,l1)
inducedMonomialOrder(M2,l2)
inducedMonomialOrder(M2,l3)
inducedMonomialOrder(M3,l1)
inducedMonomialOrder(M3,l2)
inducedMonomialOrder(M3,l3)
inducedMonomialOrder(M3,{0,2,3})
inducedMonomialOrder(M3,{1,3,4})

-- Tests of new minimal presentation (of ring) code:
restart
loadPackage "ParameterSchemes"

TEST ///
A = ZZ/101[x,y]/(y-x^3-x^5-x^7)
B = minimalPresentation A
F = A.minimalPresentationMap
G = A.minimalPresentationMapInv
assert(G*F == map(A,A,gens A))
assert(F*G == map(B,B,gens B))
assert(ideal B == 0)
assert(numgens B == 1)
///

TEST ///
R = ZZ/101[x,y]
I = ideal(y-x^3-x^5-x^7)
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(numgens ring J == 1)
assert(J == 0)
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)
///


TEST ///
A = QQ[x,y]/(y-x^3-x^5-x^7)
B = minimalPresentation A
F = A.minimalPresentationMap
G = A.minimalPresentationMapInv
assert(G*F == map(A,A,gens A))
assert(F*G == map(B,B,gens B))
assert(ideal B == 0)
assert(numgens B == 1)
///

TEST ///
A = QQ[x,y]/(2*y-x^3-x^5-x^7)
B = minimalPresentation A
F = A.minimalPresentationMap
G = A.minimalPresentationMapInv
assert(G*F == map(A,A,gens A))
assert(F*G == map(B,B,gens B))
assert(ideal B == 0)
assert(numgens B == 1)
///

TEST ///
R = QQ[x,y]
I = ideal(2*y-x^3-x^5-x^7)
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(numgens ring J == 1)
assert(J == 0)
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)
///

TEST ///  -- FAILS
A = ZZ[x,y,z]/(2*y+z-x^3-x^5-x^7, z^2)
B = minimalPresentation A
F = A.minimalPresentationMap
G = A.minimalPresentationMapInv
assert(G*F == map(A,A,gens A))
assert(F*G == map(B,B,gens B))
assert(ideal B == 0)
assert(numgens B == 1)

TEST ///
R = ZZ[x,y,z]
I = ideal(2*y+z-x^3-x^5-x^7, z^2)
J = minimalPresentation I
assert(numgens ring J == 2)
use ring J
assert(J == ideal"x14+2x12+3x10+2x8-4x7y+x6-4x5y-4x3y+4y2")
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(numgens ring J == 2)
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)
///

A = QQ[a,b,c]/(a^2-3*b,a*c-c^4*b)
I = ideal 0_A
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(numgens ring J == 2)
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)

A = QQ[a,b,c]/(a^2-3*b^2,a^3-c^4*b)
I = ideal 0_A
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)

A = ZZ/101[a,b]/(a^2+b^2)
B = A[c,d]/(a*c+b*d-1)
C = B[e,f]/(e^2-b-1)
I = ideal 0_C
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)

I = ideal presentation (flattenRing C)#0
J = minimalPresentation I
F = I.cache.minimalPresentationMap
G = I.cache.minimalPresentationMapInv
assert(target F === ring I)
assert(source F === ring J)
assert(target G === ring J)
assert(source G === ring I)

-- Examples from parameter schemes
kk = ZZ/101
R = kk[a..f]
I = ideal"ab,bc,cd,de,ea,ac"
time (J,F) = groebnerScheme(I);

time (J,F) = groebnerScheme(I, Minimize=>false);
time J = minimalPresentation J;

kk = ZZ/101
R = kk[a..d]
I = ideal borel monomialIdeal"b2c"
time (J,F) = groebnerScheme(I, Minimize=>false);
time J = minimalPresentation J;
betti J
J = trim J
see J
R = ZZ/101[x,y]/(y-x^3-x^5-x^7)
I = ideal presentation R
minimalPresentation I

-- Test of parameterRing:
restart
loadPackage "ParameterSchemes"
kk = ZZ/32003
R = kk[a..d]
B = ideal{a*c, a*b, a^2, b^4*c, b^5}
L = smallerMonomials(B)
D = parameterRing(B,L,symbol t)
#D
S = kk[t_1..t_60, Degrees=>D]
debug Core
raw S
(matrix D) * transpose matrix{{0,-7,-9,-11}}

