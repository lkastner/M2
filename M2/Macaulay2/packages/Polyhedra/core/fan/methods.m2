-- returns the matrix whose columns are a basis of the lineality space
-- this is always set by the constructor of F
linealitySpace Fan := F -> F.cache.computedLinealityBasis


-- returns the matrix whose columns are primitive generators of the 1-dim cones of the fan
rays Fan := Matrix => {} >> o -> F -> F.cache.rays ??= (
    -- TODO: could this be computed if honestMaxObjects is set?
    inputrays := try F.cache.inputRays else error "no input rays given for the fan";
    -- defined in Polyhedra/helpers.m2
    makeRaysUniqueAndPrimitive(inputrays, linealitySpace F)
)


-- returns a list of index lists for the generating Cones of the fan
maxCones = method(TypicalValue => List)
maxCones Fan := F -> F.cache.maxCones ??= (
    -- TODO: could this be computed if honestMaxObjects is set?
    cones := try F.cache.inputCones else error "no input cones given for the fan";
    if F.cache.?inputRays then (
	raysF := rays F;
	inputRaysF := F.cache.inputRays;
	linealityF := linealitySpace F;
	-- defined in Polyhedra/helpers.m2
	rc := rayCorrespondenceMap(inputRaysF, linealityF, raysF);
	cones = apply(cones, C -> select(sort apply(C, e -> rc#e), e -> e != -1)));
    cones = unique apply(cones, sort);
    result := new MutableList;
    -- TODO: simplify this
    for C in cones do (
	test := all(cones,
	    c -> (
		n := #((set c) * (set C));
		if n == #C then (
		    C == c
		) else (
		    true
		)
	    )
	);
	if test then result##result = C;
    );
    toList result
)

-- TODO: do we need this for general purposes, or can it be deprecated?
maxObjects Fan := maxCones

-- UNEXPORTED METHOD
-- returns a list of maximal cones as honest cones
computedMaxCones = method()
computedMaxCones Fan := List => F -> F.cache.computedMaxCones ??= (
    R := rays F;
    L := linealitySpace F;
    apply(maxCones F, m -> coneFromVData(R_m, L)))

-- UNEXPORTED METHOD
-- returns a hash table of maximal cones as honest cones
honestMaxObjects = method()
honestMaxObjects Fan := HashTable => F -> F.cache.honestMaxObjects ??=
    new HashTable from apply(maxCones F, computedMaxCones F, identity)


-- whether or not F is well-defined
isWellDefined Fan := F -> F.cache.isWellDefined ??= (
    cones := honestMaxObjects F;
    n := #cones;
    for i from 0 to n-1 do (
	ki := (keys cones)#i;
	Ci := cones#ki;
	if(#ki != numColumns rays Ci) then(
	    if debugLevel > 0 then << "The cone " << ki << " has redundant rays." << endl;
	    return false;
	);
	for j from i to n-1 do (
	    kj := (keys cones)#j;
	    Cj := cones#kj;
	    -- defined in Polyhedra/extended/commonFace.m2
	    if not commonFace(Ci, Cj) then (
		if debugLevel > 0 then << "The cones " << ki << " and " << kj << " do not intersect in a common face." << endl;
		return false
	    )
	)
    );
    return true
)


--   INPUT : 'F',  a Fan
--  OUTPUT : 'true' or 'false'
isPointed Fan := F -> all(computedMaxCones F, isPointed)


-- PURPOSE : Checks if the input is smooth
--   INPUT : 'F'  a Fan
--  OUTPUT : 'true' or 'false'
isSmooth Fan := {} >> o -> F -> getProperty(F, smooth)


isPolytopal = method(TypicalValue => Boolean)
isPolytopal Fan := F -> getProperty(F, polytopal)

-- PURPOSE : Giving the k dimensional Cones of the Fan
--   INPUT : (k,F)  where 'k' is a positive integer and F is a Fan 
--  OUTPUT : a List of Cones
cones = method(TypicalValue => List)
cones(ZZ,Fan) := (k,F) -> (
   d := dim F;
   faces := getProperty(F, computedFacesThroughRays);
   faces#(d-k)
)


-- PURPOSE : Returning a polytope of which the fan is the normal if the fan is polytopal
--   INPUT : 'F',  a Fan
--  OUTPUT : A Polytope of which 'F' is the normal fan
polytope = method(TypicalValue => Polyhedron)
polytope Fan := F -> getProperty(F, computedPolytope)



isPure Fan := F -> getProperty(F, pure)
isComplete Fan := F -> getProperty(F, computedComplete)

objectsOfDim(ZZ, Fan) := (k,F) -> (
	-- Checking for input errors
	if k < 0 or dim F < k then error("k must be between 0 and the dimension of the polyhedral object family.");
	L := select(computedMaxCones F, C -> dim C >= k);
	-- Collecting the 'k'-dim faces of all generating cones of dimension greater than 'k'
	unique flatten apply(L, C -> faces(dim(C)-k,C)))


facesAsCones(ZZ, Fan) := (d, F) -> (
   raysF := rays F;
   linF := linealitySpace F;
   result := faces(d, F);
   apply(result, f -> coneFromVData(raysF_f, linF))
)
