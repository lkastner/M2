-- returns the matrix whose columns are a basis of the lineality space
-- this is always set by the constructor of F
linealitySpace Fan := F -> F.cache.computedLinealityBasis ??=
    makeRaysPrimitive mingens image F.cache.inputLinealityGenerators


-- returns the matrix whose columns are primitive generators of the 1-dim cones of the fan
rays Fan := Matrix => {} >> o -> F -> F.cache.rays ??= (
    -- TODO: could this be computed if honestMaxObjects is set?
    inputrays := try F.cache.inputRays else error "no input rays given for the fan";
    -- defined in Polyhedra/helpers.m2
    makeRaysUniqueAndPrimitive(inputrays, linealitySpace F)
)


-----------------------------------------------------------------------------
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


-----------------------------------------------------------------------------
-- testing equality
Fan == Fan := (F1, F2) -> (
    R1 := rays F1;
    R2 := rays F2;
    if numrows R1 != numrows R2
    or numcols R1 != numcols R2 then return false;
    M1 := maxCones F1;
    M2 := maxCones F2;
    if #M1 != #M2 then return false;
    L1 := linealitySpace F1;
    rayMap := rayCorrespondenceMap(R1, L1, R2);
    M1Mapped := apply(M1, C -> sort apply(C, l -> rayMap#l));
    M2Mapped := apply(M2, C -> sort C);
    sort M1Mapped == sort M2Mapped
)


-----------------------------------------------------------------------------
-- returns the dimension of F, i.e. maximum of the dimension of its cones
dim Fan := F -> F.cache.dim ??= (
    R := rays F;
    L := linealitySpace F;
    max apply(maxCones F, m -> rank(R_m | L))
)


-- returns the dimension of the ambient space of F
ambientDimension Fan := F -> F.cache.ambientDimension ??= (
    for key in keys rayProperties do (
	if F.cache#?key then return numRows F.cache#key);
    -- TODO: when can this ever happen?
    error "no property available to compute ambient dimension"
)


-----------------------------------------------------------------------------
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


-- whether or not F is smooth
isSmooth Fan := {} >> o -> F -> F.cache.isSmooth ??= (
    -- TODO: could it be better to run isSimplicial first?
    try if not F.cache.isSimplicial then return false;
    -- TODO: could this be faster?
    -- if F.cache.?computedMaxCones then return all(computedMaxCones F, isSmooth);
    R := rays F;
    L := transpose linealitySpace F;
    MC := maxCones F;
    MC = apply(MC, m -> R_m);
    -- defined in Polyhedra/helpers.m2
    all(MC, r -> spanSmoothCone(transpose r, L))
)


-- whether or not F is simplicial
isSimplicial Fan := F -> F.cache.isSimplicial ??= (
    try if F.cache.isSmooth      then return true;
    if F.cache.?computedMaxCones then return all(computedMaxCones F, isSimplicial);
    R := rays F;
    L := linealitySpace F;
    all(maxCones F,
	m -> (
	    testmat := R_m | L;
	    (numColumns testmat) == (rank testmat)
	)
    )
)


-- whether or not F is pointed
isPointed Fan := F -> F.cache.isPointed ??= all(computedMaxCones F, isPointed)


-- whether or not F is pure, i.e. all maximal cones have the same dimension
isPure Fan := F -> F.cache.isPure ??= (
    d := dim F;
    if F.cache.?computedMaxCones then
	return all(computedMaxCones F, C -> dim C == d);
    R := rays F;
    L := linealitySpace F;
    all(maxCones F, m -> rank(R_m | L) == d)
)


-- whether or not F is complete, i.e. its support is the entire space
isComplete Fan := F -> F.cache.isComplete ??= (
    n := dim F;
    if n != ambientDimension F then return false;
    -- TODO: simplify this
    symmDiff := (X,Y) -> (
	summand1 := select(X, x -> position(Y, y->y==x) === null);
	summand2 := select(Y, y -> position(X, x->y==x) === null);
	flatten {summand1, summand2}
    );
    Lfaces := {};
    CFsave := {};
    scan(computedMaxCones F,
	C -> (
	    if dim C == n then (
		R := rays C;
		L := linealitySpace C;
		CFacets := toList getProperty(C, facetsThroughRayData);
		CFacets = apply(CFacets, facet -> coneFromVData(R_facet, L));
		CFsave = flatten {CFsave, {CFacets}};
		Lfaces = symmDiff(Lfaces, CFacets);
	    )
	    else return false
	)
    );
    Lfaces == {}
)


-- whether or not F is polytopal
-- TODO: move isPolytopal from extended/not_refactored.m2
isPolytopal = method(TypicalValue => Boolean)
isPolytopal Fan := F -> getProperty(F, polytopal)


-- if the fan is polytopal, returns a polytope P such that F is the normal fan of P
polytope = method(TypicalValue => Polyhedron)
polytope Fan := F -> (
    -- isPolytopal computes and caches the polytope
    if isPolytopal F then F.cache.computedPolytope
    else error "expected a polytopal fan"
)


-- returns a hash table with entries d => { C | codim C == d }
-- TODO: rename the cache key
faces Fan := HashTable => F -> F.cache.computedFacesThroughRays ??= (
    dimF := dim F;
    raysF := rays F;
    L := linealitySpace F;
    MC := computedMaxCones F;
    -- a hash table with entries C => codim C
    allcones := new MutableHashTable;
    for C in MC do (
	dimC := dim C;
	raysC := rays C;
	facesC := faces C;
	-- defined in Polyhedra/helpers.m2
	rc := rayCorrespondenceMap(raysC, L, raysF);
	for i in keys facesC do (
	    codimInF := i + dimF - dimC;
	    for C' in facesC#i do (
		raysC' := sort apply(C', e -> rc#e);
		allcones#raysC' = codimInF);
	);
    );
    partition(C -> allcones#C, sort keys allcones)
)

-- returns a list of codim k cones of F
-- compare with computedMaxCones
-- TODO: should these be cached?
facesAsCones(ZZ, Fan) := List => (k, F) -> (
    R := rays F;
    L := linealitySpace F;
    apply(faces(k, F), m -> coneFromVData(R_m, L))
)

fVector Fan := F -> F.cache.fVector ??= apply(dim F + 1, d -> #faces(dim F - d, F))

-- returns a list of d-dimensional cones of a fan
cones = method()
cones(ZZ, Fan) := List => (d, F) -> faces(dim F - d, F)

-- UNEXPORTED METHOD
-- TODO: export this, currently only used by smoothSubfan
-- returns the smooth cones of a fan as a list of index lists
smoothCones = method()
smoothCones Fan := List => F -> F.cache.smoothCones ??= (
    allcones := faces F;
    R' := transpose rays F;
    L' := transpose linealitySpace F;
    flatten for i in keys allcones list select(allcones#i,
	-- defined in Polyhedra/helpers.m2
	m -> spanSmoothCone(R'^m, L'))
)
