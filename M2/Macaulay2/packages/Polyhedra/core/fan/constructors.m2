-- Defining the new type Fan
Fan = new Type of PolyhedralObject
Fan.synonym = "polyhedral fan"
Fan.GlobalAssignHook = globalAssignFunction
Fan.GlobalReleaseHook = globalReleaseFunction
compute#Fan = new MutableHashTable;

-- properties which encode rays of a cone or fan
-- TODO: combine with the one in core/cone/constructors.m2
rayProperties = set { rays, inputRays, inputLinealityGenerators, computedLinealityBasis }

-- TODO: complete this list from symbols.m2
-- fanProperties = set { inputCones, ambientDimension,
--     computedComplete, computedFVector,
--     computedLinealityBasis, computedLinealityGenerators, computedFacesThroughRays,
--     polytopal, smooth, simplicial, pure, facetsThroughRayData, isWellDefined }

-----------------------------------------------------------------------------
-- Low level constructor of a fan

new Fan from HashTable := (Fan, data) -> constructTypeFromHash(Fan,
    applyPairs(data, (name, property) -> name =>
	if rayProperties#?name then makeRaysPrimitive property else property
	-- TODO: can we also list all possible properties of a fan for sanitization?
	--if fanProperties#?name then                   property else
	--error("encountered unexpected property ", toString name, " when constructing a fan")
    )
)

-----------------------------------------------------------------------------
-- Equality testing

Fan == Fan := (F1, F2) -> (
   r1 := rays F1;
   l1 := linealitySpace F1;
   r2 := rays F2;
   if numRows r1 != numRows r2 then return false;
   if numColumns r1 != numColumns r2 then return false;
   m1 := maxCones F1;
   m2 := maxCones F2;
   if #m1 != #m2 then return false;
   rayMap := rayCorrespondenceMap(r1, l1, r2);
   m1Mapped := apply(m1,
      maxCone -> (
         sort apply(maxCone, l -> rayMap#l)
      )
   );
   m2Mapped := apply(m2, maxCone -> sort maxCone);
   (sort m1Mapped) == (sort m2Mapped)
)

-----------------------------------------------------------------------------
-- Main constructors of a fan

-- PURPOSE : Building the Fan 'F'
--   INPUT : 'L',  a list of cones and fans in the same ambient space
--  OUTPUT : The fan of all Cones in 'L' and all Cones in of the fans in 'L' and all their faces
fan = method(TypicalValue => Fan)
fan(Matrix, Matrix, List)     :=
fan(Matrix, Matrix, Sequence) := (inputrays, linealityGens, inputcones) -> (
    inputcones = toList inputcones;
    lineality := makeRaysPrimitive(mingens image linealityGens);
    if numcols inputrays < max flatten inputcones then error "the number of indices exceeds the number of rays";
    if numrows inputrays != numrows linealityGens then error "rays and lineality must have same ambient dimension";
    new Fan from hashTable {
	inputRays              => inputrays,
	inputCones             => inputcones,
	computedLinealityBasis => lineality,
    }
)

fan(Matrix, List)     :=
fan(Matrix, Sequence) := (inputrays, inputcones) -> (
    linealityGens := map(target inputrays, (ring inputrays)^0, 0);
    fan(inputrays, linealityGens, inputcones)
)

-- returns the fan given by C and all of its faces
fan Cone := C -> (
    inputrays := rays C;
    inputcone := for i from 0 to numcols inputrays - 1 list i;
    F := fan(inputrays, linealitySpace C, { inputcone });
    F.cache.honestMaxObjects = new HashTable from { inputcone => C };
    F
)

importFrom_Core {"concatCols"}
cols = m -> apply(numcols m, j -> m_{j})

fan List := inputCones -> (
    A := apply(inputCones, C ->
	if instance(C, Cone) then cols rays C else
	if instance(C, Matrix) then cols C else
	error "Fan constructor expected a list of cones or matrices");
    B := unique flatten A;
    H := hashTable apply(toList pairs B, reverse);
    rayList := concatCols B;
    maxList := apply(A, C -> apply(C, ray -> H#ray));
    fan(rayList, -* linealityGens, *- maxList)
)

-----------------------------------------------------------------------------
-- Other methods returning a fan
-- Also see stellarSubdivision and ccRefinement in extended/fan/methods.m2

linearTransform = method()
linearTransform(Fan, Matrix) := Fan => (F, A) -> (
    newRays := A * rays F;
    newLineality := linealitySpace F;
    if rank newLineality != rank(newLineality | gens kernel A)
    then printerr "Warning: output fan may not be well-defined. Check with 'isWellDefined'";
    newLineality = mingens image(A * newLineality);
    goodNewRays := makeRaysUniqueAndPrimitive(newRays, newLineality);
    new Fan from hashTable {
	inputRays              => goodNewRays,
	maxCones               => maxCones F,
	computedLinealityBasis => newLineality
    }
)

-- PURPOSE : Computing the subfan of all smooth cones of the Fan
--   INPUT : 'F',  a Fan
--  OUTPUT : The Fan of smooth cones
smoothSubfan = method()
smoothSubfan Fan := Fan => F -> (
    new Fan from hashTable {
	inputRays              => rays F,
	inputCones             => smoothCones F,
	computedLinealityBasis => linealitySpace F
    }
)


-- PURPOSE : Computing the 'n'-skeleton of a fan
--   INPUT : (n,F),  where 'n' is a positive integer and
--                   'F' is a Fan
--  OUTPUT : the Fan consisting of the 'n' dimensional cones in 'F'
skeleton = method()
skeleton(ZZ, Fan) := Fan => (n, F) -> (
    -- Checking for input errors
    if n < 0 or dim F < n then error "expected an integer be between 0 and the dimension of the fan";
    new Fan from hashTable {
	inputRays              => rays F,
	inputCones             => cones(n, F),
	computedLinealityBasis => linealitySpace F
    }
)

-- TODO: move to gfanInterface?
fanFromGfan = method()
fanFromGfan List := Fan => gfanOutput -> (
-- 0 rays -> Matrix
-- 1 lineality -> Matrix
-- 2 cones -> List<List>
-- 3 dimension -> ZZ
-- 4 pure -> bool
-- 5 simplicial -> bool
-- 6 fVector -> List
   numberOfGfanOutputs := 7;
   if #gfanOutput != numberOfGfanOutputs then
      error("fanFromGfan was given a list with " | toString(#gfanOutput)
         | " inputs and " | toString(numberOfGfanOutputs) | " are required.");
   R := gfanOutput#0;
   L := gfanOutput#1;

   -- Perform some basic sanity checks on the fan. If the fan is empty (i.e.
   -- there are no rays and no lineality space), then all of the other values
   -- need to agree with that (there cannot be any cones, the dimension must
   -- be zero, and the f-vector must be empty).
   if ((numColumns R == 0) and (numColumns L == 0))
   and ((#(gfanOutput#2) != 0) or (gfanOutput#3 != 0) or (#(gfanOutput#6) != 0))
   then error("Inconsistent input into fanFromGfan");

   if (numColumns R == 0) then R = map(ZZ^(numRows L), ZZ^0, 0);
   if (numColumns L == 0) then L = map(ZZ^(numRows R), ZZ^0, 0);
    result := fan(R, L,              gfanOutput#2);
    result.cache.dim               = gfanOutput#3;
    result.cache.pure              = gfanOutput#4;
    result.cache.simplicial        = gfanOutput#5;
    result.cache.computedFVector   = gfanOutput#6;
    result
)

addCone = method()
addCone(Fan, Cone) := Fan => (F, C) -> (
   if ambDim F != ambDim C then error("Fan and Cone must live in same ambient space.");
   linF := linealitySpace F;
   linC := linealitySpace C;
   if image linF != image linC then error("Cannot add cone with different lineality space.");
   joinedRays := makeRaysUniqueAndPrimitive(rays F | rays C, linF);
   mc := maxCones F;
   map := rayCorrespondenceMap(rays F, linF, joinedRays);
   mc = apply(mc, c-> apply(c, e->map#e));
   map = rayCorrespondenceMap(rays C, linF, joinedRays);
   newCone := toList apply(numColumns rays C, i -> map#i);
   mc = append(mc, newCone);
   new Fan from hashTable {
      inputRays              => joinedRays,
      inputCones             => mc,
      ambientDimension       => ambDim F,
      computedLinealityBasis => linF,
   }
)

-- TODO: deprecate these
addCone(Cone, Fan) := Fan => (C, F) -> addCone(F, C)
addCone(List, Fan) := Fan => (L, F) -> (
   if not all(L, l->instance(l, Cone)) then error("List does not contain cones");
   result := F;
   for cone in L do (
      result = addCone(result, cone)
   );
   result
)
